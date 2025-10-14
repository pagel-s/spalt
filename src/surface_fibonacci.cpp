#include "surface_fibonacci.h"
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <open3d/geometry/PointCloud.h>
#include <open3d/geometry/TriangleMesh.h>

#include <Eigen/Dense>

namespace spalt {

FibonacciSurfaceGenerator::FibonacciSurfaceGenerator(const FibonacciSurfaceParams& params)
    : params_(params) {
    if (!params_.isValid()) {
        std::cerr << "Warning: Invalid FibonacciSurfaceParams provided, using defaults" << std::endl;
        params_ = FibonacciSurfaceParams();
    }
}

std::pair<std::shared_ptr<SimpleMesh>, std::shared_ptr<SimplePointCloud>>
FibonacciSurfaceGenerator::generateSurface(const Eigen::MatrixXd& centers, 
                                         const Eigen::VectorXd& radii) const {
    
    // Sample points from molecular surface
    Eigen::MatrixXd surface_points = sampleMolecularSurface(centers, radii);
    
    // Create simple point cloud
    auto point_cloud = std::make_shared<SimplePointCloud>();
    for (int i = 0; i < surface_points.rows(); ++i) {
        point_cloud->points.push_back(surface_points.row(i).transpose());
        point_cloud->normals.push_back(Eigen::Vector3d::Zero()); // Will be estimated if needed
    }

    // Create simple mesh
    auto mesh = std::make_shared<SimpleMesh>();
    
    if (surface_points.rows() == 0) {
        return std::make_pair(mesh, point_cloud);
    }

    // Convert to Open3D point cloud for mesh generation
    auto open3d_pcd = std::make_shared<open3d::geometry::PointCloud>();
    for (int i = 0; i < surface_points.rows(); ++i) {
        open3d_pcd->points_.push_back(surface_points.row(i).transpose());
    }
    
    // Estimate normals for mesh generation
    if (params_.estimate_normals) {
        open3d_pcd->EstimateNormals();
    }
    
    // Generate mesh using ball pivoting
    std::shared_ptr<open3d::geometry::TriangleMesh> open3d_mesh = nullptr;
    std::vector<double> ball_radii = params_.ball_radii;
    if (ball_radii.empty()) {
        ball_radii.push_back(params_.probe_radius);
    }
    
    try {
        open3d_mesh = open3d::geometry::TriangleMesh::CreateFromPointCloudBallPivoting(*open3d_pcd, ball_radii);
    } catch (const std::exception& e) {
        // If ball pivoting fails, create a simple vertex-only mesh
        for (int i = 0; i < surface_points.rows(); ++i) {
            mesh->vertices.push_back(surface_points.row(i).transpose());
        }
        return std::make_pair(mesh, point_cloud);
    }
    
    // Convert Open3D mesh to SimpleMesh
    for (size_t i = 0; i < open3d_mesh->vertices_.size(); ++i) {
        mesh->vertices.push_back(open3d_mesh->vertices_[i]);
    }
    
    for (size_t i = 0; i < open3d_mesh->triangles_.size(); ++i) {
        mesh->triangles.push_back(open3d_mesh->triangles_[i]);
    }
    
    return std::make_pair(mesh, point_cloud);
}

Eigen::MatrixXd FibonacciSurfaceGenerator::sampleMolecularSurface(const Eigen::MatrixXd& centers, 
                                                                const Eigen::VectorXd& radii) const {
    if (centers.rows() != radii.rows()) {
        throw std::invalid_argument("Number of centers must match number of radii");
    }
    
    int num_atoms = centers.rows();
    if (num_atoms == 0) {
        return Eigen::MatrixXd(0, 3);
    }
    
    // Calculate surface radii (vdW + probe)
    Eigen::VectorXd surface_radii = radii.array() + params_.probe_radius;
    
    // Calculate number of samples per atom (scaled by relative radius, squared)
    Eigen::VectorXd radius_ratio = (radii.array() / params_.carbon_radius).matrix();
    Eigen::VectorXi num_samples_per_atom = (params_.num_samples_per_atom * radius_ratio.array().square()).cast<int>();
    
    // Ensure minimum of 1 sample per atom
    num_samples_per_atom = num_samples_per_atom.cwiseMax(1);
    
    // Calculate total number of points
    int total_points = num_samples_per_atom.sum();
    if (total_points == 0) {
        return Eigen::MatrixXd(0, 3);
    }
    
    Eigen::MatrixXd all_sphere_points(total_points, 3);
    last_sample_atom_indices_.assign(total_points, -1);

    int point_idx = 0;

    for (int atom_idx = 0; atom_idx < num_atoms; ++atom_idx) {
        const int num_samples = num_samples_per_atom(atom_idx);
        if (num_samples <= 0) {
            continue;
        }

        Eigen::MatrixXd sphere_points = generateFibonacciSphere(num_samples);

        // Optional random rotation: disabled for determinism; uncomment if needed for randomness
        // std::mt19937 rng(atom_idx);
        // std::uniform_real_distribution<double> dist(0.0, 1.0);
        // const double theta = dist(rng) * 2.0 * M_PI;
        // const double phi = std::acos(2.0 * dist(rng) - 1.0);
        // const double z = dist(rng) * 2.0 * M_PI;
        // Eigen::AngleAxisd rot1(theta, Eigen::Vector3d::UnitZ());
        // Eigen::AngleAxisd rot2(phi, Eigen::Vector3d::UnitY());
        // Eigen::AngleAxisd rot3(z, Eigen::Vector3d::UnitZ());
        // Eigen::Matrix3d rotation = (rot1 * rot2 * rot3).toRotationMatrix();
        // sphere_points = sphere_points * rotation.transpose();

        sphere_points *= surface_radii(atom_idx);
        sphere_points.rowwise() += centers.row(atom_idx);

        all_sphere_points.block(point_idx, 0, num_samples, 3) = sphere_points;
        std::fill_n(last_sample_atom_indices_.begin() + point_idx, num_samples, atom_idx);
        point_idx += num_samples;
    }

    // Filter points to create solvent-accessible surface
    std::vector<int> valid_indices;
    valid_indices.reserve(total_points);

    for (int point_i = 0; point_i < total_points; ++point_i) {
        bool is_exposed = true;
        for (int atom_j = 0; atom_j < num_atoms; ++atom_j) {
            const double dist = (all_sphere_points.row(point_i) - centers.row(atom_j)).norm();
            const double cutoff = radii(atom_j) + params_.probe_radius - 0.01;
            if (dist < cutoff) {
                is_exposed = false;
                break;
            }
        }
        if (is_exposed) {
            valid_indices.push_back(point_i);
        }
    }

    if (valid_indices.empty()) {
        last_sample_atom_indices_.clear();
        return Eigen::MatrixXd(0, 3);
    }

    Eigen::MatrixXd filtered_points(valid_indices.size(), 3);
    for (size_t idx = 0; idx < valid_indices.size(); ++idx) {
        filtered_points.row(idx) = all_sphere_points.row(valid_indices[idx]);
    }

    // Shrink atom index map to only exposed points
    std::vector<int> exposed_atom_indices(valid_indices.size());
    for (size_t idx = 0; idx < valid_indices.size(); ++idx) {
        exposed_atom_indices[idx] = last_sample_atom_indices_[valid_indices[idx]];
    }
    last_sample_atom_indices_.swap(exposed_atom_indices);

    return filtered_points;
}

Eigen::MatrixXd FibonacciSurfaceGenerator::generateFibonacciSphere(int num_samples) const {
    // Check cache first
    auto it = sphere_cache_.find(num_samples);
    if (it != sphere_cache_.end()) {
        return it->second;
    }
    
    if (num_samples <= 0) {
        return Eigen::MatrixXd(0, 3);
    }
    
    // Generate Fibonacci spiral points on unit sphere
    double offset = 2.0 / num_samples;
    double increment = M_PI * (3.0 - std::sqrt(5.0));
    
    Eigen::MatrixXd points(num_samples, 3);
    
    for (int i = 0; i < num_samples; ++i) {
        double y = (i * offset) - 1.0 + (offset / 2.0);
        double r = std::sqrt(1.0 - y * y);
        double phi = std::fmod((i + 1), num_samples) * increment;
        
        points(i, 0) = std::cos(phi) * r;
        points(i, 1) = y;
        points(i, 2) = std::sin(phi) * r;
    }
    
    // Cache the result
    sphere_cache_[num_samples] = points;
    
    return points;
}


void FibonacciSurfaceGenerator::setParams(const FibonacciSurfaceParams& params) {
    if (params.isValid()) {
        params_ = params;
        // Clear cache when parameters change
        sphere_cache_.clear();
    } else {
        open3d::utility::LogWarning("Invalid FibonacciSurfaceParams provided, keeping current parameters");
    }
}

} // namespace spalt
