#pragma once

#include <vector>
#include <memory>
#include <Eigen/Core>
#include <open3d/geometry/TriangleMesh.h>
#include <open3d/geometry/PointCloud.h>
#include <unordered_map>

namespace spalt {

// Simple mesh structure for integration with Surface class
struct SimpleMesh {
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> triangles;
};

// Simple point cloud structure for integration with Surface class  
struct SimplePointCloud {
    std::vector<Eigen::Vector3d> points;
    std::vector<Eigen::Vector3d> normals;
};

/**
 * @brief Parameters for Fibonacci-based surface generation
 * 
 * Parameters controlling the Fibonacci sphere sampling approach for molecular surface
 * generation. This method uses Fibonacci spiral sampling on spheres around atoms
 * combined with Open3D's ball pivoting algorithm for triangulation.
 */
struct FibonacciSurfaceParams {
    double probe_radius = 1.;           ///< Probe radius in Angstroms for surface generation (default: 1.0)
    int num_samples_per_atom = 100;      ///< Base number of sample points per atom for Fibonacci sphere sampling (default: 100, max: 100)
    double cutoff = 1.4;                 ///< Cutoff distance for surface radius scaling and point filtering (default: 1.4)
    std::vector<double> ball_radii = {0.6, 1.0, 1.4}; ///< Multiple ball radii for Open3D ball pivoting algorithm triangulation (default: {0.6, 1.0, 1.4})
    bool estimate_normals = true;        ///< Whether to estimate vertex normals during surface generation (default: true)
    double carbon_radius = 1.45515;          ///< Reference carbon radius in Angstroms for atomic radius scaling (default: 1.45515)
    
    // Validation
    bool isValid() const {
        return probe_radius > 0.0 && 
               num_samples_per_atom > 0 && 
               num_samples_per_atom <= 100 &&
               cutoff > 0.0 &&
               !ball_radii.empty() &&
               carbon_radius > 0.0;
    }
};

/**
 * @brief Fibonacci-based molecular surface generator
 * 
 * This class implements an alternative surface generation method that:
 * 1. Samples points on spheres around each atom using Fibonacci spiral
 * 2. Filters points to create solvent-accessible surface
 * 3. Uses Open3D's ball pivoting algorithm to create triangular mesh
 */
class FibonacciSurfaceGenerator {
public:
    /**
     * @brief Constructor
     * @param params Surface generation parameters (default: FibonacciSurfaceParams with default values)
     */
    explicit FibonacciSurfaceGenerator(const FibonacciSurfaceParams& params = FibonacciSurfaceParams());
    
    /**
     * @brief Generate surface mesh and point cloud for a molecule
     * @param centers Atom centers matrix (N x 3) where N is the number of atoms
     * @param radii Atom radii vector (N) containing Van der Waals radii for each atom
     * @return Pair of (mesh, point_cloud) where mesh is triangulated surface and point_cloud contains sampled points
     */
    std::pair<std::shared_ptr<SimpleMesh>, std::shared_ptr<SimplePointCloud>>
    generateSurface(const Eigen::MatrixXd& centers, 
                   const Eigen::VectorXd& radii) const;

    /**
     * @brief Sample points from molecular surface using Fibonacci approach
     * @param centers Atom centers matrix (N x 3) where N is the number of atoms
     * @param radii Atom radii vector (N) containing Van der Waals radii for each atom
     * @return Sampled surface points matrix (M x 3) where M is the number of surface points
     */
    Eigen::MatrixXd sampleMolecularSurface(const Eigen::MatrixXd& centers, 
                                          const Eigen::VectorXd& radii) const;

    /**
     * @brief Get the current parameters
     * @return Current FibonacciSurfaceParams configuration
     */
    const FibonacciSurfaceParams& getParams() const { return params_; }

    /**
     * @brief Update parameters
     * @param params New FibonacciSurfaceParams configuration (will be validated)
     */
    void setParams(const FibonacciSurfaceParams& params);

    /**
     * @brief Generate points on unit sphere using Fibonacci spiral
     * @param num_samples Number of points to generate on the unit sphere (must be > 0)
     * @return Points on unit sphere matrix (num_samples x 3) with coordinates in [-1, 1]
     */
    Eigen::MatrixXd generateFibonacciSphere(int num_samples) const;

private:
    

    FibonacciSurfaceParams params_;
    
    // Cache for sphere samples to avoid regeneration
    mutable std::unordered_map<int, Eigen::MatrixXd> sphere_cache_;

    // Stores the originating atom index for each sampled surface point (updated per sample)
    mutable std::vector<int> last_sample_atom_indices_;
};

} // namespace spalt
