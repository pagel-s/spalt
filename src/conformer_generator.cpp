#include "conformer_generator.h"
#include "molecule.h" // For MMFF function declarations
#include <GraphMol/Conformer.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

// Simple K-means implementation for conformer clustering
class SimpleKMeans {
  public:
    SimpleKMeans(int k, int max_iter = 100) : k_(k), max_iter_(max_iter) {
    }

    std::vector<int> fit(const std::vector<std::vector<double>>& data) {
        int n = data.size();
        if (n == 0 || k_ > n) {
            std::vector<int> labels(n, 0);
            return labels;
        }

        // Initialize centroids randomly
        std::vector<std::vector<double>> centroids(k_);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, n - 1);

        std::vector<int> used(n, 0);
        for (int i = 0; i < k_; ++i) {
            int idx;
            do {
                idx = dis(gen);
            } while (used[idx]);
            used[idx] = 1;
            centroids[i] = data[idx];
        }

        std::vector<int> labels(n, 0);

        for (int iter = 0; iter < max_iter_; ++iter) {
            // Assign points to clusters
            for (int i = 0; i < n; ++i) {
                double min_dist = std::numeric_limits<double>::max();
                int best_cluster = 0;

                for (int j = 0; j < k_; ++j) {
                    double dist = euclideanDistance(data[i], centroids[j]);
                    if (dist < min_dist) {
                        min_dist = dist;
                        best_cluster = j;
                    }
                }
                labels[i] = best_cluster;
            }

            // Update centroids
            std::vector<std::vector<double>> new_centroids(k_);
            std::vector<int> counts(k_, 0);

            for (int i = 0; i < n; ++i) {
                int cluster = labels[i];
                if (counts[cluster] == 0) {
                    new_centroids[cluster] = std::vector<double>(data[i].size(), 0.0);
                }

                for (size_t j = 0; j < data[i].size(); ++j) {
                    new_centroids[cluster][j] += data[i][j];
                }
                counts[cluster]++;
            }

            // Normalize centroids
            bool converged = true;
            for (int i = 0; i < k_; ++i) {
                if (counts[i] > 0) {
                    for (size_t j = 0; j < new_centroids[i].size(); ++j) {
                        new_centroids[i][j] /= counts[i];
                    }

                    // Check convergence
                    if (euclideanDistance(centroids[i], new_centroids[i]) > 1e-6) {
                        converged = false;
                    }
                }
            }

            centroids = new_centroids;

            if (converged)
                break;
        }

        return labels;
    }

  private:
    int k_;
    int max_iter_;

    double euclideanDistance(const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() != b.size())
            return std::numeric_limits<double>::max();

        double sum = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return std::sqrt(sum);
    }
};

ConformerGenerator::ConformerGenerator(const ConformerParams& params) : params_(params) {
    params_.normalizeWeights();
}

std::unordered_map<int, int> ConformerGenerator::generateConformers(
    std::shared_ptr<RDKit::ROMol> mol, const ConformerParams* params) {
    const ConformerParams& p = params ? *params : params_;

    // Add hydrogens if needed
    addHydrogens(mol);

    // Get conformer counts for each method
    std::vector<int> counts = p.getConformerCounts();
    int optimized_count = counts[0];
    int etkdg_count = counts[1];
    int uniform_count = counts[2];

    bool optimization_succeeded = true;

    // Try optimized conformers first
    if (optimized_count > 0 && p.enable_optimization) {
        try {
            sampleOptimizedConformers(mol, optimized_count, p);
        } catch (const std::exception& e) {
            std::cerr << "WARNING: Optimization failed for molecule, falling back to ETKDG/Uniform "
                         "sampling: "
                      << e.what() << std::endl;
            optimization_succeeded = false;
        }
    }

    // If optimization failed, redistribute weights
    ConformerParams adjusted_params = p;
    if (!optimization_succeeded) {
        adjusted_params.applyFallbackWeights();
        std::vector<int> fallback_counts = adjusted_params.getConformerCounts();
        etkdg_count = fallback_counts[1];
        uniform_count = fallback_counts[2];
    }

    // Generate ETKDG conformers
    if (etkdg_count > 0) {
        sampleETKDGConformers(mol, etkdg_count, adjusted_params);
    }

    // Generate uniform conformers
    if (uniform_count > 0) {
        sampleUniformConformers(mol, uniform_count, adjusted_params);
    }

    // Align conformers if requested
    if (p.align_conformers) {
        alignConformers(mol);
    }

    // Cluster and select representative conformers
    return clusterConformers(mol, p.num_clusters);
}

std::pair<std::shared_ptr<RDKit::ROMol>, std::unordered_map<int, int>>
ConformerGenerator::generateConformersFromSmiles(const std::string& smiles,
                                                 const ConformerParams* params) {
    // Parse SMILES
    std::shared_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) {
        // Try removing any potential brackets
        std::string clean_smiles = smiles;
        if (smiles.length() > 2 && smiles[0] == '"' && smiles.back() == '"') {
            clean_smiles = smiles.substr(1, smiles.length() - 2);
        }
        mol.reset(RDKit::SmilesToMol(clean_smiles));
    }

    if (!mol) {
        throw std::runtime_error("Could not parse SMILES string: " + smiles);
    }

    auto cluster_map = generateConformers(mol, params);
    return std::make_pair(mol, cluster_map);
}

void ConformerGenerator::sampleOptimizedConformers(std::shared_ptr<RDKit::ROMol> mol, int num_confs,
                                                   const ConformerParams& params) {
    if (num_confs <= 0)
        return;

    // Generate conformers with ETKDG first
    RDKit::DGeomHelpers::EmbedParameters embed_params;
    embed_params.useRandomCoords = false;  // Use ETKDG
    embed_params.randomSeed = params.random_seed;
    embed_params.enforceChirality = params.enforce_chirality;
    embed_params.numThreads = params.num_threads;
    embed_params.clearConfs = params.clear_existing;

    RDKit::INT_VECT conf_ids =
        RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, num_confs, embed_params);
    int num_confs_generated = conf_ids.size();

    if (num_confs_generated == 0) {
        throw std::runtime_error("Failed to generate conformers");
    }

    // Optimize conformers with MMFF
    std::vector<std::pair<int, double>> optimization_results;
    try {
        RDKit::MMFF::MMFFOptimizeMoleculeConfs(*mol, optimization_results, params.num_threads);

        // Check for convergence issues
        bool any_not_converged = false;
        for (const auto& result : optimization_results) {
            if (result.first != 0) {  // 0 means converged
                any_not_converged = true;
                break;
            }
        }

        if (any_not_converged) {
            std::cerr << "WARNING: Not all conformers converged during MMFF optimization"
                      << std::endl;
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("MMFF optimization failed: " + std::string(e.what()));
    }
}

void ConformerGenerator::sampleETKDGConformers(std::shared_ptr<RDKit::ROMol> mol, int num_confs,
                                               const ConformerParams& params) {
    if (num_confs <= 0)
        return;

    RDKit::DGeomHelpers::EmbedParameters embed_params;
    embed_params.useRandomCoords = false;  // Use ETKDG
    embed_params.randomSeed = params.random_seed;
    embed_params.enforceChirality = params.enforce_chirality;
    embed_params.numThreads = params.num_threads;
    embed_params.clearConfs = false;  // Don't clear existing conformers

    RDKit::INT_VECT conf_ids =
        RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, num_confs, embed_params);
    int num_confs_generated = conf_ids.size();

    if (num_confs_generated == 0) {
        std::cerr << "WARNING: Failed to generate ETKDG conformers" << std::endl;
    }
}

void ConformerGenerator::sampleUniformConformers(std::shared_ptr<RDKit::ROMol> mol, int num_confs,
                                                 const ConformerParams& params) {
    if (num_confs <= 0)
        return;

    RDKit::DGeomHelpers::EmbedParameters embed_params;
    embed_params.useRandomCoords = true;  // Use random coordinates
    embed_params.randomSeed = params.random_seed;
    embed_params.enforceChirality = params.enforce_chirality;
    embed_params.numThreads = params.num_threads;
    embed_params.clearConfs = false;  // Don't clear existing conformers

    RDKit::INT_VECT conf_ids =
        RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, num_confs, embed_params);
    int num_confs_generated = conf_ids.size();

    if (num_confs_generated == 0) {
        std::cerr << "WARNING: Failed to generate uniform conformers" << std::endl;
    }
}

std::unordered_map<int, int> ConformerGenerator::clusterConformers(
    std::shared_ptr<RDKit::ROMol> mol, int num_clusters) {
    int num_confs = mol->getNumConformers();

    if (num_confs <= num_clusters) {
        // If we have fewer conformers than clusters, return all
        std::unordered_map<int, int> result;
        for (int i = 0; i < num_confs; ++i) {
            result[i] = i;
        }
        return result;
    }

    // Extract conformer positions
    std::vector<std::vector<double>> positions;
    for (int i = 0; i < num_confs; ++i) {
        const auto& conf = mol->getConformer(i);
        std::vector<double> pos;
        for (unsigned int j = 0; j < conf.getNumAtoms(); ++j) {  // Fixed: changed i to j
            const auto& atom_pos = conf.getAtomPos(j);
            pos.push_back(atom_pos.x);
            pos.push_back(atom_pos.y);
            pos.push_back(atom_pos.z);
        }
        positions.push_back(pos);
    }

    // Perform K-means clustering
    SimpleKMeans kmeans(num_clusters);
    std::vector<int> cluster_labels = kmeans.fit(positions);

    // Select representative conformer from each cluster
    return selectClusterConformers(cluster_labels, num_clusters);
}

std::unordered_map<int, int> ConformerGenerator::selectClusterConformers(
    const std::vector<int>& cluster_labels, int num_clusters) {
    std::unordered_map<int, int> result;
    std::vector<std::vector<int>> clusters(num_clusters);

    // Group conformers by cluster
    for (size_t i = 0; i < cluster_labels.size(); ++i) {
        int cluster_id = cluster_labels[i];
        if (cluster_id >= 0 && cluster_id < num_clusters) {
            clusters[cluster_id].push_back(static_cast<int>(i));
        }
    }

    // Select random conformer from each cluster
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < num_clusters; ++i) {
        if (!clusters[i].empty()) {
            std::uniform_int_distribution<> dis(0, static_cast<int>(clusters[i].size()) - 1);
            int selected_idx = dis(gen);
            result[i] = clusters[i][selected_idx];
        }
    }

    return result;
}

void ConformerGenerator::alignConformers(std::shared_ptr<RDKit::ROMol> mol) {
    int num_confs = mol->getNumConformers();
    if (num_confs <= 1)
        return;

    try {
        RDKit::MolAlign::alignMolConformers(*mol);
    } catch (const std::exception& e) {
        std::cerr << "WARNING: Failed to align conformers: " << e.what() << std::endl;
    }
}

void ConformerGenerator::addHydrogens(std::shared_ptr<RDKit::ROMol> mol) {
    RDKit::MolOps::addHs(*mol);
}
