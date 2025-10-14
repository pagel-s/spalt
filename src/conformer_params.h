#pragma once
#include <string>
#include <vector>

/**
 * @brief Parameters for advanced conformer generation and clustering
 * 
 * Comprehensive parameter structure for the advanced conformer generation algorithm
 * that uses multiple sampling strategies with K-means clustering to select diverse,
 * high-quality conformers. Based on the Python reference implementation.
 */
struct ConformerParams {
    // Basic conformer generation
    int total_conformers = 50;           ///< Total number of conformers to generate across all sampling methods (default: 50)
    int num_clusters = 5;                ///< Number of clusters for K-means clustering to select diverse conformers (default: 5)
    int num_threads = 1;                 ///< Number of threads for parallel conformer generation (default: 1)
    
    // Sampling weights (must sum to 1.0)
    double optimized_weight = 4.0/6.0;   ///< Weight for MMFF-optimized conformers (default: 4/6 ≈ 0.667)
    double etkdg_weight = 1.0/6.0;       ///< Weight for ETKDG (Experimental-Torsion-angle preference with Distance Geometry) conformers (default: 1/6 ≈ 0.167)
    double uniform_weight = 1.0/6.0;     ///< Weight for uniform random conformers (default: 1/6 ≈ 0.167)
    
    // Algorithm parameters
    int random_seed = 4;                 ///< Random seed for reproducible conformer generation (default: 4)
    bool enforce_chirality = true;       ///< Whether to enforce stereochemistry during conformer generation (default: true)
    bool align_conformers = true;        ///< Whether to align conformers to principal axes after generation (default: true)
    bool use_random_coords = false;      ///< Whether to use random initial coordinates instead of distance geometry (default: false)
    
    // Optimization parameters
    bool enable_optimization = false;    ///< Whether to run MMFF (Merck Molecular Force Field) optimization on conformers (default: false)
    bool clear_existing = false;         ///< Whether to clear existing conformers before generating new ones (default: false)
    
    // Fallback parameters (when optimization fails)
    double fallback_etkdg_weight = 0.5;  ///< ETKDG weight to use when MMFF optimization fails (default: 0.5)
    double fallback_uniform_weight = 0.5; ///< Uniform weight to use when MMFF optimization fails (default: 0.5)
    
    /**
     * @brief Validate and normalize weights
     */
    void normalizeWeights() {
        double total = optimized_weight + etkdg_weight + uniform_weight;
        if (total > 0.0) {
            optimized_weight /= total;
            etkdg_weight /= total;
            uniform_weight /= total;
        }
    }
    
    /**
     * @brief Get number of conformers for each sampling method
     */
    std::vector<int> getConformerCounts() const {
        std::vector<int> counts;
        counts.push_back(static_cast<int>(total_conformers * optimized_weight));
        counts.push_back(static_cast<int>(total_conformers * etkdg_weight));
        counts.push_back(static_cast<int>(total_conformers * uniform_weight));
        return counts;
    }
    
    /**
     * @brief Apply fallback weights (when optimization fails)
     */
    void applyFallbackWeights() {
        optimized_weight = 0.0;
        etkdg_weight = fallback_etkdg_weight;
        uniform_weight = fallback_uniform_weight;
        normalizeWeights();
    }
};
