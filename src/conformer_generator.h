#pragma once
#include "conformer_params.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Conformer.h>
#include <unordered_map>
#include <vector>
#include <memory>

/**
 * @brief Advanced conformer generation with clustering and multiple sampling strategies
 * 
 * Based on the Python reference implementation, this class provides:
 * - Multiple conformer sampling strategies (optimized, ETKDG, uniform)
 * - K-means clustering for diverse conformer selection
 * - MMFF optimization for energy minimization
 * - Configurable parameters for fine-tuning generation
 */
class ConformerGenerator {
public:
    /**
     * @brief Constructor with default parameters
     */
    ConformerGenerator() = default;
    
    /**
     * @brief Constructor with custom parameters
     * @param params Conformer generation parameters
     */
    explicit ConformerGenerator(const ConformerParams& params);
    
    /**
     * @brief Generate conformers for a molecule using multiple strategies
     * @param mol RDKit molecule (will be modified with conformers)
     * @param params Optional custom parameters (uses member params if not provided)
     * @return Map of cluster_id -> conformer_id for selected conformers
     */
    std::unordered_map<int, int> generateConformers(
        std::shared_ptr<RDKit::ROMol> mol,
        const ConformerParams* params = nullptr
    );
    
    /**
     * @brief Generate conformers from SMILES string
     * @param smiles SMILES string
     * @param params Optional custom parameters
     * @return Tuple of (molecule, cluster_map)
     */
    std::pair<std::shared_ptr<RDKit::ROMol>, std::unordered_map<int, int>> 
    generateConformersFromSmiles(
        const std::string& smiles,
        const ConformerParams* params = nullptr
    );
    
    /**
     * @brief Get current parameters
     */
    const ConformerParams& getParams() const { return params_; }
    
    /**
     * @brief Update parameters
     */
    void setParams(const ConformerParams& params) { params_ = params; }

private:
    ConformerParams params_;
    
    /**
     * @brief Sample conformers with MMFF optimization (4/6 of total)
     */
    void sampleOptimizedConformers(
        std::shared_ptr<RDKit::ROMol> mol,
        int num_confs,
        const ConformerParams& params
    );
    
    /**
     * @brief Sample conformers using ETKDG method (1/6 of total)
     */
    void sampleETKDGConformers(
        std::shared_ptr<RDKit::ROMol> mol,
        int num_confs,
        const ConformerParams& params
    );
    
    /**
     * @brief Sample conformers with uniform random coordinates (1/6 of total)
     */
    void sampleUniformConformers(
        std::shared_ptr<RDKit::ROMol> mol,
        int num_confs,
        const ConformerParams& params
    );
    
    /**
     * @brief Cluster conformers using K-means and select representative from each cluster
     * @param mol Molecule with conformers
     * @param num_clusters Number of clusters
     * @return Map of cluster_id -> conformer_id
     */
    std::unordered_map<int, int> clusterConformers(
        std::shared_ptr<RDKit::ROMol> mol,
        int num_clusters
    );
    
    /**
     * @brief Select one conformer from each cluster (random selection)
     * @param cluster_labels Cluster assignment for each conformer
     * @param num_clusters Number of clusters
     * @return Map of cluster_id -> conformer_id
     */
    std::unordered_map<int, int> selectClusterConformers(
        const std::vector<int>& cluster_labels,
        int num_clusters
    );
    
    /**
     * @brief Align all conformers to the first one
     */
    void alignConformers(std::shared_ptr<RDKit::ROMol> mol);
    
    /**
     * @brief Add hydrogens to molecule if not present
     */
    void addHydrogens(std::shared_ptr<RDKit::ROMol> mol);
};
