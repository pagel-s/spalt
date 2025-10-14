/**
 * @file fps.h
 * @brief Farthest Point Sampling (FPS) algorithm for vertex selection
 * @author Sebastian
 * @date 2024
 */

#ifndef FPS_H
#define FPS_H

#include <cstddef>

#include <vector>

/**
 * @struct Point
 * @brief Simple 3D point representation
 * 
 * Represents a point in 3D space with x, y, z coordinates.
 * Used for geometric calculations in the FPS algorithm.
 */
struct Point {
    double x; ///< X coordinate
    double y; ///< Y coordinate  
    double z; ///< Z coordinate
};

/**
 * @brief Calculate squared distance between two 3D points
 * 
 * Computes the squared Euclidean distance between two points.
 * Using squared distance avoids expensive square root operations
 * when only relative distances are needed.
 * 
 * @param a First point
 * @param b Second point
 * @return Squared distance between the points
 */
double squaredDistance(const Point& a, const Point& b);

/**
 * @brief Perform Farthest Point Sampling on a set of 3D points
 * 
 * Implements the Farthest Point Sampling algorithm to select a subset of points
 * that are maximally distributed in space. This is useful for selecting a
 * representative subset of vertices from a molecular surface.
 * 
 * Algorithm:
 * 1. Select first point arbitrarily (typically the first point)
 * 2. For each remaining point to select:
 *    - Find the point with maximum minimum distance to already selected points
 *    - Add it to the selected set
 * 
 * @param points Vector of 3D points to sample from
 * @param k Number of points to select (must be <= points.size())
 * @return Vector of indices into the input points vector
 * 
 * @note Returns empty vector if k is 0 or larger than points.size()
 * @note Time complexity: O(k * n) where n is the number of input points
 */
std::vector<int> farthestPointSampling(const std::vector<Point>& points, size_t k);

#endif