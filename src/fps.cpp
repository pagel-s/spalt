#include "fps.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>

double squaredDistance(const Point& a, const Point& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

std::vector<int> farthestPointSampling(const std::vector<Point>& points, size_t k) {
    if (points.empty() || k == 0)
        return {};

    const size_t n = points.size();
    if (k > n)
        k = n;

    std::vector<Point> sampled;
    std::vector<int> sampled_indices;
    sampled.reserve(k);
    sampled_indices.reserve(k);

    std::srand(std::time(nullptr));
    size_t first_idx = std::rand() % n;
    sampled.push_back(points[first_idx]);
    sampled_indices.push_back(static_cast<int>(first_idx));

    std::vector<double> minDist(n, std::numeric_limits<double>::max());

    while (sampled.size() < k) {
        const Point& last = sampled.back();
        for (size_t j = 0; j < n; j++) {
            double d = squaredDistance(points[j], last);
            if (d < minDist[j]) {
                minDist[j] = d;
            }
        }

        size_t farthest_idx = 0;
        double maxDist = -1.0;
        for (size_t j = 0; j < n; j++) {
            if (minDist[j] > maxDist) {
                maxDist = minDist[j];
                farthest_idx = j;
            }
        }
        sampled.push_back(points[farthest_idx]);
        sampled_indices.push_back(static_cast<int>(farthest_idx));
    }

    return sampled_indices;
}