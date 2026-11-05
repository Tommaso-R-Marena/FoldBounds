// Independent full 3D cubic SAW enumeration and leaf-only contact scoring.
// No symmetry reduction; no shared Python contact-map or incremental scorer.
// Usage: g++ -O3 -std=c++17 experiments/independent_cubic_audit.cpp -o cubic_audit
//        ./cubic_audit 12 > contact_counts.tsv
// Output rows: unsigned_contact_bitmask exact_walk_multiplicity
#include <array>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
using Point = std::array<int, 3>;
const Point directions[6] = {
    {0,0,-1}, {0,-1,0}, {-1,0,0}, {0,0,1}, {0,1,0}, {1,0,0}
};
int n;
std::vector<Point> path;
std::vector<std::pair<int, int>> pairs;
std::map<uint64_t, uint64_t> histogram;
uint64_t leaf_count = 0;
bool occupied(const Point &candidate) {
    for (const Point &old : path) if (candidate == old) return true;
    return false;
}
void visit() {
    if (int(path.size()) == n) {
        uint64_t mask = 0;
        for (size_t b = 0; b < pairs.size(); b++) {
            auto [i, j] = pairs[b];
            const Point &a = path[i], &z = path[j];
            int distance = std::abs(a[0]-z[0])
                         + std::abs(a[1]-z[1])
                         + std::abs(a[2]-z[2]);
            if (distance == 1) mask |= (uint64_t(1) << b);
        }
        histogram[mask]++;
        leaf_count++;
        return;
    }
    const Point current = path.back();
    for (const Point &d : directions) {
        Point nxt = {
            current[0] + d[0], current[1] + d[1], current[2] + d[2]
        };
        if (occupied(nxt)) continue;
        path.push_back(nxt);
        visit();
        path.pop_back();
    }
}
int main(int argc, char **argv) {
    if (argc != 2) return 2;
    n = std::atoi(argv[1]);
    if (n < 2 || n > 14) return 2;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < j - 2; i++)
            if ((j-i) % 2) pairs.emplace_back(i,j);
    if (pairs.size() > 63) return 3;
    path = {{0,0,0}, {1,0,0}};
    visit();
    std::cerr << "n=" << n << " anchored_walks=" << leaf_count
              << " distinct_contact_maps=" << histogram.size() << "\n";
    for (const auto &[mask, count] : histogram)
        std::cout << mask << " " << count << "\n";
    return 0;
}
