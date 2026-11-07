#include "quantum.h"
#include "sampling.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (argc < 5 || argc > 7) {
        std::cerr << "Usage: export_cloud OUTPUT.csv n l basis_index [points=70000] [seed=42]\n";
        return 1;
    }
    try {
        const int n = std::stoi(argv[2]), l = std::stoi(argv[3]), m = std::stoi(argv[4]);
        const int count = argc > 5 ? std::stoi(argv[5]) : 70000;
        const auto seed = argc > 6 ? static_cast<unsigned>(std::stoul(argv[6])) : 42u;
        if (n < 1 || n > 4 || l < 0 || l >= n || std::abs(m) > l || count < 1 || count > 2000000)
            throw std::invalid_argument("Expected n=1..4, l=0..n-1, |basis_index|<=l, points=1..2000000");
        std::ofstream out(argv[1]);
        if (!out) throw std::runtime_error("Cannot open output file");
        const auto points = generate_orbital(n, l, m, count, seed);
        out << "x,y,z,sign\n" << std::setprecision(9);
        for (const auto& p : points) out << p.x << ',' << p.y << ',' << p.z << ',' << p.sign << '\n';
        if (!out) throw std::runtime_error("Could not write complete cloud");
        std::cout << "Saved " << points.size() << " real-orbital points to " << argv[1] << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
