#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

#define WITHOUT_NUMPY

#include "matplotlibcpp.h"  

namespace plt = matplotlibcpp;

struct Colors {
    std::string line_left = "blue";
    std::string semi_circle_left = "orange";
    std::string semi_circle_right = "green";
    std::string line_right = "purple";
    std::string f_mapping = "red";
    std::string g_f_mapping = "cyan";
    std::string inverse_mapping = "magenta";
};

void generate_segments(std::vector<std::complex<double>>& line_left,
                       std::vector<std::complex<double>>& semi_circle_left,
                       std::vector<std::complex<double>>& semi_circle_right,
                       std::vector<std::complex<double>>& line_right)
{
    for(int i = 0; i < 100; ++i) {
        double real = -10 + (8.0 * i) / 99.0;  // from -10 to -2
        double imag = 0.01;
        line_left.emplace_back(real, imag);
    }

    for(int i = 0; i < 100; ++i) {
        double theta = M_PI - (M_PI * i) / 99.0;  // from π to 0
        double real = -1 + std::cos(theta);
        double imag = std::sin(theta) + 0.01;
        semi_circle_left.emplace_back(real, imag);
    }

    for(int i = 0; i < 100; ++i) {
        double theta = M_PI - (M_PI * i) / 99.0;  // from π to 0
        double real = 1 + std::cos(theta);
        double imag = std::sin(theta) + 0.01;
        semi_circle_right.emplace_back(real, imag);
    }

    for(int i = 0; i < 100; ++i) {
        double real = 2 + (8.0 * i) / 99.0;  // from 2 to 10
        double imag = 0.01;
        line_right.emplace_back(real, imag);
    }
}

std::vector<std::complex<double>> mapping_f(const std::vector<std::complex<double>>& points) {
    std::vector<std::complex<double>> result;
    result.reserve(points.size());
    for(const auto& p : points) {
        if(p != std::complex<double>(0,0)) {
            result.emplace_back((p - 2.0) / (2.0 * p));
        } else {
            result.emplace_back(0.0, 0.0);
        }
    }
    return result;
}

std::vector<std::complex<double>> mapping_f_inv(const std::vector<std::complex<double>>& w_points) {
    std::vector<std::complex<double>> result;
    result.reserve(w_points.size());
    for(const auto& w : w_points) {
        // Avoid division by zero if (1 - 2w) == 0
        if(w != std::complex<double>(0.5, 0)) {
            result.emplace_back(2.0 / (1.0 - 2.0 * w));
        } else {
            result.emplace_back(0.0, 0.0);
        }
    }
    return result;
}

std::vector<std::complex<double>> mapping_g(const std::vector<std::complex<double>>& w_points) {
    std::vector<std::complex<double>> result;
    result.reserve(w_points.size());
    std::complex<double> factor = M_PI / std::complex<double>(0,1);  // π / i
    for(const auto& w : w_points) {
        result.emplace_back(factor * (w - 1.0));
    }
    return result;
}

std::vector<std::complex<double>> mapping_g_inv(const std::vector<std::complex<double>>& y_points) {
    std::vector<std::complex<double>> result;
    result.reserve(y_points.size());
    std::complex<double> factor = std::complex<double>(0,1) / M_PI;
    for(const auto& y : y_points) {
        result.emplace_back(1.0 + factor * y);
    }
    return result;
}

std::vector<std::complex<double>> mapping_T(const std::vector<std::complex<double>>& z_points) {
    auto f_points = mapping_f(z_points);
    return mapping_g(f_points);
}

std::vector<std::complex<double>> mapping_T_inv(const std::vector<std::complex<double>>& w_points) {
    auto g_inv_points = mapping_g_inv(w_points);
    return mapping_f_inv(g_inv_points);
}

int main() {
    Colors colors;

    std::vector<std::complex<double>> line_left, semi_circle_left, semi_circle_right, line_right;
    generate_segments(line_left, semi_circle_left, semi_circle_right, line_right);

    std::vector<std::complex<double>> original_points;
    original_points.reserve(line_left.size() + semi_circle_left.size()
                            + semi_circle_right.size() + line_right.size());
    original_points.insert(original_points.end(), line_left.begin(), line_left.end());
    original_points.insert(original_points.end(), semi_circle_left.begin(), semi_circle_left.end());
    original_points.insert(original_points.end(), semi_circle_right.begin(), semi_circle_right.end());
    original_points.insert(original_points.end(), line_right.begin(), line_right.end());

    auto f_points = mapping_f(original_points);
    auto g_f_points = mapping_g(f_points);
    auto inverted_points = mapping_T_inv(g_f_points); // Should match original_points

    std::vector<double> orig_x, orig_y;
    std::vector<double> f_x, f_y;
    std::vector<double> gf_x, gf_y;
    std::vector<double> inv_x, inv_y;

    for (size_t i = 0; i < original_points.size(); ++i) {
        orig_x.push_back(original_points[i].real());
        orig_y.push_back(original_points[i].imag());

        f_x.push_back(f_points[i].real());
        f_y.push_back(f_points[i].imag());

        gf_x.push_back(g_f_points[i].real());
        gf_y.push_back(g_f_points[i].imag());

        inv_x.push_back(inverted_points[i].real());
        inv_y.push_back(inverted_points[i].imag());
    }

    // Plot 1: Original segments
    plt::figure();
    plt::scatter(orig_x, orig_y, 1.0, {{"c", "black"}});
    plt::title("Исходные отрезки и дуги");
    plt::xlabel("Re(z)");
    plt::ylabel("Im(z)");
    plt::xlim(-12, 12);
    plt::ylim(-2, 12);
    plt::grid(true);
    plt::axis("equal");
    plt::show();

    // Plot 2: After f(z) = (z - 2)/(2z)
    plt::figure();
    plt::scatter(f_x, f_y, 1.0, {{"c", colors.f_mapping}});
    plt::title("Отображение f(z) = (z - 2)/(2z)");
    plt::xlabel("Re(f(z))");
    plt::ylabel("Im(f(z))");
    plt::xlim(-5, 5);
    plt::ylim(-5, 15);
    plt::grid(true);
    plt::axis("equal");
    plt::show();

    // Plot 3: After composition
    plt::figure();
    plt::scatter(gf_x, gf_y, 1.0, {{"c", colors.g_f_mapping}});
    plt::title("Отображение g(f(z)) = (π / i) * (f(z) - 1)");
    plt::xlabel("Re(g(f(z)))");
    plt::ylabel("Im(g(f(z)))");
    plt::xlim(-20, 20);
    plt::ylim(-20, 20);
    plt::grid(true);
    plt::axis("equal");
    plt::show();

    // Plot 4: After inverse composition
    plt::figure();
    plt::scatter(inv_x, inv_y, 1.0, {{"c", colors.inverse_mapping}});
    plt::title("Обратное отображение T⁻¹(g(f(z)))");
    plt::xlabel("Re(T⁻¹(g(f(z))))");
    plt::ylabel("Im(T⁻¹(g(f(z))))");
    plt::xlim(-12, 12);
    plt::ylim(-2, 12);
    plt::grid(true);
    plt::axis("equal");
    plt::show();

    return 0;
}
