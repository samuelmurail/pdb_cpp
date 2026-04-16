#include "sasa.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {

constexpr float PI_F = 3.14159265358979323846f;

struct CellKey {
    int x;
    int y;
    int z;

    bool operator==(const CellKey &other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct CellKeyHash {
    std::size_t operator()(const CellKey &key) const {
        std::size_t seed = std::hash<int>{}(key.x);
        seed ^= std::hash<int>{}(key.y) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<int>{}(key.z) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

std::string trim_upper(const std::array<char, 5> &value) {
    std::string out;
    out.reserve(value.size());
    for (char ch : value) {
        if (ch == '\0' || ch == ' ') {
            continue;
        }
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
    }
    return out;
}

std::string infer_element_from_name(const std::array<char, 5> &name) {
    std::string letters;
    letters.reserve(name.size());
    for (char ch : name) {
        if (std::isalpha(static_cast<unsigned char>(ch))) {
            letters.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
        }
    }
    if (letters.empty()) {
        return "C";
    }
    if (letters.size() >= 2) {
        const std::string pair = letters.substr(0, 2);
        if (pair == "CL" || pair == "BR" || pair == "SE" || pair == "NA" || pair == "MG" ||
            pair == "ZN" || pair == "FE" || pair == "CA" || pair == "MN" || pair == "CU") {
            return pair;
        }
    }
    return letters.substr(0, 1);
}

std::string get_element_symbol(const Model &model, std::size_t index) {
    std::string elem = trim_upper(model.get_elem()[index]);
    if (!elem.empty()) {
        return elem;
    }
    return infer_element_from_name(model.get_name()[index]);
}

bool is_hydrogen_like(const std::string &element) {
    return element == "H" || element == "D" || element == "T";
}

float lookup_vdw_radius(const std::string &element) {
    static const std::unordered_map<std::string, float> radii = {
        {"H", 1.10f}, {"D", 1.10f}, {"HE", 1.40f}, {"C", 1.70f}, {"N", 1.55f},
        {"O", 1.52f}, {"F", 1.47f}, {"NE", 1.54f}, {"P", 1.80f}, {"S", 1.80f},
        {"CL", 1.75f}, {"AR", 1.88f}, {"SE", 1.90f}, {"BR", 1.85f}, {"KR", 2.02f},
        {"I", 1.98f}, {"MG", 1.73f}, {"NA", 2.27f}, {"K", 2.75f}, {"CA", 2.31f},
        {"MN", 1.73f}, {"FE", 1.72f}, {"CU", 1.40f}, {"ZN", 1.39f},
    };

    const auto found = radii.find(element);
    if (found != radii.end()) {
        return found->second;
    }
    return 1.70f;
}

std::vector<std::array<float, 3>> generate_sphere_points(std::size_t n_points) {
    std::vector<std::array<float, 3>> points;
    points.reserve(n_points);
    const float golden_angle = PI_F * (3.0f - std::sqrt(5.0f));

    for (std::size_t i = 0; i < n_points; ++i) {
        const float y = 1.0f - (2.0f * (static_cast<float>(i) + 0.5f) / static_cast<float>(n_points));
        const float r = std::sqrt(std::max(0.0f, 1.0f - y * y));
        const float theta = golden_angle * static_cast<float>(i);
        points.push_back({std::cos(theta) * r, y, std::sin(theta) * r});
    }

    return points;
}

CellKey make_cell_key(float x, float y, float z, float cell_size) {
    return CellKey{
        static_cast<int>(std::floor(x / cell_size)),
        static_cast<int>(std::floor(y / cell_size)),
        static_cast<int>(std::floor(z / cell_size)),
    };
}

} // namespace

SasaResult compute_sasa(const Model &model, float probe_radius, std::size_t n_points, bool include_hydrogen) {
    if (probe_radius < 0.0f) {
        throw std::runtime_error("Probe radius must be non-negative");
    }
    if (n_points == 0) {
        throw std::runtime_error("n_points must be greater than zero");
    }

    const std::size_t atom_count = model.size();
    SasaResult result;
    result.atom_areas.assign(atom_count, 0.0f);
    if (atom_count == 0) {
        return result;
    }

    const auto sphere_points = generate_sphere_points(n_points);
    const auto &xs = model.get_x();
    const auto &ys = model.get_y();
    const auto &zs = model.get_z();

    std::vector<bool> included(atom_count, true);
    std::vector<float> expanded_radii(atom_count, 0.0f);
    float max_expanded_radius = 0.0f;

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const std::string element = get_element_symbol(model, atom_index);
        if (!include_hydrogen && is_hydrogen_like(element)) {
            included[atom_index] = false;
            continue;
        }
        expanded_radii[atom_index] = lookup_vdw_radius(element) + probe_radius;
        max_expanded_radius = std::max(max_expanded_radius, expanded_radii[atom_index]);
    }

    if (max_expanded_radius <= 0.0f) {
        return result;
    }

    const float cell_size = 2.0f * max_expanded_radius;
    std::unordered_map<CellKey, std::vector<int>, CellKeyHash> cells;
    cells.reserve(atom_count * 2);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        if (!included[atom_index]) {
            continue;
        }
        cells[make_cell_key(xs[atom_index], ys[atom_index], zs[atom_index], cell_size)].push_back(
            static_cast<int>(atom_index));
    }

    std::vector<int> candidate_blockers;
    candidate_blockers.reserve(128);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        if (!included[atom_index]) {
            continue;
        }

        const float x_i = xs[atom_index];
        const float y_i = ys[atom_index];
        const float z_i = zs[atom_index];
        const float radius_i = expanded_radii[atom_index];
        const float radius_i_sq = radius_i * radius_i;
        const CellKey base_cell = make_cell_key(x_i, y_i, z_i, cell_size);

        candidate_blockers.clear();
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    const CellKey neighbor_cell{base_cell.x + dx, base_cell.y + dy, base_cell.z + dz};
                    const auto found = cells.find(neighbor_cell);
                    if (found == cells.end()) {
                        continue;
                    }
                    for (int blocker_index : found->second) {
                        if (static_cast<std::size_t>(blocker_index) == atom_index) {
                            continue;
                        }
                        const float radius_j = expanded_radii[blocker_index];
                        const float dx_ij = x_i - xs[blocker_index];
                        const float dy_ij = y_i - ys[blocker_index];
                        const float dz_ij = z_i - zs[blocker_index];
                        const float max_overlap = radius_i + radius_j;
                        if ((dx_ij * dx_ij) + (dy_ij * dy_ij) + (dz_ij * dz_ij) < max_overlap * max_overlap) {
                            candidate_blockers.push_back(blocker_index);
                        }
                    }
                }
            }
        }

        std::size_t accessible_points = 0;
        for (const auto &unit_point : sphere_points) {
            const float px = x_i + radius_i * unit_point[0];
            const float py = y_i + radius_i * unit_point[1];
            const float pz = z_i + radius_i * unit_point[2];
            bool blocked = false;

            for (int blocker_index : candidate_blockers) {
                const float dx = px - xs[blocker_index];
                const float dy = py - ys[blocker_index];
                const float dz = pz - zs[blocker_index];
                const float blocker_radius = expanded_radii[blocker_index];
                if ((dx * dx) + (dy * dy) + (dz * dz) < blocker_radius * blocker_radius) {
                    blocked = true;
                    break;
                }
            }

            if (!blocked) {
                ++accessible_points;
            }
        }

        const double area =
            (static_cast<double>(accessible_points) / static_cast<double>(n_points)) *
            (4.0 * static_cast<double>(PI_F) * static_cast<double>(radius_i_sq));
        result.atom_areas[atom_index] = area;
        result.total += area;
    }

    return result;
}