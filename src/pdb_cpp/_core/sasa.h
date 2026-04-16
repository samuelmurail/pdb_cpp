#ifndef SASA_H
#define SASA_H

#pragma once

#include <cstddef>
#include <vector>

#include "Model.h"

struct SasaResult {
    double total = 0.0;
    std::vector<double> atom_areas;
};

SasaResult compute_sasa(
    const Model &model,
    float probe_radius = 1.4f,
    std::size_t n_points = 960,
    bool include_hydrogen = false);

#endif // SASA_H