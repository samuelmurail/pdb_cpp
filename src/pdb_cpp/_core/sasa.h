#ifndef SASA_H
#define SASA_H

/**
 * @file sasa.h
 * @brief Solvent-accessible surface area computation.
 */

#pragma once

#include <cstddef>
#include <vector>

#include "Model.h"

struct SasaResult {
    /** Total solvent-accessible surface area. */
    double total = 0.0;
    /** Optional per-atom SASA values in input atom order. */
    std::vector<double> atom_areas;
};

/** Compute solvent-accessible surface area using a Shrake-Rupley sampler. */
SasaResult compute_sasa(
    const Model &model,
    float probe_radius = 1.4f,
    std::size_t n_points = 960,
    bool include_hydrogen = false);

#endif // SASA_H