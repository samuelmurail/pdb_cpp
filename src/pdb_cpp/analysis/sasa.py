#!/usr/bin/env python3
# coding: utf-8

"""SASA and interface-SASA helpers."""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from ..core import Coor as CoreCoor, compute_sasa as core_compute_sasa

__all__ = ["sasa", "buried_surface_area", "shape_complementarity"]


_VDW_RADII = {
	"H": 1.10,
	"D": 1.10,
	"HE": 1.40,
	"C": 1.70,
	"N": 1.55,
	"O": 1.52,
	"F": 1.47,
	"NE": 1.54,
	"P": 1.80,
	"S": 1.80,
	"CL": 1.75,
	"AR": 1.88,
	"SE": 1.90,
	"BR": 1.85,
	"KR": 2.02,
	"I": 1.98,
	"MG": 1.73,
	"NA": 2.27,
	"K": 2.75,
	"CA": 2.31,
	"MN": 1.73,
	"FE": 1.72,
	"CU": 1.40,
	"ZN": 1.39,
}


def _coor_from_model(model):
	coor = CoreCoor()
	coor.add_Model(model)
	return coor


def _residue_key(entry):
	return (
		entry["chain"],
		entry["resid"],
		entry["insertres"],
		entry["uniq_resid"],
		entry["resname"],
	)


def _infer_element_from_name(atom_name):
	letters = [char.upper() for char in atom_name if char.isalpha()]
	if not letters:
		return "C"
	if len(letters) >= 2:
		pair = "".join(letters[:2])
		if pair in {"CL", "BR", "SE", "NA", "MG", "ZN", "FE", "CA", "MN", "CU"}:
			return pair
	return letters[0]


def _atom_element_symbol(elem, atom_name):
	symbol = elem.upper() if elem else _infer_element_from_name(atom_name)
	return symbol


def _is_hydrogen_like(element):
	return element in {"H", "D", "T"}


def _lookup_vdw_radius(element):
	return _VDW_RADII.get(element, 1.70)


@lru_cache(maxsize=None)
def _generate_sphere_points(n_points):
	points = np.empty((n_points, 3), dtype=float)
	golden_angle = math.pi * (3.0 - math.sqrt(5.0))
	for point_index in range(n_points):
		y_coord = 1.0 - (2.0 * (point_index + 0.5) / n_points)
		radius = math.sqrt(max(0.0, 1.0 - y_coord * y_coord))
		theta = golden_angle * point_index
		points[point_index, 0] = math.cos(theta) * radius
		points[point_index, 1] = y_coord
		points[point_index, 2] = math.sin(theta) * radius
	return points


def _surface_dot_count(radius, dots_per_sq_angstrom):
	area = 4.0 * math.pi * radius * radius
	return max(32, int(round(area * dots_per_sq_angstrom)))


def _make_cell_key(x_coord, y_coord, z_coord, cell_size):
	return (
		int(math.floor(x_coord / cell_size)),
		int(math.floor(y_coord / cell_size)),
		int(math.floor(z_coord / cell_size)),
	)


def _included_atom_coords_and_radii(model, probe_radius, include_hydrogen):
	coords = np.asarray(model.xyz, dtype=float)
	if coords.size == 0:
		return np.empty((0, 3), dtype=float), np.empty((0,), dtype=float)

	included_mask = []
	expanded_radii = []
	for elem, atom_name in zip(model.elem_str, model.name_str):
		element = _atom_element_symbol(elem, atom_name)
		include_atom = include_hydrogen or not _is_hydrogen_like(element)
		included_mask.append(include_atom)
		expanded_radii.append(_lookup_vdw_radius(element) + probe_radius)

	included_mask = np.asarray(included_mask, dtype=bool)
	if not np.any(included_mask):
		return np.empty((0, 3), dtype=float), np.empty((0,), dtype=float)

	return coords[included_mask], np.asarray(expanded_radii, dtype=float)[included_mask]


def _model_surface_points(model, probe_radius, dots_per_sq_angstrom, include_hydrogen):
	coords, expanded_radii = _included_atom_coords_and_radii(
		model,
		probe_radius,
		include_hydrogen,
	)
	if len(coords) == 0:
		return (
			np.empty((0, 3), dtype=float),
			np.empty((0, 3), dtype=float),
			np.empty((0,), dtype=float),
		)

	max_expanded_radius = float(np.max(expanded_radii))
	cell_size = 2.0 * max_expanded_radius

	atom_cells = {}
	for atom_index, coord in enumerate(coords):
		cell_key = _make_cell_key(coord[0], coord[1], coord[2], cell_size)
		atom_cells.setdefault(cell_key, []).append(atom_index)

	surface_points = []
	normals = []
	weights = []

	for atom_index, coord in enumerate(coords):
		radius = float(expanded_radii[atom_index])
		base_cell = _make_cell_key(coord[0], coord[1], coord[2], cell_size)
		candidate_blockers = []
		for dx in (-1, 0, 1):
			for dy in (-1, 0, 1):
				for dz in (-1, 0, 1):
					neighbor_key = (base_cell[0] + dx, base_cell[1] + dy, base_cell[2] + dz)
					for blocker_index in atom_cells.get(neighbor_key, []):
						if blocker_index == atom_index:
							continue
						blocker_radius = float(expanded_radii[blocker_index])
						delta = coord - coords[blocker_index]
						max_overlap = radius + blocker_radius
						if float(np.dot(delta, delta)) < max_overlap * max_overlap:
							candidate_blockers.append(blocker_index)

		point_count = _surface_dot_count(radius, dots_per_sq_angstrom)
		unit_points = _generate_sphere_points(point_count)
		points = coord + (radius * unit_points)

		if candidate_blockers:
			blocked = np.zeros(point_count, dtype=bool)
			for blocker_index in candidate_blockers:
				blocker_delta = points - coords[blocker_index]
				blocked |= np.sum(blocker_delta * blocker_delta, axis=1) < (
					expanded_radii[blocker_index] ** 2
				)
			accessible_mask = ~blocked
		else:
			accessible_mask = np.ones(point_count, dtype=bool)

		if not np.any(accessible_mask):
			continue

		area_per_point = (4.0 * math.pi * radius * radius) / float(point_count)
		surface_points.append(points[accessible_mask])
		normals.append(unit_points[accessible_mask])
		weights.append(np.full(int(np.count_nonzero(accessible_mask)), area_per_point, dtype=float))

	if not surface_points:
		return (
			np.empty((0, 3), dtype=float),
			np.empty((0, 3), dtype=float),
			np.empty((0,), dtype=float),
		)

	return (
		np.concatenate(surface_points, axis=0),
		np.concatenate(normals, axis=0),
		np.concatenate(weights, axis=0),
	)


def _surface_interface_mask(points, blocker_coords, blocker_radii):
	if len(points) == 0 or len(blocker_coords) == 0:
		return np.zeros(len(points), dtype=bool)

	max_blocker_radius = float(np.max(blocker_radii))
	blocker_cells = _build_point_hash(blocker_coords, max_blocker_radius)
	interface_mask = np.zeros(len(points), dtype=bool)

	for point_index, point in enumerate(points):
		base_cell = _make_cell_key(point[0], point[1], point[2], max_blocker_radius)
		for dx in (-1, 0, 1):
			for dy in (-1, 0, 1):
				for dz in (-1, 0, 1):
					neighbor_key = (base_cell[0] + dx, base_cell[1] + dy, base_cell[2] + dz)
					for blocker_index in blocker_cells.get(neighbor_key, ()): 
						delta = point - blocker_coords[blocker_index]
						if float(np.dot(delta, delta)) < float(blocker_radii[blocker_index] ** 2):
							interface_mask[point_index] = True
							break
					if interface_mask[point_index]:
						break
				if interface_mask[point_index]:
					break
			if interface_mask[point_index]:
				break

	return interface_mask


def _build_point_hash(points, cell_size):
	point_cells = {}
	for point_index, point in enumerate(points):
		cell_key = _make_cell_key(point[0], point[1], point[2], cell_size)
		point_cells.setdefault(cell_key, []).append(point_index)
	return point_cells


def _nearest_surface_scores(
	source_points,
	source_normals,
	source_weights,
	target_points,
	target_normals,
	search_radius,
):
	if len(source_points) == 0 or len(target_points) == 0:
		return np.empty((0,), dtype=float), np.empty((0,), dtype=float), np.empty((0,), dtype=float)

	search_radius_sq = search_radius * search_radius
	target_cells = _build_point_hash(target_points, search_radius)
	scores = []
	weights = []
	distances = []

	for point_index, point in enumerate(source_points):
		base_cell = _make_cell_key(point[0], point[1], point[2], search_radius)
		candidate_indices = []
		for dx in (-1, 0, 1):
			for dy in (-1, 0, 1):
				for dz in (-1, 0, 1):
					neighbor_key = (base_cell[0] + dx, base_cell[1] + dy, base_cell[2] + dz)
					candidate_indices.extend(target_cells.get(neighbor_key, ()))

		if not candidate_indices:
			continue

		candidate_points = target_points[candidate_indices]
		deltas = candidate_points - point
		distance_sq = np.sum(deltas * deltas, axis=1)
		nearest_offset = int(np.argmin(distance_sq))
		nearest_distance_sq = float(distance_sq[nearest_offset])
		if nearest_distance_sq > search_radius_sq:
			continue

		nearest_index = candidate_indices[nearest_offset]
		score = float(np.dot(source_normals[point_index], -target_normals[nearest_index]))
		scores.append(max(-1.0, min(1.0, score)))
		weights.append(float(source_weights[point_index]))
		distances.append(math.sqrt(nearest_distance_sq))

	return (
		np.asarray(scores, dtype=float),
		np.asarray(weights, dtype=float),
		np.asarray(distances, dtype=float),
	)


def _weighted_mean(values, weights):
	return float(np.sum(values * weights) / np.sum(weights))


def _weighted_median(values, weights):
	order = np.argsort(values)
	sorted_values = values[order]
	sorted_weights = weights[order]
	cumulative = np.cumsum(sorted_weights)
	half_total = 0.5 * float(np.sum(sorted_weights))
	return float(sorted_values[int(np.searchsorted(cumulative, half_total, side="left"))])


def _weighted_trimmed_mean(values, weights, trim_fraction):
	if trim_fraction <= 0.0:
		return _weighted_mean(values, weights)

	order = np.argsort(values)
	sorted_values = values[order]
	sorted_weights = weights[order]
	total_weight = float(np.sum(sorted_weights))
	lower_bound = trim_fraction * total_weight
	upper_bound = (1.0 - trim_fraction) * total_weight

	kept_sum = 0.0
	kept_weight = 0.0
	running_weight = 0.0
	for value, weight in zip(sorted_values, sorted_weights):
		next_weight = running_weight + float(weight)
		keep_start = max(running_weight, lower_bound)
		keep_end = min(next_weight, upper_bound)
		if keep_end > keep_start:
			segment_weight = keep_end - keep_start
			kept_sum += float(value) * segment_weight
			kept_weight += segment_weight
		running_weight = next_weight

	if kept_weight == 0.0:
		return _weighted_mean(values, weights)
	return kept_sum / kept_weight


def _aggregate_surface_scores(values, weights, reducer, trim_fraction):
	if len(values) == 0:
		return float("nan")
	if reducer == "median":
		return _weighted_median(values, weights)
	if reducer == "trimmed_mean":
		return _weighted_trimmed_mean(values, weights, trim_fraction)
	raise ValueError("reducer must be 'trimmed_mean' or 'median'")


def _polar_atom_mask(model):
	polar_elements = {"N", "O", "S", "P", "SE"}
	elements = []
	for elem, atom_name in zip(model.elem_str, model.name_str):
		symbol = elem.upper() if elem else _infer_element_from_name(atom_name)
		elements.append(symbol)
	return np.asarray([element in polar_elements for element in elements], dtype=bool)


def _group_atom_areas_by_residue(model, atom_areas):
	polar_mask = _polar_atom_mask(model)
	residue_lookup = {}
	residue_areas = []

	for atom_index, atom_area in enumerate(atom_areas):
		key = (
			model.chain_str[atom_index],
			model.resid[atom_index],
			model.insertres_str[atom_index],
			model.uniq_resid[atom_index],
			model.resname_str[atom_index],
		)
		if key not in residue_lookup:
			residue_lookup[key] = len(residue_areas)
			residue_areas.append(
				{
					"chain": model.chain_str[atom_index],
					"resid": model.resid[atom_index],
					"insertres": model.insertres_str[atom_index],
					"uniq_resid": model.uniq_resid[atom_index],
					"resname": model.resname_str[atom_index],
					"area": 0.0,
					"polar_area": 0.0,
					"apolar_area": 0.0,
				}
			)
		residue_areas[residue_lookup[key]]["area"] += float(atom_area)
		if polar_mask[atom_index]:
			residue_areas[residue_lookup[key]]["polar_area"] += float(atom_area)
		else:
			residue_areas[residue_lookup[key]]["apolar_area"] += float(atom_area)

	return residue_areas


def _compute_native_result(model, probe_radius, n_points, include_hydrogen, by_atom, by_residue):
	result = core_compute_sasa(
		model,
		probe_radius=probe_radius,
		n_points=n_points,
		include_hydrogen=include_hydrogen,
		by_atom=True,
	)

	atom_areas = np.asarray(result["atom_areas"], dtype=float)
	polar_mask = _polar_atom_mask(model)
	result["polar"] = float(np.sum(atom_areas[polar_mask]))
	result["apolar"] = float(np.sum(atom_areas[~polar_mask]))

	if by_residue:
		result["residue_areas"] = _group_atom_areas_by_residue(model, atom_areas)

	if not by_atom and "atom_areas" in result:
		del result["atom_areas"]

	return result


def _build_residue_burial(isolated_entries, complex_lookup, partner):
	residue_burial = []
	for entry in isolated_entries:
		complex_entry = complex_lookup.get(
			_residue_key(entry),
			{"area": 0.0, "polar_area": 0.0, "apolar_area": 0.0},
		)
		complex_area = float(complex_entry["area"])
		complex_polar_area = float(complex_entry["polar_area"])
		complex_apolar_area = float(complex_entry["apolar_area"])
		isolated_area = float(entry["area"])
		isolated_polar_area = float(entry["polar_area"])
		isolated_apolar_area = float(entry["apolar_area"])
		residue_burial.append(
			{
				"partner": partner,
				"chain": entry["chain"],
				"resid": entry["resid"],
				"insertres": entry["insertres"],
				"uniq_resid": entry["uniq_resid"],
				"resname": entry["resname"],
				"isolated_area": isolated_area,
				"complex_area": complex_area,
				"buried_area": isolated_area - complex_area,
				"isolated_polar_area": isolated_polar_area,
				"complex_polar_area": complex_polar_area,
				"buried_polar_area": isolated_polar_area - complex_polar_area,
				"isolated_apolar_area": isolated_apolar_area,
				"complex_apolar_area": complex_apolar_area,
				"buried_apolar_area": isolated_apolar_area - complex_apolar_area,
			}
		)
	return residue_burial


def sasa(
	coor,
	selection=None,
	probe_radius=1.4,
	n_points=960,
	include_hydrogen=False,
	by_atom=False,
	by_residue=False,
):
	"""Compute SASA for each model in a Coor object."""
	results = []
	for frame_index in range(coor.model_num):
		if selection is None:
			subset_model = coor.models[frame_index]
			subset_coor = _coor_from_model(subset_model)
		else:
			subset_coor = coor.select_atoms(selection, frame=frame_index)
			subset_model = subset_coor.models[0]

		if subset_model.len == 0:
			raise ValueError("No atoms selected for SASA calculation")

		result = _compute_native_result(
			subset_model,
			probe_radius,
			n_points,
			include_hydrogen,
			by_atom,
			by_residue,
		)
		results.append(result)

	return results


def buried_surface_area(
	coor,
	receptor_sel,
	ligand_sel,
	probe_radius=1.4,
	n_points=960,
	include_hydrogen=False,
	by_residue=False,
):
	"""Compute buried interface surface for each model in a Coor object."""
	results = []
	for frame_index in range(coor.model_num):
		receptor_indices = set(coor.get_index_select(receptor_sel, frame=frame_index))
		ligand_indices = set(coor.get_index_select(ligand_sel, frame=frame_index))

		if not receptor_indices:
			raise ValueError("Receptor selection is empty")
		if not ligand_indices:
			raise ValueError("Ligand selection is empty")
		if receptor_indices & ligand_indices:
			raise ValueError("Receptor and ligand selections must not overlap")

		receptor_coor = coor.select_atoms(receptor_sel, frame=frame_index)
		ligand_coor = coor.select_atoms(ligand_sel, frame=frame_index)
		complex_coor = coor.select_atoms(
			f"({receptor_sel}) or ({ligand_sel})",
			frame=frame_index,
		)

		receptor_result = _compute_native_result(
			receptor_coor.models[0],
			probe_radius,
			n_points,
			include_hydrogen,
			False,
			by_residue,
		)
		ligand_result = _compute_native_result(
			ligand_coor.models[0],
			probe_radius,
			n_points,
			include_hydrogen,
			False,
			by_residue,
		)
		complex_result = _compute_native_result(
			complex_coor.models[0],
			probe_radius,
			n_points,
			include_hydrogen,
			False,
			by_residue,
		)

		receptor_sasa = receptor_result["total"]
		receptor_polar_sasa = receptor_result["polar"]
		receptor_apolar_sasa = receptor_result["apolar"]
		ligand_sasa = ligand_result["total"]
		ligand_polar_sasa = ligand_result["polar"]
		ligand_apolar_sasa = ligand_result["apolar"]
		complex_sasa = complex_result["total"]
		complex_polar_sasa = complex_result["polar"]
		complex_apolar_sasa = complex_result["apolar"]
		buried_surface = receptor_sasa + ligand_sasa - complex_sasa
		buried_polar_surface = receptor_polar_sasa + ligand_polar_sasa - complex_polar_sasa
		buried_apolar_surface = receptor_apolar_sasa + ligand_apolar_sasa - complex_apolar_sasa

		result = {
			"receptor_sasa": receptor_sasa,
			"receptor_polar_sasa": receptor_polar_sasa,
			"receptor_apolar_sasa": receptor_apolar_sasa,
			"ligand_sasa": ligand_sasa,
			"ligand_polar_sasa": ligand_polar_sasa,
			"ligand_apolar_sasa": ligand_apolar_sasa,
			"complex_sasa": complex_sasa,
			"complex_polar_sasa": complex_polar_sasa,
			"complex_apolar_sasa": complex_apolar_sasa,
			"buried_surface": buried_surface,
			"buried_polar_surface": buried_polar_surface,
			"buried_apolar_surface": buried_apolar_surface,
			"interface_area": buried_surface / 2.0,
			"interface_polar_area": buried_polar_surface / 2.0,
			"interface_apolar_area": buried_apolar_surface / 2.0,
			"probe_radius": probe_radius,
			"n_points": n_points,
		}

		if by_residue:
			complex_lookup = {
				_residue_key(entry): {
					"area": float(entry["area"]),
					"polar_area": float(entry["polar_area"]),
					"apolar_area": float(entry["apolar_area"]),
				}
				for entry in complex_result["residue_areas"]
			}
			receptor_burial = _build_residue_burial(
				receptor_result["residue_areas"],
				complex_lookup,
				partner="receptor",
			)
			ligand_burial = _build_residue_burial(
				ligand_result["residue_areas"],
				complex_lookup,
				partner="ligand",
			)
			result["receptor_residue_sasa"] = receptor_result["residue_areas"]
			result["ligand_residue_sasa"] = ligand_result["residue_areas"]
			result["complex_residue_sasa"] = complex_result["residue_areas"]
			result["residue_buried_surface"] = receptor_burial + ligand_burial

		results.append(result)

	return results


def shape_complementarity(
	coor,
	receptor_sel,
	ligand_sel,
	probe_radius=1.4,
	dots_per_sq_angstrom=12.0,
	search_radius=1.5,
	include_hydrogen=False,
	reducer="trimmed_mean",
	trim_fraction=0.1,
):
	"""Estimate Lawrence-Colman style shape complementarity for an interface.

	Surface dots are generated independently for each partner using a rolling-probe
	surface with outward normals. Each interface dot is matched to its closest dot
	on the opposite partner within ``search_radius`` and scored with the normal
	complementarity term ``dot(n_a, -n_b)``.
	"""
	if probe_radius < 0.0:
		raise ValueError("Probe radius must be non-negative")
	if dots_per_sq_angstrom <= 0.0:
		raise ValueError("dots_per_sq_angstrom must be greater than zero")
	if search_radius <= 0.0:
		raise ValueError("search_radius must be greater than zero")
	if not 0.0 <= trim_fraction < 0.5:
		raise ValueError("trim_fraction must satisfy 0 <= trim_fraction < 0.5")

	results = []
	for frame_index in range(coor.model_num):
		receptor_indices = set(coor.get_index_select(receptor_sel, frame=frame_index))
		ligand_indices = set(coor.get_index_select(ligand_sel, frame=frame_index))

		if not receptor_indices:
			raise ValueError("Receptor selection is empty")
		if not ligand_indices:
			raise ValueError("Ligand selection is empty")
		if receptor_indices & ligand_indices:
			raise ValueError("Receptor and ligand selections must not overlap")

		receptor_model = coor.select_atoms(receptor_sel, frame=frame_index).models[0]
		ligand_model = coor.select_atoms(ligand_sel, frame=frame_index).models[0]
		receptor_atom_coords, receptor_atom_radii = _included_atom_coords_and_radii(
			receptor_model,
			probe_radius,
			include_hydrogen,
		)
		ligand_atom_coords, ligand_atom_radii = _included_atom_coords_and_radii(
			ligand_model,
			probe_radius,
			include_hydrogen,
		)

		receptor_points, receptor_normals, receptor_weights = _model_surface_points(
			receptor_model,
			probe_radius,
			dots_per_sq_angstrom,
			include_hydrogen,
		)
		ligand_points, ligand_normals, ligand_weights = _model_surface_points(
			ligand_model,
			probe_radius,
			dots_per_sq_angstrom,
			include_hydrogen,
		)
		receptor_interface_mask = _surface_interface_mask(
			receptor_points,
			ligand_atom_coords,
			ligand_atom_radii,
		)
		ligand_interface_mask = _surface_interface_mask(
			ligand_points,
			receptor_atom_coords,
			receptor_atom_radii,
		)

		receptor_points = receptor_points[receptor_interface_mask]
		receptor_normals = receptor_normals[receptor_interface_mask]
		receptor_weights = receptor_weights[receptor_interface_mask]
		ligand_points = ligand_points[ligand_interface_mask]
		ligand_normals = ligand_normals[ligand_interface_mask]
		ligand_weights = ligand_weights[ligand_interface_mask]

		receptor_scores, receptor_score_weights, receptor_distances = _nearest_surface_scores(
			receptor_points,
			receptor_normals,
			receptor_weights,
			ligand_points,
			ligand_normals,
			search_radius,
		)
		ligand_scores, ligand_score_weights, ligand_distances = _nearest_surface_scores(
			ligand_points,
			ligand_normals,
			ligand_weights,
			receptor_points,
			receptor_normals,
			search_radius,
		)

		if len(receptor_scores) == 0 and len(ligand_scores) == 0:
			raise ValueError("No opposing interface surface dots found within search_radius")

		combined_scores = np.concatenate(
			[array for array in (receptor_scores, ligand_scores) if len(array) > 0]
		)
		combined_weights = np.concatenate(
			[array for array in (receptor_score_weights, ligand_score_weights) if len(array) > 0]
		)

		result = {
			"shape_complementarity": _aggregate_surface_scores(
				combined_scores,
				combined_weights,
				reducer,
				trim_fraction,
			),
			"receptor_shape_complementarity": _aggregate_surface_scores(
				receptor_scores,
				receptor_score_weights,
				reducer,
				trim_fraction,
			),
			"ligand_shape_complementarity": _aggregate_surface_scores(
				ligand_scores,
				ligand_score_weights,
				reducer,
				trim_fraction,
			),
			"receptor_interface_dots": int(len(receptor_scores)),
			"ligand_interface_dots": int(len(ligand_scores)),
			"interface_dot_pairs": int(len(combined_scores)),
			"mean_interface_distance": float(
				np.mean(
					np.concatenate(
						[array for array in (receptor_distances, ligand_distances) if len(array) > 0]
					)
				)
			),
			"probe_radius": probe_radius,
			"dots_per_sq_angstrom": dots_per_sq_angstrom,
			"search_radius": search_radius,
			"reducer": reducer,
			"trim_fraction": trim_fraction,
		}
		results.append(result)

	return results