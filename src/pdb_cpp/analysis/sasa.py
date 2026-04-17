#!/usr/bin/env python3
# coding: utf-8

"""SASA and interface-SASA helpers."""

from __future__ import annotations

import numpy as np

from ..core import Coor as CoreCoor, compute_sasa as core_compute_sasa

__all__ = ["sasa", "buried_surface_area"]


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