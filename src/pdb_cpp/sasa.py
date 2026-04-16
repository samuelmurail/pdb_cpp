#!/usr/bin/env python3
# coding: utf-8

"""Helpers for solvent-accessible surface area calculations."""

__all__ = ["buried_surface_area"]


def _residue_key(entry):
    return (
        entry["chain"],
        entry["resid"],
        entry["insertres"],
        entry["uniq_resid"],
        entry["resname"],
    )


def _build_residue_burial(isolated_entries, complex_lookup, partner):
    residue_burial = []
    for entry in isolated_entries:
        complex_area = float(complex_lookup.get(_residue_key(entry), 0.0))
        isolated_area = float(entry["area"])
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
            }
        )
    return residue_burial


def buried_surface_area(
    coor,
    receptor_sel,
    ligand_sel,
    frame=0,
    probe_radius=1.4,
    n_points=960,
    include_hydrogen=False,
    by_residue=False,
):
    """Compute buried surface for two non-overlapping selections.

    Parameters
    ----------
    coor : Coor
        Structure containing both partners.
    receptor_sel : str
        Selection string for the receptor partner.
    ligand_sel : str
        Selection string for the ligand partner.
    frame : int, default=0
        Model frame to evaluate.
    probe_radius : float, default=1.4
        Solvent probe radius in Angstrom.
    n_points : int, default=960
        Number of Shrake-Rupley sphere points per atom.
    include_hydrogen : bool, default=False
        Whether to include hydrogen-like atoms.
    by_residue : bool, default=False
        Whether to include per-residue SASA and burial breakdowns.

    Returns
    -------
    dict
        SASA totals for receptor, ligand, and complex, plus buried surface.
    """

    receptor_indices = set(coor.get_index_select(receptor_sel, frame=frame))
    ligand_indices = set(coor.get_index_select(ligand_sel, frame=frame))

    if not receptor_indices:
        raise ValueError("Receptor selection is empty")
    if not ligand_indices:
        raise ValueError("Ligand selection is empty")
    if receptor_indices & ligand_indices:
        raise ValueError("Receptor and ligand selections must not overlap")

    receptor = coor.select_atoms(receptor_sel, frame=frame).models[0]
    ligand = coor.select_atoms(ligand_sel, frame=frame).models[0]
    complex_model = coor.select_atoms(
        f"({receptor_sel}) or ({ligand_sel})",
        frame=frame,
    ).models[0]

    receptor_result = receptor.sasa(
        probe_radius=probe_radius,
        n_points=n_points,
        include_hydrogen=include_hydrogen,
        by_residue=by_residue,
    )
    ligand_result = ligand.sasa(
        probe_radius=probe_radius,
        n_points=n_points,
        include_hydrogen=include_hydrogen,
        by_residue=by_residue,
    )
    complex_result = complex_model.sasa(
        probe_radius=probe_radius,
        n_points=n_points,
        include_hydrogen=include_hydrogen,
        by_residue=by_residue,
    )
    receptor_sasa = receptor_result["total"]
    ligand_sasa = ligand_result["total"]
    complex_sasa = complex_result["total"]
    buried_surface = receptor_sasa + ligand_sasa - complex_sasa

    result = {
        "receptor_sasa": receptor_sasa,
        "ligand_sasa": ligand_sasa,
        "complex_sasa": complex_sasa,
        "buried_surface": buried_surface,
        "interface_area": buried_surface / 2.0,
        "probe_radius": probe_radius,
        "n_points": n_points,
    }

    if by_residue:
        complex_lookup = {
            _residue_key(entry): float(entry["area"])
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

    return result