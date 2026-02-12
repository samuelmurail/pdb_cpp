#!/usr/bin/env python3
# coding: utf-8

import numpy as np


def rmsd(coor_1, coor_2, selection="name CA", index_list=None, frame_ref=0):
    """Compute RMSD between two sets of coordinates.

    Parameters
    ----------
    coor_1 : Coor
        First set of coordinates.
    coor_2 : Coor
        Second set of coordinates.
    selection : str, optional
        Selection string used when index_list is not provided.
    index_list : list, optional
        Pair of index lists [index_1, index_2].
    frame_ref : int, optional
        Reference frame index in coor_2.

    Returns
    -------
    list[float]
        RMSD values for each model in coor_1.
    """

    if frame_ref < 0 or frame_ref >= coor_2.model_num:
        raise ValueError(
            "Reference frame index is larger than the number of frames in the reference structure"
        )

    if index_list is None:
        index_1 = coor_1.get_index_select(selection)
        index_2 = coor_2.get_index_select(selection)
    else:
        index_1 = index_list[0]
        index_2 = index_list[1]

    if len(index_1) == 0 or len(index_2) == 0:
        raise ValueError("No atoms selected for RMSD calculation")
    if len(index_1) != len(index_2):
        raise ValueError("Index lists must have the same length")

    ref_model = coor_2.models[frame_ref]
    ref_x = np.asarray(ref_model.get_x())
    ref_y = np.asarray(ref_model.get_y())
    ref_z = np.asarray(ref_model.get_z())

    index_1 = np.asarray(index_1, dtype=int)
    index_2 = np.asarray(index_2, dtype=int)

    rmsd_list = []
    for model in coor_1.models:
        x = np.asarray(model.get_x())
        y = np.asarray(model.get_y())
        z = np.asarray(model.get_z())

        dx = x[index_1] - ref_x[index_2]
        dy = y[index_1] - ref_y[index_2]
        dz = z[index_1] - ref_z[index_2]
        rmsd_list.append(float(np.sqrt((dx * dx + dy * dy + dz * dz).mean())))

    return rmsd_list
