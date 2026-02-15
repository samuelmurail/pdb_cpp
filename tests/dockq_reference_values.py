#!/usr/bin/env python3
# coding: utf-8

"""Reference DockQ outputs captured from DockQ CLI."""

DOCKQ_REFERENCES = {
    "1rxz_colabfold_vs_1rxz": {
        "DockQ": 0.934,
        "Fnat": 0.963,
        "Fnonnat": 0.088,
        "iRMS": 0.618,
        "LRMS": 1.050,
    },
    "model_vs_native": {
        "DockQ": 0.700,
        "Fnat": 0.533,
        "Fnonnat": 0.238,
        "iRMS": 1.232,
        "LRMS": 1.516,
    },
    "1jd4_vs_5m6n": {
        "DockQ": 0.015,
        "Fnat": 0.000,
        "Fnonnat": 1.000,
        "iRMS": 10.453,
        "LRMS": 52.493,
    },
}
