#!/usr/bin/env python3
# coding: utf-8

from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import Coor, analysis


def _parse_chain_map(value: str) -> dict[str, str]:
    chain_map = {}
    if not value.strip():
        raise argparse.ArgumentTypeError("chain map cannot be empty")

    for item in value.split(","):
        native_model = item.strip()
        if not native_model:
            continue
        if ":" not in native_model:
            raise argparse.ArgumentTypeError(
                "chain map must use native:model pairs, for example A:B,B:A,C:C"
            )
        native_chain, model_chain = (part.strip() for part in native_model.split(":", 1))
        if not native_chain or not model_chain:
            raise argparse.ArgumentTypeError(
                "chain map entries must use non-empty native:model pairs"
            )
        chain_map[native_chain] = model_chain

    if not chain_map:
        raise argparse.ArgumentTypeError("chain map cannot be empty")
    return chain_map


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute DockQ with pdb_cpp using analysis.dockQ_multimer(), "
            "including automatic native->model chain mapping."
        )
    )
    parser.add_argument("model", help="Path to the model structure (.pdb/.cif/.pqr/.gro)")
    parser.add_argument("native", help="Path to the native structure (.pdb/.cif/.pqr/.gro)")
    parser.add_argument(
        "--chain-map",
        type=_parse_chain_map,
        help="Explicit native:model chain mapping, for example A:B,B:A,C:C",
    )
    parser.add_argument(
        "--n-cpu",
        type=int,
        default=1,
        help="Number of CPUs passed to analysis.dockQ_multimer()",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print machine-readable JSON instead of a plain-text report",
    )
    return parser


def _float_or_none(value: float | None) -> float | None:
    if value is None:
        return None
    return float(value)


def _serialise_result(result: dict) -> dict:
    interfaces = {}
    for (native_chain_1, native_chain_2), iface_result in result["interfaces"].items():
        key = f"{native_chain_1}-{native_chain_2}"
        if iface_result is None:
            interfaces[key] = None
            continue

        interfaces[key] = {
            "model_rec_chain": iface_result["model_rec_chain"],
            "model_lig_chain": iface_result["model_lig_chain"],
            "DockQ": [_float_or_none(v) for v in iface_result["DockQ"]],
            "Fnat": [_float_or_none(v) for v in iface_result["Fnat"]],
            "Fnonnat": [_float_or_none(v) for v in iface_result["Fnonnat"]],
            "LRMS": [_float_or_none(v) for v in iface_result["LRMS"]],
            "iRMS": [_float_or_none(v) for v in iface_result["iRMS"]],
            "rRMS": [_float_or_none(v) for v in iface_result["rRMS"]],
            "clashes": [int(v) for v in iface_result["clashes"]],
        }

    return {
        "chain_map": dict(result["chain_map"]),
        "GlobalDockQ": [_float_or_none(v) for v in result["GlobalDockQ"]],
        "interfaces": interfaces,
    }


def _format_value(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.3f}"


def _print_report(result: dict, model_path: str, native_path: str) -> None:
    print(f"model: {Path(model_path)}")
    print(f"native: {Path(native_path)}")
    print(
        "chain_map: "
        + ", ".join(f"{native}:{model}" for native, model in result["chain_map"].items())
    )

    if result["GlobalDockQ"]:
        print(f"GlobalDockQ: {_format_value(result['GlobalDockQ'][0])}")
    else:
        print("GlobalDockQ: NA")

    print("")
    print("native_if	model_if	DockQ	Fnat	Fnonnat	LRMS	iRMS	rRMS	clashes")
    for (native_chain_1, native_chain_2), iface_result in result["interfaces"].items():
        if iface_result is None:
            print(f"{native_chain_1}-{native_chain_2}	NA	NA	NA	NA	NA	NA	NA	NA")
            continue

        model_if = f"{iface_result['model_rec_chain']}-{iface_result['model_lig_chain']}"
        print(
            f"{native_chain_1}-{native_chain_2}\t"
            f"{model_if}\t"
            f"{_format_value(iface_result['DockQ'][0])}\t"
            f"{_format_value(iface_result['Fnat'][0])}\t"
            f"{_format_value(iface_result['Fnonnat'][0])}\t"
            f"{_format_value(iface_result['LRMS'][0])}\t"
            f"{_format_value(iface_result['iRMS'][0])}\t"
            f"{_format_value(iface_result['rRMS'][0])}\t"
            f"{int(iface_result['clashes'][0])}"
        )


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.n_cpu < 1:
        parser.error("--n-cpu must be >= 1")

    model = Coor(args.model)
    native = Coor(args.native)
    result = analysis.dockQ_multimer(
        model,
        native,
        chain_map=args.chain_map,
        n_cpu=args.n_cpu,
    )

    if args.json:
        print(json.dumps(_serialise_result(result), indent=2, sort_keys=True))
    else:
        _print_report(result, args.model, args.native)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())