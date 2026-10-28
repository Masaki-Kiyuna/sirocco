#!/usr/bin/env python3
"""Generate production neutral-H free-free inverse-CDF banks.

The audited kernel and adaptive temperature nodes live in the immutable SSH
handoff directory.  Output is a compact little-endian binary table consumed by
the private SIROCCO implementation.  Each row stores segment, temperature,
physical total power coefficient, mean photon energy, and uniformly spaced
inverse-CDF energy quantiles.
"""

import argparse
import importlib.util
from pathlib import Path
import struct
import sys
import types

import numpy as np

# XD2000's system Python is intentionally minimal.  The audited kernel imports
# SciPy for an attachment-table helper that this generator never calls, so a
# fail-fast stub is sufficient here; all neutral-H ff operations use NumPy.
if not hasattr(np, "trapezoid"):
    np.trapezoid = np.trapz
if "scipy.interpolate" not in sys.modules:
    scipy = types.ModuleType("scipy")
    interpolate = types.ModuleType("scipy.interpolate")

    def _unused_pchip(*args, **kwargs):
        raise RuntimeError("PchipInterpolator is not used by the neutral-H ff generator")

    interpolate.PchipInterpolator = _unused_pchip
    scipy.interpolate = interpolate
    sys.modules["scipy"] = scipy
    sys.modules["scipy.interpolate"] = interpolate


MAGIC = b"HNFFCDF1"
VERSION = 1
EV = 1.602176634e-12
H = 6.62607015e-27
C = 2.99792458e10
POWER_FACTOR = 8.0 * np.pi * EV**4 / (H**3 * C**2)


def load_design(path: Path):
    spec = importlib.util.spec_from_file_location("neutral_h_ff_design", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_nodes(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return [(int(row["segment"]), float(row["T_K"])) for row in data]


def write_bank(output: Path, design, mode: str, rows, nq: int):
    if mode == "john88_guarded":
        local_modes = ("cloudy_segment", "john_segment", "cloudy_segment")
    elif mode == "cloudy_c25_allT":
        local_modes = ("cloudy_c25_allT",)
    else:
        raise ValueError(mode)

    u = np.linspace(0.0, 1.0, nq)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as stream:
        stream.write(MAGIC)
        stream.write(struct.pack("<IIII", VERSION, len(rows), nq, len(local_modes)))
        for segment, temperature in rows:
            pdf, cdf, power_shape, mean = design.distribution(
                float(np.log10(temperature)), local_modes[segment]
            )
            keep = np.r_[True, np.diff(cdf) > 0.0]
            energy = np.interp(u, cdf[keep], design.ENERGY[keep])
            power = power_shape * POWER_FACTOR
            stream.write(struct.pack("<iddd", segment, temperature, power, mean))
            stream.write(np.asarray(energy, dtype="<f8").tobytes())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--handoff",
        type=Path,
        default=Path("/work/kiyunams/sirocco_run/ssh_handoff_hminus_h2_h2plus_20261027"),
    )
    parser.add_argument("--output", type=Path, default=Path("/work/kiyunams/sirocco/data/hminus"))
    args = parser.parse_args()

    design_path = args.handoff / "tools/design_neutral_h_ff_full_tables.py"
    design = load_design(design_path)
    diag = args.handoff / "microphysics/materials/data/ssh_diagnostics"
    jobs = (
        (
            "john88_guarded",
            8071,
            diag / "neutral_H_ff_john88_guarded_temperature_nodes.csv",
        ),
        (
            "cloudy_c25_allT",
            8063,
            diag / "neutral_H_ff_cloudy_c25_allT_temperature_nodes.csv",
        ),
    )
    for mode, nq, node_path in jobs:
        rows = read_nodes(node_path)
        output = args.output / f"neutral_H_ff_{mode}.bin"
        write_bank(output, design, mode, rows, nq)
        print(f"{mode}: {len(rows)} x {nq} -> {output} ({output.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
