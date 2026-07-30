#!/usr/bin/env python3
"""Compare SpeedCHEM real-PLOG ignition history with Cantera and solver peers."""

from __future__ import annotations

import argparse
from pathlib import Path

import cantera as ct
import numpy as np

SPECIES = ("NC16H34", "O2", "OH", "CO", "CO2", "H2O")


def read_history(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("HIST,"):
            rows.append([float(value) for value in line.split(",")[1:]])
    if len(rows) < 3:
        raise ValueError(f"{path}: expected at least 3 HIST rows")
    return np.asarray(rows)


def idt50(history: np.ndarray) -> float:
    threshold = history[0, 1] + 50.0
    crossing = np.flatnonzero(history[:, 1] >= threshold)
    if not len(crossing):
        return float("nan")
    i = int(crossing[0])
    return float(
        np.interp(threshold, history[i - 1 : i + 1, 1], history[i - 1 : i + 1, 0])
    )


def max_dtdt(history: np.ndarray) -> float:
    return float(np.max(np.diff(history[:, 1]) / np.diff(history[:, 0])))


def cantera_history(mechanism: Path, times: np.ndarray) -> np.ndarray:
    gas = ct.Solution(str(mechanism))
    gas.TP = 1000.0, 2.0e6
    gas.set_equivalence_ratio(1.0, "NC16H34:1", "O2:1,N2:3.76")
    reactor = ct.IdealGasReactor(gas, energy="on", clone=False)
    network = ct.ReactorNet([reactor])
    network.rtol = 1.0e-9
    network.atol = 1.0e-15
    indices = [gas.species_index(name) for name in SPECIES]
    rows = []
    for time in times:
        if time > 0.0:
            network.advance(float(time))
        rows.append([time, reactor.T, reactor.thermo.P, *reactor.thermo.X[indices]])
    return np.asarray(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", required=True, type=Path)
    parser.add_argument("--numeric", required=True, type=Path)
    parser.add_argument("--analytic", type=Path)
    args = parser.parse_args()

    numeric = read_history(args.numeric)
    reference = cantera_history(args.mechanism, numeric[:, 0])
    failures: list[str] = []

    speed_idt = idt50(numeric)
    cantera_idt = idt50(reference)
    idt_rel = abs(speed_idt / cantera_idt - 1.0)
    speed_peak = max_dtdt(numeric)
    cantera_peak = max_dtdt(reference)
    peak_rel = abs(speed_peak / cantera_peak - 1.0)
    t_rel_l2 = np.linalg.norm(numeric[:, 1] - reference[:, 1]) / np.linalg.norm(
        reference[:, 1] - reference[0, 1]
    )
    final_t_error = abs(numeric[-1, 1] - reference[-1, 1])

    print(
        "CANTERA idt50_speed_s={:.9e} idt50_ref_s={:.9e} rel={:.3e}".format(
            speed_idt, cantera_idt, idt_rel
        )
    )
    print(
        "CANTERA max_dTdt_speed={:.9e} max_dTdt_ref={:.9e} rel={:.3e}".format(
            speed_peak, cantera_peak, peak_rel
        )
    )
    print(
        f"CANTERA T_history_rel_l2={t_rel_l2:.3e} "
        f"T_final_abs_K={final_t_error:.3e}"
    )

    if not np.isfinite(idt_rel) or idt_rel > 1.0e-2:
        failures.append(f"+50 K IDT relative error {idt_rel:.3e} > 1e-2")
    if not np.isfinite(peak_rel) or peak_rel > 2.5e-1:
        failures.append(f"max dT/dt relative error {peak_rel:.3e} > 0.25")
    if t_rel_l2 > 3.0e-2:
        failures.append(f"temperature-history relative L2 {t_rel_l2:.3e} > 0.03")
    if final_t_error > 1.0:
        failures.append(f"final-temperature error {final_t_error:.3e} K > 1 K")

    for column, name in enumerate(SPECIES, start=3):
        rel_l2 = np.linalg.norm(numeric[:, column] - reference[:, column]) / max(
            np.linalg.norm(reference[:, column]), np.finfo(float).tiny
        )
        final_abs = abs(numeric[-1, column] - reference[-1, column])
        print(
            f"CANTERA species={name} rel_l2={rel_l2:.3e} "
            f"final_abs={final_abs:.3e}"
        )
        if rel_l2 > 1.0e-1:
            failures.append(f"{name} history relative L2 {rel_l2:.3e} > 0.1")
        if final_abs > 1.0e-5:
            failures.append(f"{name} final absolute error {final_abs:.3e} > 1e-5")

    if args.analytic:
        analytic = read_history(args.analytic)
        interpolated = np.column_stack(
            [
                analytic[:, 0],
                *(
                    np.interp(analytic[:, 0], numeric[:, 0], numeric[:, column])
                    for column in range(1, numeric.shape[1])
                ),
            ]
        )
        t_max = float(np.max(np.abs(analytic[:, 1] - interpolated[:, 1])))
        p_rel = float(
            np.max(
                np.abs(analytic[:, 2] - interpolated[:, 2]) / interpolated[:, 2]
            )
        )
        species_max = float(
            np.max(np.abs(analytic[:, 3:9] - interpolated[:, 3:9]))
        )
        print(
            f"SOLVERS T_max_abs_K={t_max:.3e} P_max_rel={p_rel:.3e} "
            f"major_species_max_abs={species_max:.3e}"
        )
        if t_max > 0.1 or p_rel > 1.0e-5 or species_max > 1.0e-6:
            failures.append("numeric/analytic history comparison exceeded tolerance")

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("RESULT: PASS - real PLOG history agrees with independent Cantera reference")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
