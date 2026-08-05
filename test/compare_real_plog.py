#!/usr/bin/env python3
"""Independent Cantera checks for the pinned public C3Mech PLOG mechanism."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

import cantera as ct
import numpy as np

HISTORY_SPECIES = ("H2", "O2", "H", "OH", "HO2", "H2O", "AR")
STATIC_COMPOSITION = {"H2": 2.0, "O2": 1.0, "AR": 7.0, "H": 1.0e-5, "HO2": 1.0e-8}
HISTORY_COMPOSITION = {"H2": 2.0, "O2": 1.0, "AR": 7.0}


def is_plog_reaction(reaction: ct.Reaction) -> bool:
    return (
        getattr(reaction, "reaction_type", "") == "pressure-dependent-Arrhenius"
        or type(reaction.rate).__name__ == "PlogRate"
    )


def plog_indices(gas: ct.Solution) -> list[int]:
    return [i for i, reaction in enumerate(gas.reactions()) if is_plog_reaction(reaction)]


def read_rate_rows(path: Path) -> list[tuple[int, int, float, float, float]]:
    rows: list[tuple[int, int, float, float, float]] = []
    for line in path.read_text().splitlines():
        if not line.startswith("RATE,"):
            continue
        fields = next(csv.reader([line]))
        rows.append(
            (int(fields[1]), int(fields[2]), float(fields[3]), float(fields[4]), float(fields[5]))
        )
    if not rows:
        raise ValueError(f"{path}: no RATE rows")
    return rows


def compare_rate_grid(
    gas: ct.Solution,
    published: ct.Solution,
    rows: Iterable[tuple[int, int, float, float, float]],
) -> list[str]:
    failures: list[str] = []
    known_plog = set(plog_indices(gas))
    published_plog = set(plog_indices(published))
    max_speed_rel = 0.0
    max_published_rel = 0.0
    worst: tuple[int, float, float] | None = None
    count = 0

    for packed, reaction_one_based, temperature, pressure, speed_k_ck in rows:
        del packed
        reaction_index = reaction_one_based - 1
        if reaction_index not in known_plog:
            failures.append(
                f"SpeedCHEM marked reaction {reaction_one_based} as PLOG, but Cantera did not"
            )
            continue
        if reaction_index not in published_plog:
            failures.append(
                f"published C3Mech YAML lost PLOG reaction {reaction_one_based}"
            )
            continue

        reaction = gas.reaction(reaction_index)
        published_reaction = published.reaction(reaction_index)
        if reaction.equation != published_reaction.equation:
            failures.append(
                f"reaction-order mismatch at {reaction_one_based}: "
                f"{reaction.equation!r} != {published_reaction.equation!r}"
            )
            continue

        orders = dict(getattr(reaction, "orders", {}))
        forward_order = (
            float(sum(orders.values()))
            if orders
            else float(sum(reaction.reactants.values()))
        )
        # SpeedCHEM/CHEMKIN uses mol/cm^3; Cantera's Python API uses kmol/m^3.
        # k_SI = k_CK * 1000**(1-order).
        speed_k_si = speed_k_ck * 1000.0 ** (1.0 - forward_order)

        gas.TP = temperature, pressure
        published.TP = temperature, pressure
        cantera_k = float(gas.forward_rate_constants[reaction_index])
        published_k = float(published.forward_rate_constants[reaction_index])

        speed_rel = abs(speed_k_si - cantera_k) / max(abs(cantera_k), np.finfo(float).tiny)
        published_rel = abs(published_k - cantera_k) / max(
            abs(cantera_k), np.finfo(float).tiny
        )
        if speed_rel > max_speed_rel:
            max_speed_rel = speed_rel
            worst = (reaction_one_based, temperature, pressure)
        max_published_rel = max(max_published_rel, published_rel)
        count += 1

    print(
        f"PLOG_RATE rows={count} max_speedchem_rel={max_speed_rel:.3e} "
        f"max_published_yaml_rel={max_published_rel:.3e} worst={worst}"
    )
    if max_speed_rel > 5.0e-4:
        failures.append(f"PLOG rate-grid relative error {max_speed_rel:.3e} > 5e-4")
    if max_published_rel > 5.0e-12:
        failures.append(
            f"published/regenerated Cantera PLOG mismatch {max_published_rel:.3e} > 5e-12"
        )
    return failures


def read_state_grid(path: Path) -> tuple[list[str], np.ndarray]:
    species: list[str] | None = None
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        if line.startswith("SPECIES,"):
            species = next(csv.reader([line]))[1:]
        elif line.startswith("STATE,"):
            fields = next(csv.reader([line]))
            rows.append([float(value) for value in fields[1:]])
    if species is None or not rows:
        raise ValueError(f"{path}: missing SPECIES or STATE rows")
    return species, np.asarray(rows)


def cantera_const_v_rhs(
    gas: ct.Solution, temperature: float, pressure: float
) -> tuple[float, float, np.ndarray]:
    gas.TPX = temperature, pressure, STATIC_COMPOSITION
    density = float(gas.density)
    production = gas.net_production_rates
    dydt = production * gas.molecular_weights / density
    dtdt = -float(np.dot(gas.partial_molar_int_energies, production)) / (
        density * float(gas.cv_mass)
    )
    return density, dtdt, dydt


def compare_state_grid(gas: ct.Solution, species: list[str], rows: np.ndarray) -> list[str]:
    failures: list[str] = []
    indices = np.asarray([gas.species_index(name) for name in species], dtype=int)
    max_rho_rel = 0.0
    max_t_rel = 0.0
    max_species_l2 = 0.0
    max_species_scaled = 0.0
    worst_case = -1

    expected_columns = 4 + 1 + len(species)
    if rows.shape[1] != expected_columns:
        return [
            f"state-grid column count {rows.shape[1]} != expected {expected_columns}"
        ]

    for row in rows:
        case = int(row[0])
        temperature, pressure, speed_rho = row[1:4]
        speed_dtdt = row[4]
        speed_dydt = row[5:]
        ref_rho, ref_dtdt, ref_dydt_all = cantera_const_v_rhs(
            gas, temperature, pressure
        )
        ref_dydt = ref_dydt_all[indices]

        rho_rel = abs(speed_rho / ref_rho - 1.0)
        t_rel = abs(speed_dtdt - ref_dtdt) / max(
            abs(ref_dtdt), abs(speed_dtdt), 1.0e-30
        )
        species_l2 = float(
            np.linalg.norm(speed_dydt - ref_dydt)
            / max(np.linalg.norm(ref_dydt), np.finfo(float).tiny)
        )
        floor = max(float(np.max(np.abs(ref_dydt))) * 1.0e-10, 1.0e-30)
        species_scaled = float(
            np.max(
                np.abs(speed_dydt - ref_dydt)
                / np.maximum(np.maximum(np.abs(speed_dydt), np.abs(ref_dydt)), floor)
            )
        )

        if max(rho_rel, t_rel, species_l2, species_scaled) > max(
            max_rho_rel, max_t_rel, max_species_l2, max_species_scaled
        ):
            worst_case = case
        max_rho_rel = max(max_rho_rel, rho_rel)
        max_t_rel = max(max_t_rel, t_rel)
        max_species_l2 = max(max_species_l2, species_l2)
        max_species_scaled = max(max_species_scaled, species_scaled)

    print(
        f"STATE_GRID cases={len(rows)} rho_max_rel={max_rho_rel:.3e} "
        f"dTdt_max_rel={max_t_rel:.3e} species_max_rel_l2={max_species_l2:.3e} "
        f"species_max_scaled={max_species_scaled:.3e} worst_case={worst_case}"
    )
    if max_rho_rel > 2.0e-4:
        failures.append(f"state-grid density error {max_rho_rel:.3e} > 2e-4")
    if max_t_rel > 1.0e-2:
        failures.append(f"state-grid dT/dt error {max_t_rel:.3e} > 1e-2")
    if max_species_l2 > 1.0e-2:
        failures.append(f"state-grid species relative L2 {max_species_l2:.3e} > 1e-2")
    if max_species_scaled > 1.0e-1:
        failures.append(
            f"state-grid species maximum scaled error {max_species_scaled:.3e} > 0.1"
        )
    return failures


def read_history(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
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
    if i == 0:
        return float(history[0, 0])
    return float(
        np.interp(threshold, history[i - 1 : i + 1, 1], history[i - 1 : i + 1, 0])
    )


def max_dtdt(history: np.ndarray) -> float:
    return float(np.max(np.diff(history[:, 1]) / np.diff(history[:, 0])))


def cantera_history(mechanism: Path, times: np.ndarray, t0: float, p0: float) -> np.ndarray:
    gas = ct.Solution(str(mechanism))
    gas.TPX = t0, p0, HISTORY_COMPOSITION
    reactor = ct.IdealGasReactor(gas, energy="on", clone=True)
    network = ct.ReactorNet([reactor])
    network.rtol = 1.0e-10
    network.atol = 1.0e-18
    indices = [reactor.phase.species_index(name) for name in HISTORY_SPECIES]
    rows: list[list[float]] = []
    for time in times:
        if time > network.time:
            network.advance(float(time))
        rows.append(
            [time, reactor.T, reactor.phase.P, *reactor.phase.X[indices], np.sum(reactor.phase.Y)]
        )
    return np.asarray(rows)


def compare_history(
    mechanism: Path, numeric: np.ndarray, analytic: np.ndarray | None
) -> list[str]:
    failures: list[str] = []
    reference = cantera_history(
        mechanism, numeric[:, 0], float(numeric[0, 1]), float(numeric[0, 2])
    )

    speed_idt = idt50(numeric)
    cantera_idt = idt50(reference)
    idt_rel = abs(speed_idt / cantera_idt - 1.0) if np.isfinite(speed_idt) and np.isfinite(cantera_idt) else float("nan")
    speed_peak = max_dtdt(numeric)
    cantera_peak = max_dtdt(reference)
    peak_rel = abs(speed_peak / cantera_peak - 1.0)
    t_rel_l2 = float(
        np.linalg.norm(numeric[:, 1] - reference[:, 1])
        / max(np.linalg.norm(reference[:, 1] - reference[0, 1]), np.finfo(float).tiny)
    )
    final_t_error = abs(numeric[-1, 1] - reference[-1, 1])

    print(
        "CANTERA_HISTORY idt50_speed_s={:.9e} idt50_ref_s={:.9e} rel={:.3e}".format(
            speed_idt, cantera_idt, idt_rel
        )
    )
    print(
        "CANTERA_HISTORY max_dTdt_speed={:.9e} max_dTdt_ref={:.9e} rel={:.3e}".format(
            speed_peak, cantera_peak, peak_rel
        )
    )
    print(
        f"CANTERA_HISTORY T_rel_l2={t_rel_l2:.3e} T_final_abs_K={final_t_error:.3e}"
    )

    if not np.isfinite(idt_rel) or idt_rel > 2.0e-2:
        failures.append(f"+50 K IDT relative error {idt_rel:.3e} > 0.02")
    if not np.isfinite(peak_rel) or peak_rel > 3.0e-1:
        failures.append(f"max dT/dt relative error {peak_rel:.3e} > 0.30")
    if t_rel_l2 > 5.0e-2:
        failures.append(f"temperature-history relative L2 {t_rel_l2:.3e} > 0.05")
    if final_t_error > 3.0:
        failures.append(f"final-temperature error {final_t_error:.3e} K > 3 K")

    for column, name in enumerate(HISTORY_SPECIES, start=3):
        rel_l2 = float(
            np.linalg.norm(numeric[:, column] - reference[:, column])
            / max(np.linalg.norm(reference[:, column]), np.finfo(float).tiny)
        )
        final_abs = abs(numeric[-1, column] - reference[-1, column])
        print(
            f"CANTERA_HISTORY species={name} rel_l2={rel_l2:.3e} final_abs={final_abs:.3e}"
        )
        if rel_l2 > 2.0e-1:
            failures.append(f"{name} history relative L2 {rel_l2:.3e} > 0.2")
        if final_abs > 2.0e-4:
            failures.append(f"{name} final absolute error {final_abs:.3e} > 2e-4")

    if analytic is not None:
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
                np.abs(analytic[:, 2] - interpolated[:, 2])
                / np.maximum(np.abs(interpolated[:, 2]), np.finfo(float).tiny)
            )
        )
        species_max = float(
            np.max(np.abs(analytic[:, 3:10] - interpolated[:, 3:10]))
        )
        print(
            f"SOLVERS T_max_abs_K={t_max:.3e} P_max_rel={p_rel:.3e} "
            f"major_species_max_abs={species_max:.3e}"
        )
        if t_max > 0.1 or p_rel > 1.0e-5 or species_max > 1.0e-6:
            failures.append("numeric/analytic history comparison exceeded tolerance")

    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", required=True, type=Path)
    parser.add_argument("--published", required=True, type=Path)
    parser.add_argument("--rates", required=True, type=Path)
    parser.add_argument("--states", required=True, type=Path)
    parser.add_argument("--numeric", required=True, type=Path)
    parser.add_argument("--analytic", type=Path)
    args = parser.parse_args()

    gas = ct.Solution(str(args.mechanism))
    published = ct.Solution(str(args.published))
    failures: list[str] = []

    gas_plog = plog_indices(gas)
    published_plog = plog_indices(published)
    print(
        f"MECHANISM cantera={ct.__version__} species={gas.n_species} "
        f"reactions={gas.n_reactions} plog={len(gas_plog)}"
    )
    if gas.species_names != published.species_names:
        failures.append("regenerated and published C3Mech species order differs")
    if [r.equation for r in gas.reactions()] != [r.equation for r in published.reactions()]:
        failures.append("regenerated and published C3Mech reaction order differs")
    if gas_plog != published_plog or len(gas_plog) != 1:
        failures.append(
            f"unexpected PLOG topology regenerated={gas_plog} published={published_plog}"
        )

    failures.extend(compare_rate_grid(gas, published, read_rate_rows(args.rates)))
    species, state_rows = read_state_grid(args.states)
    failures.extend(compare_state_grid(gas, species, state_rows))
    numeric = read_history(args.numeric)
    analytic = read_history(args.analytic) if args.analytic else None
    failures.extend(compare_history(args.mechanism, numeric, analytic))

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1
    print("RESULT: PASS - public C3Mech PLOG mechanism agrees with Cantera")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
