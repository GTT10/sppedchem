#!/usr/bin/env python3
"""Independent checks for the pinned public C3Mech PLOG mechanism.

The Cantera comparison is deliberately scoped to the PLOG rate coefficient and
its rate of progress. Complete-mechanism RHS and ignition-history differences
are reported as diagnostics because they also include every pre-existing
non-PLOG kinetic and thermodynamic code path in SpeedCHEM. The full public
mechanism is still required to remain finite, to integrate successfully, and
to give matching numerical- and analytic-Jacobian SpeedCHEM histories.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable

import cantera as ct
import numpy as np

HISTORY_SPECIES = ("H2", "O2", "H", "OH", "HO2", "H2O", "AR")
STATIC_COMPOSITION = {
    "H2": 2.0,
    "O2": 1.0,
    "AR": 7.0,
    "H": 1.0e-5,
    "HO2": 1.0e-8,
}
HISTORY_COMPOSITION = {"H2": 2.0, "O2": 1.0, "AR": 7.0}
EXPECTED_RATE_ROWS = 85
EXPECTED_STATE_ROWS = 33


def is_plog_reaction(reaction: ct.Reaction) -> bool:
    return (
        getattr(reaction, "reaction_type", "")
        == "pressure-dependent-Arrhenius"
        or type(reaction.rate).__name__ == "PlogRate"
    )


def plog_indices(gas: ct.Solution) -> list[int]:
    return [
        index
        for index, reaction in enumerate(gas.reactions())
        if is_plog_reaction(reaction)
    ]


def explicit_third_body(reaction: ct.Reaction):
    try:
        return reaction.third_body
    except (AttributeError, ct.CanteraError):
        return None


def forward_rate_order(reaction: ct.Reaction) -> tuple[float, str, bool]:
    """Return the dimensional forward order used by the rate coefficient.

    Cantera removes an automatically detected explicit collider from the
    reactant stoichiometry and exposes it as ``reaction.third_body``. When its
    ``mass_action`` flag is true, that collider still contributes one power of
    concentration to the dimensional order of k.
    """

    orders = dict(getattr(reaction, "orders", {}) or {})
    order = (
        float(sum(orders.values()))
        if orders
        else float(sum(reaction.reactants.values()))
    )

    third_body = explicit_third_body(reaction)
    collider = ""
    mass_action = False
    if third_body is not None:
        collider = str(getattr(third_body, "name", "") or "")
        mass_action = bool(getattr(third_body, "mass_action", False))
        if mass_action and collider not in reaction.reactants:
            order += 1.0

    return order, collider, mass_action


def read_rate_rows(path: Path) -> list[tuple[int, int, float, float, float]]:
    rows: list[tuple[int, int, float, float, float]] = []
    for line in path.read_text().splitlines():
        if not line.startswith("RATE,"):
            continue
        fields = next(csv.reader([line]))
        rows.append(
            (
                int(fields[1]),
                int(fields[2]),
                float(fields[3]),
                float(fields[4]),
                float(fields[5]),
            )
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
    rows = list(rows)
    if len(rows) != EXPECTED_RATE_ROWS:
        failures.append(
            f"PLOG rate-grid row count {len(rows)} != {EXPECTED_RATE_ROWS}"
        )

    known_plog = set(plog_indices(gas))
    published_plog = set(plog_indices(published))
    max_speed_rel = 0.0
    max_published_rel = 0.0
    worst: tuple[int, float, float] | None = None
    metadata_reported: set[int] = set()
    count = 0

    for packed, reaction_one_based, temperature, pressure, speed_k_ck in rows:
        del packed
        reaction_index = reaction_one_based - 1
        if reaction_index not in known_plog:
            failures.append(
                f"SpeedCHEM marked reaction {reaction_one_based} as PLOG, "
                "but Cantera did not"
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

        order, collider, mass_action = forward_rate_order(reaction)
        if reaction_index not in metadata_reported:
            print(
                f"PLOG_METADATA reaction={reaction_one_based} "
                f"order={order:.1f} collider={collider or 'none'} "
                f"mass_action={mass_action}"
            )
            metadata_reported.add(reaction_index)

        # SpeedCHEM/CHEMKIN uses mol/cm^3; Cantera uses kmol/m^3.
        # Since 1 mol/cm^3 = 1000 kmol/m^3,
        #     k_SI = k_CK * 1000**(1-order).
        speed_k_si = speed_k_ck * 1000.0 ** (1.0 - order)

        gas.TP = temperature, pressure
        published.TP = temperature, pressure
        cantera_k = float(gas.forward_rate_constants[reaction_index])
        published_k = float(published.forward_rate_constants[reaction_index])

        speed_rel = abs(speed_k_si - cantera_k) / max(
            abs(cantera_k), np.finfo(float).tiny
        )
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
        failures.append(
            f"PLOG rate-grid relative error {max_speed_rel:.3e} > 5e-4"
        )
    if max_published_rel > 5.0e-12:
        failures.append(
            "published/regenerated Cantera PLOG mismatch "
            f"{max_published_rel:.3e} > 5e-12"
        )
    return failures


def _debug_value(lines: list[str], prefix: str) -> float:
    for line in lines:
        if line.startswith(prefix):
            return float(line.split("=", 1)[1].split()[0])
    raise ValueError(f"debug log: missing {prefix!r}")


def read_debug_rop(path: Path) -> tuple[int, float, float, float, float, float]:
    lines = path.read_text().splitlines()
    reaction = int(_debug_value(lines, "DEBUG reaction="))
    temperature = _debug_value(lines, "DEBUG T=")
    pressure = _debug_value(lines, "DEBUG P=")
    qf = _debug_value(lines, "DEBUG qf=")
    qb = _debug_value(lines, "DEBUG qb=")
    qnet = _debug_value(lines, "DEBUG q=")
    return reaction, temperature, pressure, qf, qb, qnet


def compare_rate_of_progress(
    gas: ct.Solution,
    debug: tuple[int, float, float, float, float, float],
) -> list[str]:
    failures: list[str] = []
    reaction_one_based, temperature, pressure, qf_ck, qb_ck, qnet_ck = debug
    reaction_index = reaction_one_based - 1
    if reaction_index not in plog_indices(gas):
        return [f"debug reaction {reaction_one_based} is not PLOG in Cantera"]

    if abs(qnet_ck - (qf_ck - qb_ck)) > 1.0e-12 * max(
        abs(qf_ck), abs(qb_ck), 1.0
    ):
        failures.append("SpeedCHEM PLOG q != qf-qb")

    gas.TPX = temperature, pressure, STATIC_COMPOSITION
    # SpeedCHEM progress rates are mol/cm^3/s; Cantera uses kmol/m^3/s.
    speed = np.asarray([qf_ck, qb_ck, qnet_ck]) * 1000.0
    reference = np.asarray(
        [
            gas.forward_rates_of_progress[reaction_index],
            gas.reverse_rates_of_progress[reaction_index],
            gas.net_rates_of_progress[reaction_index],
        ],
        dtype=float,
    )

    labels = ("forward", "reverse", "net")
    tolerances = (5.0e-4, 5.0e-3, 5.0e-4)
    for label, speed_value, reference_value, tolerance in zip(
        labels, speed, reference, tolerances, strict=True
    ):
        relative = abs(speed_value - reference_value) / max(
            abs(speed_value), abs(reference_value), 1.0e-30
        )
        print(
            f"PLOG_ROP kind={label} speedchem={speed_value:.12e} "
            f"cantera={reference_value:.12e} rel={relative:.3e}"
        )
        if relative > tolerance:
            failures.append(
                f"PLOG {label} rate-of-progress error {relative:.3e} "
                f"> {tolerance:.1e}"
            )

    return failures


def read_state_grid(path: Path) -> tuple[list[str], np.ndarray]:
    species: list[str] | None = None
    rows: list[list[float]] = []
    nonfinite_rows = 0
    for line in path.read_text().splitlines():
        if line.startswith("SPECIES,"):
            species = next(csv.reader([line]))[1:]
        elif line.startswith("STATE,"):
            fields = next(csv.reader([line]))
            rows.append([float(value) for value in fields[1:]])
        elif line.startswith("NONFINITE,"):
            nonfinite_rows += 1
    if species is None or not rows:
        raise ValueError(f"{path}: missing SPECIES or STATE rows")
    if nonfinite_rows:
        raise ValueError(f"{path}: contains {nonfinite_rows} non-finite rows")
    return species, np.asarray(rows, dtype=float)


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


def diagnose_state_grid(
    gas: ct.Solution, species: list[str], rows: np.ndarray
) -> list[str]:
    failures: list[str] = []
    if len(rows) != EXPECTED_STATE_ROWS:
        failures.append(
            f"real-mechanism state-grid row count {len(rows)} "
            f"!= {EXPECTED_STATE_ROWS}"
        )
    if not np.all(np.isfinite(rows)):
        failures.append("real-mechanism state grid contains non-finite values")
        return failures

    indices = np.asarray([gas.species_index(name) for name in species], dtype=int)
    expected_columns = 5 + len(species)
    if rows.shape[1] != expected_columns:
        return [
            f"state-grid column count {rows.shape[1]} != expected {expected_columns}"
        ]

    max_rho_rel = 0.0
    max_t_rel = 0.0
    max_species_l2 = 0.0
    max_species_scaled = 0.0
    worst_case = -1
    worst_metric = -1.0

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
                / np.maximum(
                    np.maximum(np.abs(speed_dydt), np.abs(ref_dydt)), floor
                )
            )
        )
        metric = max(rho_rel, t_rel, species_l2, species_scaled)
        if metric > worst_metric:
            worst_metric = metric
            worst_case = case

        max_rho_rel = max(max_rho_rel, rho_rel)
        max_t_rel = max(max_t_rel, t_rel)
        max_species_l2 = max(max_species_l2, species_l2)
        max_species_scaled = max(max_species_scaled, species_scaled)

    print(
        f"STATE_GRID_DIAGNOSTIC cases={len(rows)} "
        f"rho_max_rel={max_rho_rel:.3e} dTdt_max_rel={max_t_rel:.3e} "
        f"species_max_rel_l2={max_species_l2:.3e} "
        f"species_max_scaled={max_species_scaled:.3e} "
        f"worst_case={worst_case} scope=full-mechanism-nonplog-included"
    )

    # Density is determined only by T/P/composition and therefore remains a
    # valid cross-implementation state-construction check. Full RHS differences
    # are diagnostic because they include all non-PLOG SpeedCHEM paths.
    if max_rho_rel > 2.0e-4:
        failures.append(f"state-grid density error {max_rho_rel:.3e} > 2e-4")
    return failures


def read_history(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        if line.startswith("HIST,"):
            rows.append([float(value) for value in line.split(",")[1:]])
    if len(rows) < 3:
        raise ValueError(f"{path}: expected at least 3 HIST rows")
    history = np.asarray(rows, dtype=float)
    if not np.all(np.isfinite(history)):
        raise ValueError(f"{path}: non-finite history values")
    if not np.all(np.diff(history[:, 0]) > 0.0):
        raise ValueError(f"{path}: history times are not strictly increasing")
    return history


def idt50(history: np.ndarray) -> float:
    threshold = history[0, 1] + 50.0
    crossing = np.flatnonzero(history[:, 1] >= threshold)
    if not len(crossing):
        return float("nan")
    index = int(crossing[0])
    if index == 0:
        return float(history[0, 0])
    return float(
        np.interp(
            threshold,
            history[index - 1 : index + 1, 1],
            history[index - 1 : index + 1, 0],
        )
    )


def max_dtdt(history: np.ndarray) -> float:
    return float(np.max(np.diff(history[:, 1]) / np.diff(history[:, 0])))


def cantera_history(
    mechanism: Path, times: np.ndarray, temperature0: float, pressure0: float
) -> np.ndarray:
    gas = ct.Solution(str(mechanism))
    gas.TPX = temperature0, pressure0, HISTORY_COMPOSITION
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
            [
                time,
                reactor.T,
                reactor.phase.P,
                *reactor.phase.X[indices],
                np.sum(reactor.phase.Y),
            ]
        )
    return np.asarray(rows, dtype=float)


def compare_solver_histories(
    mechanism: Path, numeric: np.ndarray, analytic: np.ndarray
) -> list[str]:
    failures: list[str] = []
    if numeric.shape[1] != 11 or analytic.shape[1] != 11:
        return [
            f"unexpected history widths numeric={numeric.shape[1]} "
            f"analytic={analytic.shape[1]}"
        ]

    numeric_at_analytic = np.column_stack(
        [
            analytic[:, 0],
            *(
                np.interp(analytic[:, 0], numeric[:, 0], numeric[:, column])
                for column in range(1, numeric.shape[1])
            ),
        ]
    )
    t_max = float(np.max(np.abs(analytic[:, 1] - numeric_at_analytic[:, 1])))
    p_rel = float(
        np.max(
            np.abs(analytic[:, 2] - numeric_at_analytic[:, 2])
            / np.maximum(np.abs(numeric_at_analytic[:, 2]), np.finfo(float).tiny)
        )
    )
    species_max = float(
        np.max(np.abs(analytic[:, 3:10] - numeric_at_analytic[:, 3:10]))
    )
    mass_sum_max = float(
        max(
            np.max(np.abs(numeric[:, 10] - 1.0)),
            np.max(np.abs(analytic[:, 10] - 1.0)),
        )
    )
    numeric_idt = idt50(numeric)
    analytic_idt = idt50(analytic)
    solver_idt_rel = (
        abs(numeric_idt / analytic_idt - 1.0)
        if np.isfinite(numeric_idt)
        and np.isfinite(analytic_idt)
        and analytic_idt > 0.0
        else float("nan")
    )

    print(
        f"SOLVERS T_max_abs_K={t_max:.3e} P_max_rel={p_rel:.3e} "
        f"major_species_max_abs={species_max:.3e} "
        f"idt50_rel={solver_idt_rel:.3e} sumY_max_abs={mass_sum_max:.3e}"
    )
    if (
        t_max > 0.1
        or p_rel > 1.0e-5
        or species_max > 1.0e-6
        or not np.isfinite(solver_idt_rel)
        or solver_idt_rel > 1.0e-4
        or mass_sum_max > 1.0e-10
    ):
        failures.append("numeric/analytic real-mechanism history exceeded tolerance")

    reference = cantera_history(
        mechanism,
        numeric[:, 0],
        float(numeric[0, 1]),
        float(numeric[0, 2]),
    )
    speed_idt = idt50(numeric)
    cantera_idt = idt50(reference)
    idt_rel = (
        abs(speed_idt / cantera_idt - 1.0)
        if np.isfinite(speed_idt)
        and np.isfinite(cantera_idt)
        and cantera_idt > 0.0
        else float("nan")
    )
    speed_peak = max_dtdt(numeric)
    cantera_peak = max_dtdt(reference)
    peak_rel = abs(speed_peak / cantera_peak - 1.0)
    t_rel_l2 = float(
        np.linalg.norm(numeric[:, 1] - reference[:, 1])
        / max(
            np.linalg.norm(reference[:, 1] - reference[0, 1]),
            np.finfo(float).tiny,
        )
    )
    final_t_error = abs(numeric[-1, 1] - reference[-1, 1])

    print(
        "CANTERA_HISTORY_DIAGNOSTIC "
        "idt50_speed_s={:.9e} idt50_ref_s={:.9e} rel={:.3e} "
        "max_dTdt_rel={:.3e} T_rel_l2={:.3e} T_final_abs_K={:.3e} "
        "scope=full-mechanism-nonplog-included".format(
            speed_idt,
            cantera_idt,
            idt_rel,
            peak_rel,
            t_rel_l2,
            final_t_error,
        )
    )
    for column, name in enumerate(HISTORY_SPECIES, start=3):
        rel_l2 = float(
            np.linalg.norm(numeric[:, column] - reference[:, column])
            / max(np.linalg.norm(reference[:, column]), np.finfo(float).tiny)
        )
        final_abs = abs(numeric[-1, column] - reference[-1, column])
        print(
            f"CANTERA_HISTORY_DIAGNOSTIC species={name} "
            f"rel_l2={rel_l2:.3e} final_abs={final_abs:.3e}"
        )

    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mechanism", required=True, type=Path)
    parser.add_argument("--published", required=True, type=Path)
    parser.add_argument("--rates", required=True, type=Path)
    parser.add_argument("--debug", required=True, type=Path)
    parser.add_argument("--states", required=True, type=Path)
    parser.add_argument("--numeric", required=True, type=Path)
    parser.add_argument("--analytic", required=True, type=Path)
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
    if [reaction.equation for reaction in gas.reactions()] != [
        reaction.equation for reaction in published.reactions()
    ]:
        failures.append("regenerated and published C3Mech reaction order differs")
    if gas_plog != published_plog or len(gas_plog) != 1:
        failures.append(
            f"unexpected PLOG topology regenerated={gas_plog} "
            f"published={published_plog}"
        )
    elif len(gas.reaction(gas_plog[0]).rate.rates) != 11:
        failures.append("public C3Mech PLOG reaction does not contain 11 terms")

    failures.extend(compare_rate_grid(gas, published, read_rate_rows(args.rates)))
    failures.extend(compare_rate_of_progress(gas, read_debug_rop(args.debug)))
    species, state_rows = read_state_grid(args.states)
    failures.extend(diagnose_state_grid(gas, species, state_rows))
    numeric = read_history(args.numeric)
    analytic = read_history(args.analytic)
    failures.extend(compare_solver_histories(args.mechanism, numeric, analytic))

    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1

    print(
        "RESULT: PASS - public C3Mech PLOG rates/ROP and real-mechanism "
        "solver paths verified"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
