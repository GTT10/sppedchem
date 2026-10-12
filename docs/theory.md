# Governing equations and numerical solvers

For a shorter, informal introduction in Japanese, see
[SpeedCHEMの中で何が起きているか](guide-ja.md).

SpeedCHEM integrates the temperature and composition of a homogeneous,
ideal-gas mixture undergoing chemical reactions at fixed density. This
document derives that initial-value problem and explains how the numerical
solvers and analytical Jacobian act on it. The companion
[PLOG document](plog.md) develops the pressure-dependent rate calculation
and its derivatives in detail.

The equations below describe the constant-volume path implemented in this
repository. Source links identify the routines that implement each part;
the papers in [References](#references) explain the original SpeedCHEM
methods. Solver interfaces present in the source are described here without
asserting that every solver has been qualified for every mechanism.

## Contents

- [1. Model and state variables](#1-model-and-state-variables)
- [2. Reaction rates and species equations](#2-reaction-rates-and-species-equations)
- [3. Thermodynamics and the temperature equation](#3-thermodynamics-and-the-temperature-equation)
- [4. Why the ODE is stiff](#4-why-the-ode-is-stiff)
- [5. Implicit integration and solver choices](#5-implicit-integration-and-solver-choices)
- [6. Structure of the analytical Jacobian](#6-structure-of-the-analytical-jacobian)
- [7. Evaluation sequence and implementation map](#7-evaluation-sequence-and-implementation-map)

## 1. Model and state variables

The reactor is closed and adiabatic, with fixed volume and no transport
terms. Its mass and volume are constant, so its density is constant during
one chemistry integration interval. In a CFD application, the host code
supplies the initial state and density for that interval and handles flow,
diffusion, spray sources, and other processes separately.

For $N_s$ species, the public integration state is

$$
\boldsymbol y=(T,Y_1,\ldots,Y_{N_s})^{\mathsf T},
\qquad
\frac{d\boldsymbol y}{dt}=\boldsymbol f(\boldsymbol y;\rho).
$$

Here $T$ is temperature, $Y_i$ is a species mass fraction, and $\rho$ is a
parameter held fixed during integration. `chemistry_ODE_integrate` requires
`neq = ns + 1`. All species mass fractions are integrated; there is no
dependent last species eliminated from the state. Internal permutations
can move temperature to the last position and reorder species to reduce
linear-system fill-in. The public input/output order remains
`[T, Y_1, ..., Y_ns]`.

We use SI units in the derivation:

| Symbol | Meaning | Units |
| --- | --- | --- |
| $T$ | Temperature | K |
| $\rho$ | Mixture density | kg/m³ |
| $Y_i$ | Species mass fraction | dimensionless |
| $W_i$ | Species molecular weight | kg/mol |
| $C_i=\rho Y_i/W_i$ | Species molar concentration | mol/m³ |
| $q_r$ | Net rate of progress of reaction $r$ | mol/(m³ s) |
| $\dot\omega_i$ | Species molar production rate | mol/(m³ s) |
| $\bar u_i,\bar h_i$ | Molar internal energy and enthalpy | J/mol |
| $c_v$ | Mixture mass-specific heat capacity | J/(kg K) |

Define $S=\sum_i Y_i/W_i$. The ideal-gas equation gives

$$
\bar W=\frac{1}{S},
\qquad P=\rho R_u T S,
\qquad R_u=8.31446261815324\;\mathrm{J/(mol\,K)}.
$$

The expression for $\bar W$ assumes $\sum_iY_i=1$. The code evaluates
pressure directly from $\sum_iY_i/W_i$, without normalizing the state.
Pressure evolves as temperature and composition evolve, even though
density is fixed. It is an algebraic quantity rather than an additional
ODE unknown. See `pressurerhoT` in [SCmodule.f90](../src/SCmodule.f90).

### Units used inside the implementation

The rate routines use CHEMKIN-style mol–cm–s units, while density and
mixture heat capacities use SI. `SCMW` is numerically in g/mol
(equivalently kg/kmol). Thus

$$
W_i=10^{-3}\,\mathrm{SCMW}_i,
\qquad
C_{i,\mathrm{code}}=10^{-6}C_i
=10^{-3}\rho\frac{Y_i}{\mathrm{SCMW}_i}.
$$

`SCdwdt` returns mol/(cm³ s), so the implementation converts its result as

$$
\frac{dY_i}{dt}
=\frac{10^3\,\mathrm{SCMW}_i}{\rho}
\dot\omega_{i,\mathrm{code}},
\qquad
\frac{dT}{dt}
=-\frac{10^6}{\rho c_v}
\sum_i\bar u_i\dot\omega_{i,\mathrm{code}}.
$$

These factors account for the `milli`, `kilo`, and `urho` operations in
[SC_conV](../src/SCconV.f90). In particular, the rate-constant units must
match the concentration units used in the mass-action products.

## 2. Reaction rates and species equations

Write reaction $r$ as

$$
\sum_i\nu'_{ir}X_i\rightleftharpoons\sum_i\nu''_{ir}X_i,
\qquad \nu_{ir}=\nu''_{ir}-\nu'_{ir}.
$$

For ordinary mass-action kinetics,

$$
q_{f,r}=k_{f,r}\prod_i C_i^{\nu'_{ir}},
\qquad
q_{b,r}=k_{b,r}\prod_i C_i^{\nu''_{ir}},
\qquad q_r=q_{f,r}-q_{b,r}.
$$

For an irreversible reaction, $q_{b,r}=0$. This runtime path uses
stoichiometric mass-action orders. Although the interpreter recognizes
modified-order `FORD/RORD` syntax, the SpeedCHEM binary loader rejects
it as unsupported. The production rates use the stoichiometric
difference:

$$
\dot\omega_i=\sum_r\nu_{ir}q_r,
\qquad
\boxed{\frac{dY_i}{dt}=\frac{W_i}{\rho}\dot\omega_i}.
$$

The sparse stoichiometric matrices perform this final reaction-to-species
mapping. For a balanced mechanism,
$\sum_iW_i\nu_{ir}=0$, so $\sum_i dY_i/dt=0$. Elemental abundances are
conserved for the same reason. Numerical integration introduces finite
errors; conservation alone does not guarantee positive mass fractions.

### Arrhenius and reverse rate constants

The modified Arrhenius expression is

$$
k_f(T)=A T^b\exp\left(-\frac{E_a}{R_uT}\right),
\qquad
\frac{d\ln k_f}{dT}=\frac bT+\frac{E_a/R_u}{T^2}.
$$

The interpreter converts the activation-energy units specified by the
mechanism into an $E_a/R_u$ value in kelvin. The numerical value and units
of $A$ depend on concentration units, reaction order, and the convention
of inserting temperature in kelvin.

For a reversible reaction without explicit reverse Arrhenius data,
thermodynamics supplies the reverse rate:

$$
\Delta_r\bar g^\circ=\sum_i\nu_{ir}
\left(\bar h_i^\circ-T\bar s_i^\circ\right),
\qquad
K_{p,r}=\exp\left(-\frac{\Delta_r\bar g^\circ}{R_uT}\right),
$$

$$
\Delta\nu_r=\sum_i\nu_{ir},
\qquad
K_{c,r}=K_{p,r}\left(\frac{p^\circ}{R_uT}\right)^{\Delta\nu_r},
\qquad k_{b,r}=\frac{k_{f,r}}{K_{c,r}}.
$$

$K_p$ uses dimensionless ideal-gas activities $p_i/p^\circ$. $K_c$ here
is the concentration-form constant with units consistent with the
mass-action products. The implementation uses $p^\circ=101325$ Pa and
includes the concentration-volume conversion for mol/cm³ in `uequilC`.
Its stored quantity is $1/K_c$. General reactions may instead provide
explicit reverse Arrhenius coefficients; PLOG combined with explicit
`REV` is rejected by this version.

### Third-body and falloff models

An effective collider concentration is

$$
[M]_r=\sum_i\alpha_{ir}C_i,
$$

where $\alpha_{ir}$ is a species efficiency. A simple third-body
reaction multiplies its base rate coefficient by $[M]_r$.
Lindemann falloff blends low- and high-pressure coefficients:

$$
P_{r,\mathrm{red}}=\frac{k_{0,r}[M]_r}{k_{\infty,r}},
\qquad
k_{f,r}=k_{\infty,r}
\frac{P_{r,\mathrm{red}}}{1+P_{r,\mathrm{red}}}F_r.
$$

Lindemann uses $F_r=1$; Troe supplies a temperature- and
reduced-pressure-dependent broadening factor. The source variable `Pr`
in the constant-volume Jacobian stores the bounded fraction
$P_{r,\mathrm{red}}/(1+P_{r,\mathrm{red}})$, rather than the conventional
reduced pressure itself. This distinction matters when reading the
derivatives in [SCconV.f90](../src/SCconV.f90).

PLOG instead specifies a table of pressure-dependent coefficients, as
derived in [PLOG rates and their derivatives](plog.md). It does not add
another third-body multiplier in this implementation.

## 3. Thermodynamics and the temperature equation

### NASA-7 species properties

Each temperature range uses seven coefficients $a_1,\ldots,a_7$:

$$
\frac{\bar c_{p,i}^\circ}{R_u}
=a_1+a_2T+a_3T^2+a_4T^3+a_5T^4,
$$

$$
\frac{\bar h_i^\circ}{R_uT}
=a_1+\frac{a_2}{2}T+\frac{a_3}{3}T^2
+\frac{a_4}{4}T^3+\frac{a_5}{5}T^4+\frac{a_6}{T},
$$

$$
\frac{\bar s_i^\circ}{R_u}
=a_1\ln T+a_2T+\frac{a_3}{2}T^2
+\frac{a_4}{3}T^3+\frac{a_5}{4}T^4+a_7.
$$

These are the conventional
[NASA-7 property relations](https://cantera.org/stable/reference/thermo/species-thermo.html#the-nasa-7-coefficient-polynomial-parameterization),
written with one-based coefficient indices. For an ideal gas,

$$
\bar u_i=\bar h_i-R_uT,
\qquad \bar c_{v,i}=\bar c_{p,i}-R_u,
\qquad
c_v(T,\boldsymbol Y)=\sum_iY_i\frac{\bar c_{v,i}(T)}{W_i}.
$$

Species properties use their own low/high temperature switch. In direct
evaluation, the high-temperature coefficients are selected at equality
with the switch temperature. The thermo file supplies the formation
enthalpy information as well as sensible enthalpy; changing it changes
both energy release and equilibrium-derived reverse rates.

### Energy conservation

The mixture mass-specific internal energy is

$$
u(T,\boldsymbol Y)=\sum_iY_i\frac{\bar u_i(T)}{W_i}.
$$

For a closed, adiabatic, fixed-volume reactor, $du/dt=0$. Differentiating
the mixture energy gives

$$
0=c_v\frac{dT}{dt}
+\sum_i\frac{\bar u_i}{W_i}\frac{dY_i}{dt}.
$$

Substituting the species equations yields

$$
\boxed{\frac{dT}{dt}
=-\frac{1}{\rho c_v}\sum_i\bar u_i\dot\omega_i}.
$$

The reaction energy is already represented by the change of composition
in the species internal energies. Adding a second chemical heat-release
term to this equation would count that energy twice. Constant-pressure
reactors use an enthalpy balance and a different density evolution;
their temperature equation is a different model.

### Direct evaluation and tabulation

`accurate_scthermo = .true.` selects direct evaluation of
temperature-dependent properties. With tabulation enabled, the library
precomputes thermodynamic, equilibrium, and ordinary Arrhenius data and
interpolates in the allowed interior of the temperature table. Near or
outside that interior, `use_table` selects direct evaluation. These
temperature bounds control evaluation mode; they do not clamp the
integrated temperature or establish physical validity outside a
species' fitted range.

Tabulation introduces an additional approximation beyond the ODE error
control. PLOG terms are evaluated directly by `plog_kinf_eval`, including
their pressure and temperature derivatives, even when the legacy
temperature tables are enabled.

## 4. Why the ODE is stiff

Combustion mechanisms contain processes with widely separated time
scales. Fast radical relaxation can occur while fuel consumption or heat
release evolves much more slowly. Locally, the Jacobian eigenvalues
describe these coupled modes.

For the scalar decay $z'=-\lambda z$, explicit Euler gives
$z_{n+1}=(1-h\lambda)z_n$ and requires $h\lambda\leq2$ for stability
when $\lambda>0$. A very fast, stable chemical mode can therefore force
tiny explicit steps even when the observable of interest changes slowly.
An implicit method treats those fast modes through a coupled solve.
It still needs small steps when accuracy demands them, for example
during rapid ignition.

## 5. Implicit integration and solver choices

### BDF and Newton iteration

A backward differentiation formula (BDF) step can be written in normalized
form as

$$
\boldsymbol y_{n+1}-\boldsymbol a_n
-h\beta_n\boldsymbol f(\boldsymbol y_{n+1})=0,
$$

where $\boldsymbol a_n$ contains previous states and $\beta_n$ depends on
order and step history. Newton-type corrections solve

$$
\left(I-h\beta_nJ\right)\delta\boldsymbol y
=-\boldsymbol R,
\qquad J=\frac{\partial\boldsymbol f}{\partial\boldsymbol y}.
$$

The integrator can reuse a Jacobian or factorization over several
iterations or steps. Its actual update strategy belongs to the selected
solver. A sparse analytical Jacobian reduces derivative-evaluation cost
and lets a sparse linear solver exploit the reaction network's structure.
The cost also depends on matrix fill-in and ordering; a sparse input
matrix does not imply that its factors are equally sparse.

### Implicit Runge–Kutta and Rosenbrock methods

An implicit Runge–Kutta method defines coupled stage derivatives

$$
\boldsymbol k_i=\boldsymbol f\left(
\boldsymbol y_n+h\sum_{j=1}^{s}a_{ij}\boldsymbol k_j\right),
\qquad
\boldsymbol y_{n+1}=\boldsymbol y_n+h\sum_{i=1}^{s}b_i\boldsymbol k_i.
$$

RADAU5 uses the three-stage Radau IIA method, with order five. The
nonlinear stage equations require iteration and Jacobian-based linear
solves. Its collocation construction provides strong damping of fast
decaying modes. Compared with a multistep BDF method, the method organizes
the implicit work into stages within each step rather than a formula
using several previously accepted states.

For an autonomous system, a Rosenbrock-type method instead has linear
stage equations of the schematic form

$$
\left(I-h\gamma J_n\right)\boldsymbol k_i
=\boldsymbol f\left(\boldsymbol y_n
+h\sum_{j<i}a_{ij}\boldsymbol k_j\right)
+hJ_n\sum_{j<i}\gamma_{ij}\boldsymbol k_j.
$$

Here the stage derivatives have units of $\boldsymbol y/t$; individual
implementations may scale their stages differently. The coefficients
belong to the selected method. The key operation is a linear solve
with a Jacobian-based matrix, avoiding a nonlinear solve for each stage.
RODAS follows this family. ROWMAP combines a Rosenbrock–W formulation
with Krylov subspaces, using Jacobian-vector information to carry out
the stage computation. These methods integrate the same reactor RHS;
they change its numerical treatment.

### Solver interfaces in this repository

Selection and workspace allocation are implemented in
[chemistry_input.f90](../src/chemistry_input.f90). The following names are
the strings accepted by that dispatch:

| Selection | Time-integration family | Jacobian path |
| --- | --- | --- |
| `LSODES`, `LSODESJAC` | BDF with general sparse direct linear algebra | Internal differences / supplied sparse Jacobian |
| `LSODE`, `LSODEJAC` | BDF | Internal differences / supplied full Jacobian |
| `LSODA`, `LSODAJAC` | Automatic Adams/BDF switching | Internal differences / supplied full Jacobian |
| `VODE`, `VODEJAC` | BDF in the configured mode | Internal differences / supplied full Jacobian |
| `VODES`, `VODESJAC` | Fortran-90 VODE with sparse linear algebra | Internal differences / supplied sparse Jacobian |
| `RADAU5`, `RADAU5JAC` | Fifth-order Radau IIA implicit Runge–Kutta | Internal differences / supplied Jacobian through the repository's sparse adapter |
| `RODAS`, `RODASJAC` | Rosenbrock method | Internal differences / supplied Jacobian |
| `ROWMAP`, `ROWMAPJAC` | Rosenbrock–Krylov method | Internal path / supplied Jacobian-vector product |
| `GAM`, `GAMJAC` | Generalized Adams method | Internal differences / supplied Jacobian |
| `DASPK`, `DASPKJAC` | Implicit residual formulation | Internal differences / supplied residual Jacobian |
| `MEBDFJAC` | Modified extended BDF with sparse direct linear algebra | Supplied Jacobian; `MEBDF` without `JAC` stops as unsupported |

The [ODEPACK overview](https://computing.llnl.gov/projects/odepack)
describes the LSODE family. In this wrapper, `LSODESJAC` uses method flag
21 with a supplied sparsity pattern and Jacobian; `LSODES` uses flag 222
and lets DLSODES determine these internally. The method flag, not merely
the solver family name, identifies its configured mode.

[Hairer's solver descriptions](https://www.unige.ch/~hairer/software.html)
identify RADAU5 as Radau IIA of order five and RODAS as a Rosenbrock
method. Radau solves coupled implicit stage equations. Rosenbrock methods
instead form linear systems using a Jacobian at a step's reference state;
ROWMAP uses a Krylov approach for the stage computation.

Although DASPK accepts differential-algebraic problems, the chemistry
adapter here defines

$$
\boldsymbol F(\boldsymbol y,\dot{\boldsymbol y})
=\dot{\boldsymbol y}-\boldsymbol f(\boldsymbol y)=0.
$$

This is a residual representation of the same ODE, with no extra
algebraic chemical-equilibrium constraint. The analytic DASPK adapter
forms $c_jI-J$, consistent with this residual convention.

Without `itapeChem`, `solver_setup` chooses `LSODESJAC`. The `SC_SOLVER`
environment variable can override the solver string before workspace
allocation, for example `SC_SOLVER=LSODESJAC`. This choice does not
establish an optimal solver for every mechanism or state.

### Tolerances

For the ODEPACK paths configured with a scalar relative tolerance and
vector absolute tolerances, error weights have the form

$$
\mathrm{EWT}_i=\mathrm{ATOL}_i+\mathrm{RTOL}|y_i|.
$$

The solver controls a norm of estimated local error divided by these
weights. Other backends use their own error norms and adapters, and some
receive only `atol(1)`, so the same settings are not interchangeable
across every solver.

The setup defaults are `TOLR = 1e-4`, `YTOLA = 1e-15`, and
`TTOLA = 1e-2`, with `separate_tols = .false.`. Thus `TTOLA` is not
automatically the active temperature tolerance: initialization fills
`atol` with `YTOLA` and replaces the temperature entry with `TTOLA` only
when separate tolerances are enabled. The public integration call takes
its actual tolerance arguments. Temperature and mass fractions have
different units, making that distinction relevant.

Local error control is not a bound on ignition-delay error or agreement
with experiments. Those also depend on the mechanism, thermodynamic
data, reactor assumptions, and any CFD coupling.

## 6. Structure of the analytical Jacobian

With temperature first, the Jacobian has the block form

$$
J=\begin{pmatrix}
J_{TT}&\boldsymbol J_{TY}\\
\boldsymbol J_{YT}&J_{YY}
\end{pmatrix}.
$$

The derivatives are taken at fixed density. Each $Y_j$ is an independent
coordinate of the full solver state; the callback does not renormalize
mass fractions during differentiation. Enforcing $\sum_iY_i=1$ while
perturbing a species would define a different, constrained derivative.

For a plain Arrhenius mass-action reaction with positive concentrations,

$$
\left(\frac{\partial q_r}{\partial Y_j}\right)_{T,\rho}
=\frac{\nu'_{jr}q_{f,r}-\nu''_{jr}q_{b,r}}{Y_j}.
$$

This expression contains only species participating in that reaction's
rate products. Third-body, falloff, and PLOG coefficients add their own
composition derivatives. The formula with division by $Y_j$ is a convenient
positive-state expression; the original product derivative gives the
appropriate limit at zero concentration.

Let $\boldsymbol f_Y$ be the species portion of the RHS and define the
mass-specific species energies $e_i=\bar u_i/W_i$ and heat capacities
$c_{v,i}=\bar c_{v,i}/W_i$. Then

$$
(J_{YY})_{ij}=\frac{W_i}{\rho}
\sum_r\nu_{ir}\frac{\partial q_r}{\partial Y_j},
\qquad
(J_{YT})_i=\frac{W_i}{\rho}
\sum_r\nu_{ir}\frac{\partial q_r}{\partial T}.
$$

Differentiating $f_T=-\boldsymbol e^{\mathsf T}\boldsymbol f_Y/c_v$
also differentiates the energy and heat capacity:

$$
(J_{TY})_j
=-\frac{1}{c_v}
\left[\sum_i e_i(J_{YY})_{ij}+f_Tc_{v,j}\right],
$$

$$
J_{TT}=-\frac{1}{c_v}
\left[\sum_i c_{v,i}(f_Y)_i
+\sum_i e_i(J_{YT})_i
+f_T\frac{\partial c_v}{\partial T}\right].
$$

These are the underlying direct-property relations. The tabulation path
uses the corresponding stored/interpolated properties and derivatives.
The temperature row and column couple many species even when most
reaction products are sparse.

The legacy `simplified_for_sparsity` option omits the default-collider
concentration contribution from third-body composition derivatives,
retaining an approximation with fewer nonzeros. When any PLOG reaction
is present, initialization disables this option. PLOG's exact additional
coupling is derived in
[The constant-volume Jacobian](plog.md#5-the-constant-volume-jacobian).

For very small species values, the RHS/Jacobian routines substitute a
signed small value locally to avoid singular derivative expressions.
These numerical guards do not enforce positivity or project the solution
onto a conserved-composition manifold. At such guards, pressure knots,
or property-switch boundaries, smooth-state derivative formulas need
their branch conventions taken into account.

## 7. Evaluation sequence and implementation map

During one constant-volume RHS evaluation, the code:

1. Obtains $T$ and $\boldsymbol Y$ in mechanism order.
2. Computes $P=\rho R_uT\sum_iY_i/W_i$ from that state.
3. Evaluates thermodynamic properties and rate constants.
4. Replaces PLOG coefficients using the current $T$ and $P$.
5. Forms concentration products, reverse rates, and net progress rates.
6. Maps progress rates to species production rates.
7. Returns the mass-fraction and energy-conservation derivatives.

The Jacobian callback recomputes pressure from its own input state.
It therefore does not depend on an RHS callback having run immediately
before it at the same state.

| Mathematical task | Implementation |
| --- | --- |
| State contract, solver dispatch, tolerances | `chemistry_ODE_integrate`, `ODE_solver_speedchem_init` in [chemistry_input.f90](../src/chemistry_input.f90) |
| Setup defaults and environment override | `chemistry_setup::solver_setup` in [SCmodule.f90](../src/SCmodule.f90) |
| Constant-volume species and energy RHS | `SC_conV` in [SCconV.f90](../src/SCconV.f90) |
| Density-setting integration wrapper | `conV_integrate` in [SCconV.f90](../src/SCconV.f90) |
| Species and mixture thermodynamics, pressure | `SCspeciesthermo`, `SCmixturethermo` in [SCmodule.f90](../src/SCmodule.f90) |
| Equilibrium constants and PLOG | `reacpar` in [SCmodule.f90](../src/SCmodule.f90) |
| Rates of progress and production rates | `kinetics_mod::SCdwdt` in [SCmodule.f90](../src/SCmodule.f90) |
| Analytical Jacobian and backend adapters | `constV_jac_sparse`, `constV_jac_*` in [SCconV.f90](../src/SCconV.f90) |
| Sparse matrix representation and operations | [SCsparse_definitions.f90](../src/SCsparse_definitions.f90), [SCsparse.f90](../src/SCsparse.f90) |

This document explains the current equations and implementation. It does
not supply a reaction mechanism or an experimental validation claim.
The compiler-interface limitation associated with LSODES workspace
aliasing is tracked separately in [issue #12](https://github.com/GTT10/sppedchem/issues/12).

## References

- F. Perini, E. Galligani, and R. D. Reitz (2012),
  [An Analytical Jacobian Approach to Sparse Reaction Kinetics for Computationally Efficient Combustion Modeling with Large Reaction Mechanisms](https://doi.org/10.1021/ef300747n),
  *Energy & Fuels* 26, 4804–4822. Original SpeedCHEM analytical-Jacobian method.
- F. Perini, E. Galligani, and R. D. Reitz (2014),
  [A Study of Direct and Krylov Iterative Sparse Solver Techniques to Approach Linear Scaling of the Integration of Chemical Kinetics with Detailed Combustion Mechanisms](https://doi.org/10.1016/j.combustflame.2013.11.017),
  *Combustion and Flame* 161, 1180–1195. Sparse linear-solver techniques.
- [ODEPACK project documentation, LLNL](https://computing.llnl.gov/projects/odepack)
  and [original source distribution, Netlib](https://netlib.org/odepack/).
- [E. Hairer's ODE software descriptions](https://www.unige.ch/~hairer/software.html).
- [Cantera: NASA-7 species thermodynamic relations](https://cantera.org/stable/reference/thermo/species-thermo.html#the-nasa-7-coefficient-polynomial-parameterization).

The PLOG extension described in the companion document belongs to this
community-maintained version. The original SpeedCHEM references above
are not evidence that the extension appeared in those papers.
