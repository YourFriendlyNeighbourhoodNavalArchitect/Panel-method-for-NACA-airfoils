# Panel Method for NACA Airfoils

This repository contains a two-dimensional panel method for airfoils in incompressible, inviscid flow, written in Julia. It models the airfoil with constant-strength source and vortex panels (Hess–Smith method) and computes the pressure distribution, the lift and moment coefficients, and the surrounding velocity field.

## Project Overview

The airfoil contour is read from a Selig-format `.dat` file or generated for any NACA 4-digit section, then re-panelled with cosine spacing that clusters panels at the leading and trailing edges. Each panel carries its own source strength, and all panels share a single vortex strength. Flow tangency at every panel midpoint, together with the Kutta condition at the trailing edge, forms a linear system of size N + 1. From its solution the code derives the surface pressure coefficient, the lift coefficient (from both the Kutta–Joukowski theorem and pressure integration), the quarter-chord moment coefficient, and the velocity at any point in the flow.

## Repository Structure

- `src/PanelMethod.jl`: Module definition and public API.
- `src/geometry.jl`: Coordinate input, NACA 4-digit generator and panel discretisation.
- `src/influence.jl`: Closed-form panel influence kernel and influence matrices.
- `src/solver.jl`: Linear system of flow tangency and Kutta conditions.
- `src/postprocess.jl`: Force and moment coefficients, and the velocity field.
- `examples/common.jl`: Shared report and plots (panels, pressure distribution and velocity field) for the examples.
- `examples/naca2412.jl`: Worked example for NACA 2412.
- `examples/naca2412.dat`: NACA 2412 coordinates used by the example.
- `examples/karmanTrefftz.jl`: Worked example for a thick, highly cambered Kármán–Trefftz hydrofoil, compared with its exact solution.
- `examples/output/`: Plots produced by the examples.

## References

1. J. D. Anderson, *Fundamentals of Aerodynamics*, McGraw-Hill.
2. A. M. Kuethe and C.-Y. Chow, *Foundations of Aerodynamics*, Wiley.
