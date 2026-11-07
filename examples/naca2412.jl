# Pressure distribution and velocity field for NACA 2412.
# Plots lives in the separate `examples` environment. From the repository root, set it up once with
#   julia --project=examples -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
# and run
#   julia --project=examples examples/naca2412.jl

include("common.jl")

const vInf = 1.0
const alphaDeg = 5.0
const panelCount = 100

# Geometry: NACA 2412 coordinates in Selig format.
x, y = readAirfoil(joinpath(@__DIR__, "naca2412.dat"))

# Reference solution: none in closed form.
clExact, cpExact = nothing, nothing

panels = discretizeAirfoil(x, y, panelCount)
sol = solvePanels(panels; alphaDeg, vInf)

report(sol; clExact)
plotCase(joinpath(@__DIR__, "output", "naca2412"), "NACA 2412, α = $(alphaDeg)°", x, y, sol; cpExact)
