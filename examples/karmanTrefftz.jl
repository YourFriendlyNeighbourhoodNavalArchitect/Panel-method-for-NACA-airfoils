# Pressure distribution and velocity field for a thick, highly cambered Kármán–Trefftz hydrofoil.
# Plots lives in the separate `examples` environment. From the repository root, set it up once with
#   julia --project=examples -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
# and run
#   julia --project=examples examples/karmanTrefftz.jl

include("common.jl")

const vInf = 1.0
const alphaDeg = 2.0
const panelCount = 100

# Geometry: the circle `ζ = ζ0 + a·exp(iθ)` through `ζ = 1`, mapped by `z = n((ζ+1)^n + (ζ-1)^n) / ((ζ+1)^n - (ζ-1)^n)`.
# `μx` sets the thickness, `μy` the camber and `n` the trailing-edge angle `(2 - n)·180°`.
const μx, μy, n = 0.12, 0.12, 1.92
const ζ0 = complex(-μx, μy)
const a = abs(1 - ζ0)
const β = atan(μy, 1 + μx)
ratio(ζ) = (ζ - 1) / (ζ + 1)
mapping(ζ) = n * (1 + ratio(ζ)^n) / (1 - ratio(ζ)^n)
mappingDerivative(ζ) = 4n^2 * ratio(ζ)^(n - 1) / ((ζ + 1)^2 * (1 - ratio(ζ)^n)^2)
θ = range(-β, 2π - β; length = 801)
z = mapping.(ζ0 .+ a .* cis.(θ))
x, y = real.(z), imag.(z)
x[1] = x[end] = n
y[1] = y[end] = 0

# Reference solution: the circulation that puts the rear stagnation point on the trailing edge, and the
# surface velocity `dW/dz = (dW/dζ) / (dz/dζ)`, which is singular at the trailing edge itself.
α = deg2rad(alphaDeg)
Γ = 4π * a * vInf * sin(α + β)
clExact = 2Γ / (vInf * (maximum(x) - minimum(x)))
ζs = ζ0 .+ a .* cis.(θ[2:end-1])
dWdζ = @. vInf * (cis(-α) - a^2 * cis(α) / (ζs - ζ0)^2) + im * Γ / (2π * (ζs - ζ0))
cpExact = (x[2:end-1], @. 1 - abs2(dWdζ / mappingDerivative(ζs)) / vInf^2)

panels = discretizeAirfoil(x, y, panelCount)
sol = solvePanels(panels; alphaDeg, vInf)

report(sol; clExact)
plotCase(joinpath(@__DIR__, "output", "karmanTrefftz"), "Kármán–Trefftz hydrofoil, α = $(alphaDeg)°", x, y, sol; cpExact)
