# Source strengths `λ`, vortex strength `γ`, surface tangential velocity `vt` and `cp` at the control points.
struct PanelSolution{T<:AbstractFloat}
    panels::Panels{T}
    vInf::T
    alpha::T
    λ::Vector{T}
    γ::T
    vt::Vector{T}
    cp::Vector{T}
end

function Base.show(io::IO, s::PanelSolution)
    print(io, "PanelSolution(", nPanels(s.panels), " panels, α = ", round(rad2deg(s.alpha); digits = 6),
          "°, Cl = ", round(liftCoefficient(s); digits = 4), ")")
end

# Solve flow tangency at every control point plus the Kutta condition at the trailing edge.
function solvePanels(p::Panels{T}; alphaDeg::Real = 0, vInf::Real = 1) where {T}
    vInf > 0 || throw(ArgumentError("freestream speed must be positive, got $vInf"))
    N = nPanels(p)
    α, V = deg2rad(T(alphaDeg)), T(vInf)
    sα, cα = sincos(α)
    twoπ = 2 * T(π)

    I, J = influenceMatrices(p)
    rowSumI = vec(sum(I; dims = 2))
    rowSumJ = vec(sum(J; dims = 2))
    vnInf = @. V * (cα * p.nx + sα * p.ny)
    vtInf = @. V * (cα * p.tx + sα * p.ty)

    A = Matrix{T}(undef, N + 1, N + 1)
    A[1:N, 1:N] .= I
    A[1:N, N+1] .= rowSumJ
    A[N+1, 1:N] .= @views J[1, :] .+ J[N, :]
    A[N+1, N+1] = -(rowSumI[1] + rowSumI[N])
    b = [-twoπ .* vnInf; -twoπ * (vtInf[1] + vtInf[N])]

    x = A \ b
    λ, γ = x[1:N], x[N+1]
    vt = vtInf .+ (J * λ .- γ .* rowSumI) ./ twoπ
    cp = @. 1 - (vt / V)^2
    return PanelSolution(p, V, α, λ, γ, vt, cp)
end
