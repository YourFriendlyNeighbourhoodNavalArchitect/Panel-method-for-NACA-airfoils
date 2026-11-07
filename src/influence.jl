# Velocity induced at `(x, y)` by a unit source sheet on panel `j`.
# A clockwise unit vortex sheet induces the same velocity rotated by -90°.

@inline function sourceVelocity(p::Panels, j::Integer, x, y)
    @inbounds begin
        dx, dy = x - p.xe[j], y - p.ye[j]
        tx, ty, ℓ = p.tx[j], p.ty[j], p.len[j]
    end
    ξ = dx * tx + dy * ty
    η = dy * tx - dx * ty
    uξ = log(hypot(ξ, η) / hypot(ξ - ℓ, η))
    uη = atan(η, ξ - ℓ) - atan(η, ξ)
    return uξ * tx - uη * ty, uξ * ty + uη * tx
end

# Normal (`I`) and tangential (`J`) source influence matrices at the control points.
# Vortex coefficients follow as `J` (normal) and `-I` (tangential).

function influenceMatrices(p::Panels{T}) where {T}
    N = nPanels(p)
    I = Matrix{T}(undef, N, N)
    J = Matrix{T}(undef, N, N)
    @inbounds for j in 1:N, i in 1:N
        if i == j
            I[i, j] = π
            J[i, j] = 0
        else
            u, v = sourceVelocity(p, j, p.xc[i], p.yc[i])
            I[i, j] = u * p.nx[i] + v * p.ny[i]
            J[i, j] = u * p.tx[i] + v * p.ty[i]
        end
    end
    return I, J
end
