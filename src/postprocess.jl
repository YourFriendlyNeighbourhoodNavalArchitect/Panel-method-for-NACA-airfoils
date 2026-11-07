# Lift coefficient from the Kutta–Joukowski theorem.
liftCoefficient(s::PanelSolution) = 2 * s.γ * sum(s.panels.len) / (s.vInf * chordLength(s.panels))

# Lift coefficient from surface pressure integration.
function pressureLiftCoefficient(s::PanelSolution)
    p = s.panels
    cx = -sum(@. s.cp * p.len * p.nx)
    cy = -sum(@. s.cp * p.len * p.ny)
    sα, cα = sincos(s.alpha)
    return (cy * cα - cx * sα) / chordLength(p)
end

# Pitching-moment coefficient about `(xRef, yRef)`, default quarter chord, positive nose-up.
function momentCoefficient(s::PanelSolution; xRef = nothing, yRef = 0)
    p = s.panels
    c = chordLength(p)
    xr = something(xRef, minimum(p.xe) + c / 4)
    return sum(@. s.cp * p.len * ((p.xc - xr) * p.ny - (p.yc - yRef) * p.nx)) / c^2
end

# Flow velocity `(u, v)` at a point off the airfoil surface.
function flowVelocity(s::PanelSolution{T}, x::Real, y::Real) where {T}
    p = s.panels
    su, sv = zero(T), zero(T)
    gu, gv = zero(T), zero(T)
    @inbounds for j in 1:nPanels(p)
        mu, mv = sourceVelocity(p, j, x, y)
        su += s.λ[j] * mu
        sv += s.λ[j] * mv
        gu += mu
        gv += mv
    end
    sα, cα = sincos(s.alpha)
    twoπ = 2 * T(π)
    u = s.vInf * cα + (su + s.γ * gv) / twoπ
    v = s.vInf * sα + (sv - s.γ * gu) / twoπ
    return u, v
end

# Velocity components on the grid `xs × ys` as `(length(ys), length(xs))` matrices, `NaN` inside the airfoil and on its surface.
function velocityField(s::PanelSolution{T}, xs::AbstractVector, ys::AbstractVector) where {T}
    U = fill(T(NaN), length(ys), length(xs))
    V = copy(U)
    for (j, x) in enumerate(xs), (i, y) in enumerate(ys)
        (isInside(s.panels, x, y) || onSurface(s.panels, x, y)) && continue
        U[i, j], V[i, j] = flowVelocity(s, x, y)
    end
    return U, V
end

# Whether `(x, y)` lies inside the panel polygon.
function isInside(p::Panels, x::Real, y::Real)
    inside = false
    xe, ye = p.xe, p.ye
    k1 = lastindex(xe)
    @inbounds for k in eachindex(xe)
        if (ye[k] > y) != (ye[k1] > y) &&
           x < xe[k] + (y - ye[k]) * (xe[k1] - xe[k]) / (ye[k1] - ye[k])
            inside = !inside
        end
        k1 = k
    end
    return inside
end

# Whether `(x, y)` lies on a panel within `tol`, where the panel velocities are singular.
function onSurface(p::Panels{T}, x::Real, y::Real; tol = sqrt(eps(T)) * chordLength(p)) where {T}
    @inbounds for j in 1:nPanels(p)
        dx, dy = x - p.xe[j], y - p.ye[j]
        ξ = dx * p.tx[j] + dy * p.ty[j]
        η = dy * p.tx[j] - dx * p.ty[j]
        abs(η) ≤ tol && -tol ≤ ξ ≤ p.len[j] + tol && return true
    end
    return false
end
