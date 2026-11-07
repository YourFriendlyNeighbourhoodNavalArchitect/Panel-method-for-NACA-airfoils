# Read Selig-format coordinates, skipping header lines.
function readAirfoil(path::AbstractString)
    x = Float64[]
    y = Float64[]
    for line in eachline(path)
        fields = split(line)
        length(fields) ≥ 2 || continue
        xi = tryparse(Float64, fields[1])
        yi = tryparse(Float64, fields[2])
        (xi === nothing || yi === nothing) && continue
        push!(x, xi)
        push!(y, yi)
    end
    length(x) ≥ 3 || throw(ArgumentError("`$path` does not contain airfoil coordinates"))
    return x, y
end

# NACA 4-digit airfoil of unit chord in Selig order, with `n` cosine-spaced stations per surface.
function naca4(code::AbstractString; n::Integer = 101, closedTE::Bool = true)
    (length(code) == 4 && all(isdigit, code)) ||
        throw(ArgumentError("NACA 4-digit code must be four digits, got \"$code\""))
    n ≥ 3 || throw(ArgumentError("need at least 3 stations per surface, got $n"))
    m = parse(Int, code[1]) / 100
    p = parse(Int, code[2]) / 10
    t = parse(Int, code[3:4]) / 100
    a4 = closedTE ? -0.1036 : -0.1015

    xc = [(1 - cospi(k / (n - 1))) / 2 for k in 0:n-1]
    yt = @. 5t * (0.2969 * sqrt(xc) - 0.1260xc - 0.3516xc^2 + 0.2843xc^3 + a4 * xc^4)
    yc, slope = naca4Camber(m, p, xc)
    θ = atan.(slope)

    xu = @. xc - yt * sin(θ);  yu = @. yc + yt * cos(θ)
    xl = @. xc + yt * sin(θ);  yl = @. yc - yt * cos(θ)
    return [reverse(xu); xl[2:end]], [reverse(yu); yl[2:end]]
end

function naca4Camber(m, p, x)
    yc = zero(x)
    slope = zero(x)
    (m == 0 || p == 0) && return yc, slope
    for (k, xk) in pairs(x)
        if xk < p
            yc[k] = m / p^2 * (2p * xk - xk^2)
            slope[k] = 2m / p^2 * (p - xk)
        else
            yc[k] = m / (1 - p)^2 * (1 - 2p + 2p * xk - xk^2)
            slope[k] = 2m / (1 - p)^2 * (p - xk)
        end
    end
    return yc, slope
end

# Read-only view of a vector, so that `Panels` data cannot change after construction.
struct ReadOnlyVector{T} <: AbstractVector{T}
    data::Vector{T}
end

Base.size(v::ReadOnlyVector) = size(v.data)
Base.IndexStyle(::Type{<:ReadOnlyVector}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::ReadOnlyVector, i::Int) = v.data[i]
Base.setindex!(::ReadOnlyVector, _, ::Int) =
    throw(ArgumentError("Panels data are read-only; construct a new Panels instead"))

# Linear panels between consecutive counter-clockwise nodes. Control points at midpoints, with normal vectors outward.
struct Panels{T<:AbstractFloat}
    xe::ReadOnlyVector{T}
    ye::ReadOnlyVector{T}
    xc::ReadOnlyVector{T}
    yc::ReadOnlyVector{T}
    len::ReadOnlyVector{T}
    tx::ReadOnlyVector{T}
    ty::ReadOnlyVector{T}
    nx::ReadOnlyVector{T}
    ny::ReadOnlyVector{T}

    function Panels(xe::AbstractVector{<:Real}, ye::AbstractVector{<:Real})
        length(xe) == length(ye) || throw(DimensionMismatch("x and y must have equal length"))
        length(xe) ≥ 4 || throw(ArgumentError("need at least 3 panels"))
        T = float(promote_type(eltype(xe), eltype(ye)))
        xe, ye = collect(T, xe), collect(T, ye)
        signedArea(xe, ye) > 0 ||
            throw(ArgumentError("nodes must be ordered counter-clockwise (TE → upper → LE → lower → TE)"))
        dx, dy = diff(xe), diff(ye)
        len = hypot.(dx, dy)
        all(>(0), len) || throw(ArgumentError("coincident consecutive nodes produce zero-length panels"))
        tx, ty = dx ./ len, dy ./ len
        xc = (xe[1:end-1] .+ xe[2:end]) ./ 2
        yc = (ye[1:end-1] .+ ye[2:end]) ./ 2
        return new{T}(map(ReadOnlyVector, (xe, ye, xc, yc, len, tx, ty, ty, -tx))...)
    end
end

nPanels(p::Panels) = length(p.len)
Base.length(p::Panels) = nPanels(p)
Base.show(io::IO, p::Panels) = print(io, "Panels(", nPanels(p), " panels, chord = ", chordLength(p), ")")

function chordLength(p::Panels)
    xMin, xMax = extrema(p.xe)
    return xMax - xMin
end

leadingEdge(p::Panels) = argmin(p.xe)
upperSurface(p::Panels) = 1:leadingEdge(p)-1
lowerSurface(p::Panels) = leadingEdge(p):nPanels(p)

# Shoelace formula.
function signedArea(x, y)
    a = zero(eltype(x))
    @inbounds for k in eachindex(x)
        k1 = k == lastindex(x) ? firstindex(x) : k + 1
        a += x[k] * y[k1] - x[k1] * y[k]
    end
    return a / 2
end

# Re-panel `(x, y)` with `n` cosine-spaced panels, linearly interpolating each surface.
function discretizeAirfoil(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, n::Integer)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have equal length"))
    (n ≥ 4 && iseven(n)) || throw(ArgumentError("number of panels must be even and ≥ 4, got $n"))
    x, y = float.(x), float.(y)
    if signedArea(x, y) < 0
        x, y = reverse(x), reverse(y)
    end
    iLE = argmin(x)
    upper = surface(reverse(x[1:iLE]), reverse(y[1:iLE]), "upper")
    lower = surface(x[iLE:end], y[iLE:end], "lower")

    half = n ÷ 2
    s = [(1 - cospi(k / half)) / 2 for k in 0:half]
    xu = upper.x[1] .+ (upper.x[end] - upper.x[1]) .* s
    xl = lower.x[1] .+ (lower.x[end] - lower.x[1]) .* s
    yu = interp.(Ref(upper), xu)
    yl = interp.(Ref(lower), xl)
    return Panels([reverse(xu); xl[2:end]], [reverse(yu); yl[2:end]])
end

function surface(x, y, name)
    all(>(0), diff(x)) ||
        throw(ArgumentError("$name surface abscissae must increase monotonically from the leading edge"))
    return (; x, y)
end

function interp(s, xq)
    k = clamp(searchsortedlast(s.x, xq), firstindex(s.x), lastindex(s.x) - 1)
    w = (xq - s.x[k]) / (s.x[k+1] - s.x[k])
    return (1 - w) * s.y[k] + w * s.y[k+1]
end
