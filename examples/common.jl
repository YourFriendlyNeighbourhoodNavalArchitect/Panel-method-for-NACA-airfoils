# Shared output for the examples: printed coefficients and the panel, pressure and velocity-field plots.

using PanelMethod
using Plots

# Print the force and moment coefficients, with the exact lift coefficient when one is known.
function report(sol::PanelSolution; clExact = nothing)
    println(sol)
    println("Cl (Kutta–Joukowski)   = ", round(liftCoefficient(sol); digits = 4))
    println("Cl (pressure integral) = ", round(pressureLiftCoefficient(sol); digits = 4))
    isnothing(clExact) || println("Cl (exact)             = ", round(clExact; digits = 4))
    println("Cm (quarter chord)     = ", round(momentCoefficient(sol); digits = 4))
    println("Cp (minimum)           = ", round(minimum(sol.cp); digits = 4))
end

# Write the panel, pressure and velocity-field plots to `outDir`. `cpExact = (x, cp)` adds the exact surface Cp.
function plotCase(outDir, title, x, y, sol::PanelSolution; cpExact = nothing)
    mkpath(outDir)
    panels = sol.panels
    c = chordLength(panels)
    xMin, xMax = extrema(panels.xe)
    yMin, yMax = extrema(panels.ye)
    xOverC(xs) = (xs .- xMin) ./ c
    common = (; aspect_ratio = :equal, size = (1080, 540), framestyle = :box, margin = 5Plots.mm, title)

    # Panels, control points and outward normals on every other panel.
    k = 1:2:nPanels(panels)
    fig = plot(x, y; label = "Input geometry", color = :gray, linewidth = 2,
               xlims = (xMin - 0.1c, xMax + 0.1c), ylims = (yMin - 0.15c, yMax + 0.15c), common...)
    plot!(fig, panels.xe, panels.ye; label = "Panels", color = :black, marker = (:circle, 2))
    quiver!(fig, panels.xc[k], panels.yc[k]; quiver = (0.06c .* panels.nx[k], 0.06c .* panels.ny[k]),
            arrow = arrow(:simple, :head), color = :red)
    savefig(fig, joinpath(outDir, "panels.png"))

    # Pressure coefficient.
    upper, lower = upperSurface(panels), lowerSurface(panels)
    fig = plot(xOverC(panels.xc[upper]), sol.cp[upper]; label = "Upper surface", marker = (:circle, 3),
               yflip = true, xlabel = "x / c", ylabel = "Cp", size = (1080, 540), margin = 5Plots.mm, title)
    plot!(fig, xOverC(panels.xc[lower]), sol.cp[lower]; label = "Lower surface", marker = (:circle, 3))
    isnothing(cpExact) ||
        plot!(fig, xOverC(cpExact[1]), cpExact[2]; label = "Exact", color = :black, linestyle = :dash)
    savefig(fig, joinpath(outDir, "pressure.png"))

    # Speed contours and velocity vectors.
    xs = range(xMin - 0.5c, xMax + 0.5c; length = 241)
    ys = range(-0.4c, 0.4c; length = 121)
    U, V = velocityField(sol, xs, ys)
    fig = contourf(xs, ys, hypot.(U, V) ./ sol.vInf; levels = 30, color = :viridis,
                   colorbar_title = "|V| / V∞", common...)
    plot!(fig, Shape(panels.xe, panels.ye); color = :white, linecolor = :black, label = false)
    stride = 12
    X = [xi for _ in ys, xi in xs][1:stride:end, 1:stride:end]
    Y = [yi for yi in ys, _ in xs][1:stride:end, 1:stride:end]
    Us, Vs = U[1:stride:end, 1:stride:end], V[1:stride:end, 1:stride:end]
    keep = .!isnan.(Us)
    scale = 0.06c / sol.vInf
    quiver!(fig, X[keep], Y[keep]; quiver = (scale .* Us[keep], scale .* Vs[keep]), color = :white)
    savefig(fig, joinpath(outDir, "velocityField.png"))

    println("Figures written to ", outDir)
end
