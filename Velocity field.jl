begin
	#Defines the grid and ignores points that lie inside the hydrofoil.
	streamlineGrid = meshgrid(LinRange(-0.5, 1.5, 25), LinRange(-0.25, 0.25, 20));
	hydrofoilPolygon = Point.(xHydrofoil,yHydrofoil);
	gridPoints = vec(Point.(streamlineGrid[:, 1], streamlineGrid[:, 2]));
	for i = 1:length(gridPoints)
		if isinside(gridPoints[i], hydrofoilPolygon; allowonedge = false) == true
			streamlineGrid[i, :] .= 0;
		end
	end
	streamlineGrid = streamlineGrid[(streamlineGrid[:, 1] .!= 0) .|| (streamlineGrid[:, 2] .!= 0), :];
	uGrid = uInfinity*cos(convertAngle)*ones(size(streamlineGrid, 1));
	vGrid = uInfinity*sin(convertAngle)*ones(size(streamlineGrid, 1));
	#Calculates the cartesian components of flow velocity at each grid point, by addition of the freestream, source and vortex contributions.
	for i = 1:(length(xPanel) - 1)
		for j = 1:size(streamlineGrid, 1)
			uGrid[j] = uGrid[j] + sourceHorizontalIntegral(streamlineGrid[j, 1], streamlineGrid[j, 2], i)*solutionMatrix[i]/(2pi) + vortexHorizontalIntegral(streamlineGrid[j, 1], streamlineGrid[j, 2], i)*solutionMatrix[length(solutionMatrix)]/(2pi);
			vGrid[j] = vGrid[j] + sourceVerticalIntegral(streamlineGrid[j, 1], streamlineGrid[j, 2], i)*solutionMatrix[i]/(2pi) + vortexVerticalIntegral(streamlineGrid[j, 1], streamlineGrid[j, 2], i)*solutionMatrix[length(solutionMatrix)]/(2pi);
		end
	end
	#Length of quiver arrow is proportional to the speed at each point.
	quiverScale = sqrt.(uGrid.^2 .+ vGrid.^2)./uInfinity;
	P4 = quiver(streamlineGrid[:, 1], streamlineGrid[:, 2], 
		quiver = 0.05.*(quiverScale.*uGrid, quiverScale.*vGrid),
		color = :red)
	P4 = plot!(xHydrofoil, yHydrofoil,
		aspect_ratio = :equal,
		color = :black,
		linewidth = 2,
		label = :false,
		xlim = (-0.5, 1.5),
		ylim = (-0.25, 0.25),
		size = (1080, 720))
end
