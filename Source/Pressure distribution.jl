begin
	#Calculates the pressure coefficient for any given panel.
	function pressureCoefficient(controlPanel)
		return 1 - (tangentialVelocity(controlPanel)/uInfinity)^2
	end
	#Plots the pressure coefficient values at each point.
	xUpper = zeros(length(xPanel) - 1);
	xLower = zeros(length(xPanel) - 1);
	pressureValuesUpper = zeros(length(xPanel) - 1);
	pressureValuesLower = zeros(length(xPanel) - 1);
	for i = 1:(length(xPanel) - 1)
		if yM[i] >= 0
			xUpper[i] = xM[i];
			pressureValuesUpper[i] = pressureCoefficient(i);
		else
			xLower[i] = xM[i]
			pressureValuesLower[i] = pressureCoefficient(i);
		end
	end
	xUpper = filter(!iszero, xUpper);
	xLower = filter(!iszero, xLower);
	pressureValuesUpper = filter(!iszero, pressureValuesUpper);
	pressureValuesLower = filter(!iszero, pressureValuesLower);
	P3 = plot(xUpper, pressureValuesUpper,
		title = "Pressure coefficient of NACA 2412",
		label = "Upper surface of airfoil",
		marker = :circle,
		markersize = 4,
		markercolor = :blue,
		xlims = (-0.2, 1.2),
		xticks = (0:0.2:1),
		aspect_ratio = 0.2,
		ylims = (minimum(pressureValuesUpper) - 0.5, maximum(pressureValuesLower) + 0.5),
		size = (1080, 720))
	P3 = plot!(xLower, pressureValuesLower,
		label = "Lower surface of airfoil",
		marker = :circle,
		markersize = 4,
    markercolor = :red)
end
