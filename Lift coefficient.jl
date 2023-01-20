begin
	function liftCoefficient()
		#Sum of element values of length panel array.
		dummyOne = sum(L);
		#Vortex strength, as calculated above. It is the last element of the solution array.
		dummyTwo = solutionMatrix[length(solutionMatrix)];
		#Chord length.
		dummyThree = maximum(xPanel) - minimum(xPanel);
		return -2*dummyOne*dummyTwo/(dummyThree*uInfinity)
	end
	liftCoefficient()
end
