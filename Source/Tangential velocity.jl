begin
	function tangentialVelocity(controlPanel)
		#Seperation of terms for readability.
		freestreamTerm = uInfinity*sin(convertAngle - φ[controlPanel]);
		sourceTerm = 0;
		vortexTerm = solutionMatrix[length(xPanel)]/2;
		for i = 1:(length(xPanel) - 1)
			if i == controlPanel
				continue
			else
				sourceTerm = sourceTerm + solutionMatrix[i]/(2pi)*sourceTangentialIntegral(controlPanel, i);
				vortexTerm = vortexTerm + solutionMatrix[length(xPanel)]/(2pi)*vortexTangentialIntegral(controlPanel, i);
			end
		end
		return freestreamTerm + sourceTerm + vortexTerm
	end
end
