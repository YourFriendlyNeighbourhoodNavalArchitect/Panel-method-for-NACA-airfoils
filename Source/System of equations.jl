begin
	#First we populate the LHS matrix of known quantities.
	convertAngle = angleOfAttack*pi/180;
	A = zeros(length(xPanel), length(xPanel));
	for i = 1:(length(xPanel) - 1)
		#Populates the submatrix defined by the first N rows and the N first columns.
		for j = 1:(length(xPanel) - 1)
			if i == j
				A[i, j] = pi;
			else
				A[i, j] = sourceNormalIntegral(i, j);
			end
		end
		#Populates the first 100 elements of the 101th column.
		for j = 1:(length(xPanel) - 1)
			if i == j
				continue
			else 
				A[i, length(xPanel)] = A[i, length(xPanel)] + vortexNormalIntegral(i, j);
			end
		end
		#Populates the first 100 elements of the 101th row.
		if i == 1
			A[length(xPanel), i] = sourceTangentialIntegral((length(xPanel) - 1), 1);
		elseif i == (length(xPanel) - 1)
			A[length(xPanel), i] = sourceTangentialIntegral(1, (length(xPanel) - 1));
		else
			A[length(xPanel), i] = sourceTangentialIntegral(1, i) + sourceTangentialIntegral((length(xPanel) - 1), i);
		end
		#Populates the remaining element.
		if i == 1 || i == (length(xPanel) - 1)
			continue
		else
			A[length(xPanel), length(xPanel)] = A[length(xPanel), length(xPanel)] + vortexTangentialIntegral(1, i) + vortexTangentialIntegral((length(xPanel) - 1), i);
		end
	end
	A[length(xPanel), length(xPanel)] = A[length(xPanel), length(xPanel)] + vortexTangentialIntegral((length(xPanel) - 1), 1) + vortexTangentialIntegral(1, (length(xPanel) - 1)) + 2pi;
	#Then we populate the RHS matrix of known quantities.
	B = zeros(length(xPanel), 1);
	for i = 1:(length(xPanel) - 1)
		B[i, 1] = -2pi*uInfinity*cos(convertAngle - φ[i]);
	end
	B[length(xPanel), 1] = -2pi*uInfinity*(sin(convertAngle - φ[1]) + sin(convertAngle - φ[(length(xPanel) - 1)]));
	solutionMatrix = A\B;
end
