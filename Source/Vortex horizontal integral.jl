begin
	function vortexHorizontalIntegral(x, y, sourcePanel)
		lengthSource = L[sourcePanel];
		xStart = xPanel[sourcePanel];
		yStart = yPanel[sourcePanel];
		tangentSource = ω[sourcePanel];
		A = (xStart - x)*cos(tangentSource) + (yStart - y)*sin(tangentSource);
		B = (x - xStart)^2 + (y - yStart)^2;
		C = sin(tangentSource);
		D = (yStart - y);
		E = sqrt(B - A^2);
		dummyOne = atan(lengthSource + A, E);
		dummyTwo = atan(A, E);
		if (isnan(C/2*log((lengthSource^2 + 2*A*lengthSource + B)/B) + (D - A*C)/E*(dummyOne - dummyTwo)) == false) && (isinf(C/2*log((lengthSource^2 + 2*A*lengthSource + B)/B) + (D - A*C)/E*(dummyOne - dummyTwo)) == false) && ((lengthSource^2 + 2*A*lengthSource + B)/B > 0)
			QX = C/2*log((lengthSource^2 + 2*A*lengthSource + B)/B) + (D - A*C)/E*(dummyOne - dummyTwo);
		else
			QX = 0;
		end
		return QX
	end
end
