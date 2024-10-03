begin
	function sourceNormalIntegral(controlPanel, sourcePanel)
		#controlPoint refers to the i-pointer, whereas sourcePanel refers to the j-pointer. Function evaluates the effect of the j-th panel on the i-th Panel.
		lengthSource = L[sourcePanel];
		xStart = xPanel[sourcePanel];
		yStart = yPanel[sourcePanel];
		tangentControl = ω[controlPanel];
		tangentSource = ω[sourcePanel];
		xControlPoint = xM[controlPanel];
		yControlPoint = yM[controlPanel];
		A = (xStart - xControlPoint)*cos(tangentSource) + (yStart - yControlPoint)*sin(tangentSource);
		B = (xControlPoint - xStart)^2 + (yControlPoint - yStart)^2;
		C = sin(tangentSource - tangentControl);
		D = (xControlPoint - xStart)*sin(tangentControl) + (yStart - yControlPoint)*cos(tangentControl);
		if B - A^2 >= 0
			E = sqrt(B - A^2);
		else
			E = 0;
		end
		dummyOne = atan(lengthSource + A, E);
		dummyTwo = atan(A, E);
		if (isnan(C/2*log((lengthSource^2 + 2*A*lengthSource + B)/B) + (D - A*C)/E*(dummyOne - dummyTwo)) == false) && (isinf(C/2*log((lengthSource^2 + 2*A*lengthSource + B)/B) + (D - A*C)/E*(dummyOne - dummyTwo)) == false) && ((lengthSource^2 + 2*A*lengthSource + B)/B > 0) 
			I = C/2*log((lengthSource^2 + 2*A*lengthSource + B)/B) + (D - A*C)/E*(dummyOne - dummyTwo);
		else
			I = 0;
		end
		return I
	end
end
