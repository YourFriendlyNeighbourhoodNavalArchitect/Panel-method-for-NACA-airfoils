# Panel-method-for-NACA-airfoils
The entire project was carried out in Julia, a high-performance dynamic programming language, suited for computational science.
The code discretises airfoil into linear panels and utilizes the source/vortex method to model the flow around the airfoil.
User inputs the shape of the airfoil as a .dat file and values for the freestream velocity and angle of attack.
The mathematical background of the code relies heavily on fluid mechanics and aerodynamics, particularly on potential flow theory. For further reading, refer to [1] and [2].
Code calculates the pressure coefficient at each point on the airfoil and plots the respective distribution.
Subsequently it constructs the velocity field of the flow around the airfoil.

[1] Fundamentals of aerodynamics, J. Anderson

[2] Foundations of aerodynamics, A. M. Kuethe
