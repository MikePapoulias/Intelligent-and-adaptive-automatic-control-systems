Design of reference model adaptive controllers for systems with unknown parameters using the Matlab environment.

# Description

In this project, given a system with constant but unknown parameters, the goal is to design reference model adaptive controllers to ensure that the state variable vector of the original system asymptotically tracks the state variable vector of a desired reference model. The desired reference model is required to exhibit a specific damping ratio and natural frequency.

The control objective is achieved through a direct reference model adaptive controller, initially applied to the linearized version of the system with output feedback, and subsequently to its canonical form with state feedback.

Next, the system is subjected to bounded external pulse-like disturbances, requiring an investigation of its robustness for varying disturbance amplitudes.

The analysis and theoretical findings are validated through appropriate plots generated in Matlab, where the entire problem is simulated.
