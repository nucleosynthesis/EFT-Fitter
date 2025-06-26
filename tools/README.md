## RBF Splines

The spline makes use of a radial basis decomposition to produce a continous $N \to 1$ map (function) from $M$ provided sample points. The function of the $N$ variables $\vec{x}$
is assumed to be of the form,

$$
f(\vec{x}) = \sum_{i=1}^{M}w_{i}\phi(||\vec{x}-\vec{x}_{i}||),
$$

where for example $\phi(||\vec{z}||) = e^{-\dfrac{||\vec{z}||}{\epsilon^{2}}}$. 

The distance $||.||$ between two points is given by,

$$
||\vec{x}-\vec{y}||  = \sum_{j=1}^{N}(x_{j}-y_{j})^{2},
$$

Given the sample points, it is possible to determine the weights $w_{i}$ as the solution of the set of equations,

$$
\sum_{i=1}^{M}w_{i}\phi(||\vec{x}_{j}-\vec{x}_{i}||) = f(\vec{x}_{j}).
$$


The solution is obtained by inverting this matrix equation. 