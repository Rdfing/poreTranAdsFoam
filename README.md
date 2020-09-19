# transportAdsorptionFoam

![alt text](https://github.com/Rdfing/transportAdsorptionFoam/blob/master/sample_1k_sphere_DOF_10M.png?raw=true)

A solver for pore-scale species transport and surface adsorption simulations.

For the bulk transport of the chemical species, 

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial Y_{j}}{\partial t}%2B\nabla\cdot\left(\mathbf{u}Y_{j}\right)=\nabla\cdot\left(D_{j}\nabla Y_{j}\right).">

The surface reaction of adsorption (<img src="https://render.githubusercontent.com/render/math?math=N_{ads}">) and desorption (<img src="https://render.githubusercontent.com/render/math?math=N_{des}">) is considered through a partially lumped dynamic Langmuir-like model in the boundary condition.

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial Y_{j}^{ads}}{\partial t}=N_{ads}-N_{des},">

<img src="https://render.githubusercontent.com/render/math?math=N_{ads}=k_{ads}Y_{j}\left(1-\theta\right)\quad and\quad\theta=\frac{Y_{j}^{ads}}{\Gamma_{j}^{ads}},">

<img src="https://render.githubusercontent.com/render/math?math=N_{des}=k_{des}Y_{j}^{ads}.">

Sucess of the application of the model requires inverse modeling of the kinetic parameters of the adsorption.

The solver is tested on OpenFOAM-v1912

Further detailes are provided in the incoming manuscript.

Lumped version of the solver is coming too
 
