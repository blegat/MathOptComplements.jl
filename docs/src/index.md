
# MathOptComplements

MathOptComplements is a [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl/) extension for complementarity constraints.

## Motivation

MOI implements a set for mixed-complementarity constraints through [`MathOptInterface.Complements`](@extref).
However, few solvers support mixed-complementarity out of the box, often requiring
users to apply manual reformulations.

MathOptComplements provides a systematic way to handle the complementarity
constraints in JuMP and MathOptInterface by introducing a new set [`MathOptComplements.ComplementsWithSetType`](@ref).
The package provides a rich collection of tools for manipulating complementarity constraints, including:
- Equivalent **reformulations** in forms better suited to the target solver;
- Automatic **relaxations** that reformulate complementarity constraints as nonlinear constraints.

Under the hood, MathOptComplements extends [the bridge system](https://jump.dev/MathOptInterface.jl/stable/submodules/Bridges/overview/) implemented in MOI to optimally reformulate
the complementarity constraints within the model.



## Funding
We acknowledge support from the [Fondation Mathématiques Jacques Hadamard](https://www.fondation-hadamard.fr/fr/)
which has funded the PGMO-IROE project "A new optimization suite for large-scale market equilibrium".

