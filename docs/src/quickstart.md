```@meta
CurrentModule = MathOptComplements
```
## Quickstart

The following code shows how to solve a simple Mathematical Program with Complementarity Constraints (MPCC) with Ipopt. This instance is a JuMP translation
of `scholtes4.mod` in [MacMPEC](https://www.mcs.anl.gov/~leyffer/macmpec/comments.html).

### Basic usage

We start by writing the model with JuMP:

```@example quickstart
using JuMP
z0 = [0, 1]
model = Model()
@variable(model, z[i=1:2] >= 0.0, start=z0[i])
@variable(model, z3, start=0.0)
@objective(model, Min, z[1] + z[2] - z3)
@constraint(model, -4 * z[1] + z3 <= 0)
@constraint(model, -4 * z[2] + z3 <= 0)
@constraint(model, [z[1], z[2]] ∈ MOI.Complements(2))
model
```

Solving this instance with Ipopt simply amounts to:
```@example quickstart
using MathOptComplements
using Ipopt
MathOptComplements.Bridges.add_all_bridges(model)
set_optimizer(model, Ipopt.Optimizer)
JuMP.optimize!(model)
println("Solution: ", JuMP.value.(model[:z]))
```
Under the hood, MathOptComplements takes the complementarity
constraints and reformulate it as a nonlinear constraint using
the [`ScholtesRelaxation`](@ref) method.

!!! note
    We recommend setting the following options in Ipopt for optimal performance:
    ```julia
    JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
    JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
    JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)
    ```


### Changing the relaxation

The [`ScholtesRelaxation`](@ref) is used by default, but the user
has the freedom to use any of the relaxations implemented in MathOptComplements.
Replacing the [`ScholtesRelaxation`](@ref) by the classical [`FischerBurmeisterRelaxation`](@ref) simply amounts to
```@example quickstart
JuMP.set_optimizer(model, () -> MathOptComplements.Optimizer(Ipopt.Optimizer()))
MOI.set(model, MathOptComplements.DefaultComplementarityReformulation(), MathOptComplements.FischerBurmeisterRelaxation(1e-8))
JuMP.optimize!(model)

println("Solution: ", JuMP.value.(model[:z]))
```

Observe that the solution is here closer to the true solution `(0, 0)`
than the solution returned by the Scholtes relaxation.

!!! warning
    MPCCs are nonconvex problems and they rarely have a unique solution.
    In general, changing the relaxation method can yield a different local solution.

