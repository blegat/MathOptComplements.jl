# MathOptComplements

| **Build Status** |
|:----------------:|
| [![Build Status][build-img]][build-url] [![Codecov branch][codecov-img]][codecov-url] |

This package provides a set of routines to reformulate as a nonlinear program any JuMP's problem with
[complementarity constraints](https://jump.dev/JuMP.jl/v1.29.3/moi/reference/standard_form/#MathOptInterface.Complements).
Once reformulated, the problem can be solved using a nonlinear programming solver like [Ipopt](https://coin-or.github.io/Ipopt/).

## Quick-start

We aim at solving the following mathematical program with complementarity constraints (MPCC):
```julia
using JuMP

model = Model()
@variable(model, z[1:2])
set_lower_bound(z[2], 0)
@objective(model, Min, (z[1] - 1)^2 + z[2]^2)
@constraint(model, [z[2] - z[1], z[2]] ∈ MOI.Complements(2))
```
MathOptComplements reformulates automatically all the complementarity
constraints using [MOI.Bridges](https://jump.dev/MathOptInterface.jl/stable/submodules/Bridges/overview/).
Solving the problem with Ipopt and MathOptComplements just amounts to either do:
```julia
using Ipopt
using MathOptComplements

MathOptComplements.Bridges.add_all_bridges(model)
set_optimizer(model, Ipopt.Optimizer)
```
or
```julia
using Ipopt
using MathOptComplements

set_optimizer(model, () -> MathOptComplements.Optimizer(Ipopt.Optimizer()))
```
before you call
```julia
JuMP.optimize!(model)
```

We recommend setting the following options in Ipopt to improve the performance:
```julia
JuMP.set_optimizer_attribute(model, "mu_strategy", "adaptive")
JuMP.set_optimizer_attribute(model, "bound_push", 1e-1)
JuMP.set_optimizer_attribute(model, "bound_relax_factor", 0.0)

```

## Supported reformulations

You can change the reformulation by using the optimizer attribute `MathOptComplements.DefaultComplementarityReformulation`:
```julia
MOI.set(model, MathOptComplements.DefaultComplementarityReformulation(), MathOptComplements.ScholtesRelaxation(0.0))
```

> [!note]
> The `MathOptComplements.DefaultComplementarityReformulation` attribute only works if you used
> `set_optimizer(model, () -> MathOptComplements.Optimizer(...))`, not `MathOptComplements.Bridges.add_all_bridges(model)`.
> That is the only difference between the two though so if you are not using setting
> `MathOptComplements.DefaultComplementarityReformulation` because you don't change the default reformulation
> or because you set it constraint-wise with the constraint attribute `MathOptComplements.ComplementarityReformulation`
> then `MathOptComplements.Bridges.add_all_bridges(model)` will work.

MathOptComplements supports the following reformulations:
- `MathOptComplements.ScholtesRelaxation(tau)` (**default**): reformulates the complementarity `0 ≤ a ⟂ b ≥ 0` as `0 ≤ (a, b)` and `a b ≤ tau`. For `tau = 0`, the reformulation is exact and leads to the formulation of a degenerate nonlinear program. The higher the parameter `tau`, the better the behavior in Ipopt.
- `MathOptComplements.FischerBurmeisterRelaxation(tau)`: reformulates the complementarity `0 ≤ a ⟂ b ≥ 0` as `0 ≤ (a, b)` and `a + b - sqrt((a+b)^2 + tau) ≤ 0`.
- `MathOptComplements.LiuFukushimaRelaxation(tau)`: reformulates the complementarity `0 ≤ a ⟂ b ≥ 0` as `a b ≤ tau^2` and `(a + tau)(b + tau) ≥ tau^2`.
- `MathOptComplements.KanzowSchwarzRelaxation(tau)`: reformulates the complementarity `0 ≤ a ⟂ b ≥ 0` as `0 ≤ (a, b)` and `ϕ(a, b) ≤ 0`, with `ϕ(a, b) = (a - tau)(b - tau)` if `a + b > 2tau`, `-0.5((a -tau)^2 + (b - tau)^2)` otherwise.

Most reformulations are not equivalent to the original problem, explaining why they are not activated by default.
You can find [here](https://arxiv.org/html/2312.11022v2) a recent benchmark comparing the different reformulations on [MacMPEC](https://www.mcs.anl.gov/~leyffer/macmpec/).

## Funding
We acknowledge support from the [Fondation Mathématiques Jacques Hadamard](https://www.fondation-hadamard.fr/fr/)
which has funded the PGMO-IROE project "A new optimization suite for large-scale market equilibrium".

[build-img]: https://github.com/jump-dev/MathOptComplements.jl/actions/workflows/ci.yml/badge.svg?branch=main
[build-url]: https://github.com/jump-dev/MathOptComplements.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/jump-dev/MathOptComplements.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jump-dev/MathOptComplements.jl/branch/main
