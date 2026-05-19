```@meta
CurrentModule = MathOptComplements
```
## Solving an equilibrium problem with MathOptComplements

In this tutorial, we show how to solve an equilbrium problem
using the different methods implemented in MathOptComplements.

We take an example from the [JuMP documentation](https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/complementarity/#Electricity-consumption),
implementing a risk neutral competitive equilibrium between a producer and a consumer of electricity.

```@example tutorial_equilibrium
using JuMP

function equilibrium_model()
    I = 90_000                     # Annualized capital cost
    C = 60                         # Operation cost per MWh
    τ = 8_760                      # Hours per year
    θ = [0.2, 0.2, 0.2, 0.2, 0.2]  # Scenario probabilities
    A = [300, 350, 400, 450, 500]  # Utility function coefficients
    B = 1                          # Utility function coefficients
    model = Model()
    @variable(model, x >= 0, start = 1)           # Installed capacity
    @variable(model, Q[ω=1:5] >= 0, start = 1)  # Consumption
    @variable(model, Y[ω=1:5] >= 0, start = 1)  # Production
    @variable(model, P[ω=1:5], start = 1)       # Electricity price
    @variable(model, μ[ω=1:5] >= 0, start = 1)  # Capital scarcity margin
    # Unit investment cost equals annualized scarcity profit or investment is 0
    @constraint(model, I - τ * θ' * μ ⟂ x)
    # Difference between price and scarcity margin is equal to operation cost
    @constraint(model, [ω = 1:5], C - (P[ω] - μ[ω]) ⟂ Y[ω])
    # Price is equal to consumer's marginal utility
    @constraint(model, [ω = 1:5], P[ω] - (A[ω] - B * Q[ω]) ⟂ Q[ω])
    # Production is equal to consumption
    @constraint(model, [ω = 1:5], Y[ω] - Q[ω] ⟂ P[ω])
    # Production does not exceed capacity
    @constraint(model, [ω = 1:5], x - Y[ω] ⟂ μ[ω])
    return model
end
```

This instance is featuring mixed-complementarity constraints, and as such
is a good demo for MathOptComplements' capabilities.

As a reference, we use the solution returned by the [PATH solver](https://pages.cs.wisc.edu/~ferris/path.html):
```@example tutorial_equilibrium
using PATHSolver
model = equilibrium_model()
JuMP.set_optimizer(model, PATHSolver.Optimizer)
JuMP.optimize!(model)
nothing

```
The solution returned by PATH is:
```@example tutorial_equilibrium
JuMP.value(model[:x]) # production in MWh
```

### Solution with a nonlinear solver

We replace the solver PATH by Ipopt. MathOptComplements takes care of
reformulating the problem automatically with appropriate nonlinear constraints.
```@example tutorial_equilibrium
using MathOptComplements
using Ipopt
model = equilibrium_model()
MathOptComplements.Bridges.add_all_bridges(model)
set_optimizer(model, Ipopt.Optimizer)
JuMP.optimize!(model)
nothing
```
The solution returned by Ipopt is:
```@example tutorial_equilibrium
JuMP.value(model[:x]) # production in MWh
```



