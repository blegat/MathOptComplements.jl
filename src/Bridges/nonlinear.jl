"""
    NonlinearBridge{T,S} <: MOI.Bridges.Constraint.AbstractBridge

`NonlinearBridge` implements the following reformulation:

  * `[x₁, x₂]` in [`ComplementsWithSetType{S}`](@ref) into nonlinear
    constraints using a complementarity relaxation method

The relaxation method is determined by the
[`AbstractComplementarityRelaxation`](@ref) attribute (default:
[`ScholtesRelaxation`](@ref)). The set type `S` determines the bound case:

  * `S <: Union{MOI.Nonnegatives, MOI.GreaterThan}`: lower-bound relaxation
  * `S <: Union{MOI.Nonpositives, MOI.LessThan}`: upper-bound relaxation
  * `S <: MOI.Interval`: range relaxation
  * `S == MOI.Zeros`: equality (trivial complementarity)

## Source node

`NonlinearBridge` supports:

  * [`MathOptInterface.VectorOfVariables`](@extref) in [`ComplementsWithSetType{S}`](@ref)

## Target nodes

`NonlinearBridge` creates nonlinear constraints depending on the relaxation
method (for example, quadratic or `ScalarNonlinearFunction` constraints).

"""
mutable struct NonlinearBridge{T,S} <: MOI.Bridges.Constraint.AbstractBridge
    constraints::Vector
    func::MOI.VectorOfVariables
    set::ComplementsWithSetType{S}
    reformulation::AbstractComplementarityRelaxation
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{NonlinearBridge{T,S}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::ComplementsWithSetType{S},
    # Only `MathOptComplements.Optimizer` supports setting custom bridge arguments
    # so the default will be used when the bridge is added for instance to
    # `MOI.Bridges.LazyBridgeOptimizer`
    reformulation::AbstractComplementarityRelaxation = ScholtesRelaxation(zero(T)),
) where {T,S}
    # Delay reformulation until `final_touch` so that per-constraint
    # `ComplementarityReformulation` attributes can override it first.
    return NonlinearBridge{T,S}([], func, set, reformulation)
end

function MOI.supports_constraint(
    ::Type{<:NonlinearBridge},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:ComplementsWithSetType},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:NonlinearBridge{T}},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {T,S}
    return NonlinearBridge{T,S}
end

MOI.supports(::MOI.ModelLike, ::ComplementarityReformulation, ::Type{<:NonlinearBridge}) =
    true

function MOI.set(
    model::MOI.ModelLike,
    ::ComplementarityReformulation,
    bridge::NonlinearBridge,
    value::AbstractComplementarityRelaxation,
)
    bridge.reformulation = value
    # Delete old constraints so that final_touch re-reformulates
    for ci in bridge.constraints
        MOI.delete(model, ci)
    end
    empty!(bridge.constraints)
    return
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:NonlinearBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(::Type{<:NonlinearBridge})
    return Tuple{Type,Type}[]
end

function MOI.get(bridge::NonlinearBridge, ::MOI.NumberOfConstraints{F,S})::Int64 where {F,S}
    return count(ci -> ci isa MOI.ConstraintIndex{F,S}, bridge.constraints)
end

function MOI.get(bridge::NonlinearBridge, ::MOI.ListOfConstraintIndices{F,S}) where {F,S}
    return MOI.ConstraintIndex{F,S}[
        ci for ci in bridge.constraints if ci isa MOI.ConstraintIndex{F,S}
    ]
end

function MOI.delete(model::MOI.ModelLike, bridge::NonlinearBridge)
    for ci in bridge.constraints
        MOI.delete(model, ci)
    end
    return
end

MOI.Bridges.needs_final_touch(::NonlinearBridge) = true

function MOI.Bridges.final_touch(bridge::NonlinearBridge, model::MOI.ModelLike)
    if !isempty(bridge.constraints)
        return
    end
    append!(
        bridge.constraints,
        reformulate_as_nonlinear_program!(
            model,
            bridge.reformulation,
            bridge.func,
            bridge.set,
        ),
    )
    return
end

# Bound helpers: extract the relevant bounds based on the set type S.
_complementarity_bounds(::Type{MOI.Nonnegatives}, model, ::Type{T}, x2) where {T} =
    (zero(T), T(Inf))
_complementarity_bounds(::Type{MOI.Nonpositives}, model, ::Type{T}, x2) where {T} =
    (T(-Inf), zero(T))
_complementarity_bounds(::Type{MOI.Zeros}, model, ::Type{T}, x2) where {T} =
    (zero(T), zero(T))
function _complementarity_bounds(::Type{<:MOI.GreaterThan}, model, ::Type{T}, x2) where {T}
    return (MOIU.get_bounds(model, T, x2)[1], T(Inf))
end
function _complementarity_bounds(::Type{<:MOI.LessThan}, model, ::Type{T}, x2) where {T}
    return (T(-Inf), MOIU.get_bounds(model, T, x2)[2])
end
function _complementarity_bounds(::Type{<:MOI.Interval}, model, ::Type{T}, x2) where {T}
    return MOIU.get_bounds(model, T, x2)
end

"""
    reformulate_as_nonlinear_program!(model, relaxation, fun, set::ComplementsWithSetType{S})

Reformulate complementarity constraints as a nonlinear program using the given
relaxation. The set type `S` determines which bound case to use.

"""
function reformulate_as_nonlinear_program!(
    model::MOI.ModelLike,
    relaxation::AbstractComplementarityRelaxation,
    fun,
    set::ComplementsWithSetType{S},
) where {S}
    n_comp = div(set.dimension, 2)
    ind_cc = []
    for cc in 1:n_comp
        x1 = fun.variables[cc]
        x2 = fun.variables[cc+n_comp]
        lb2, ub2 = _complementarity_bounds(S, model, Float64, x2)
        if isinf(ub2)
            idc = _relax_complementarity_lower_bound!(model, relaxation, x1, x2, lb2, ub2)
        elseif isinf(lb2)
            idc = _relax_complementarity_upper_bound!(model, relaxation, x1, x2, lb2, ub2)
        else
            idc = _relax_complementarity_range!(model, relaxation, x1, x2, lb2, ub2)
        end
        append!(ind_cc, idc)
    end
    return ind_cc
end


"""
    ScholtesRelaxation <: AbstractComplementarityRelaxation

Implement the Scholtes relaxation.
For `tau ≥ 0`, the complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
0 ≤ a
0 ≤ b
a . b ≤ tau
```

If `tau` is equal to 0, then the relaxation is exact.

"""
struct ScholtesRelaxation{T} <: AbstractComplementarityRelaxation
    tau::T
end

# x1 ⟂ (lb <= x2)   ≡  0 <= x1 ; lb <= x2 ; x1 ( x2 - lb) <= 0
function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    relaxation::ScholtesRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, _ = MOIU.get_bounds(model, Float64, x1)
    if isinf(lb1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert lb1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if ub1 is finite?
    end
    idc = MOI.add_constraint(model, x1 * (x2 - lb2), MOI.LessThan(relaxation.tau))
    return [idc]
end

# x1 ⟂ (x2 <= ub)   ≡  x1 <= 0 ; lb <= x2 ; x1 ( x2 - ub) <= 0
function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    relaxation::ScholtesRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    _, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(ub1)
        MOI.add_constraint(model, x1, MOI.LessThan(0.0))
    else
        @assert ub1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if lb1 is finite?
    end
    idc = MOI.add_constraint(model, x1 * (x2 - ub2), MOI.LessThan(relaxation.tau))
    return [idc]
end

# x1 ⟂ (lb <= x2 <= ub) ≡  lb <= x2 <= ub ; x1 (x2 - lb)  <= 0 ; x1 (x2 - ub) <= 0
function _relax_complementarity_range!(
    model::MOI.ModelLike,
    relaxation::ScholtesRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    idc1 = MOI.add_constraint(model, x1 * (x2 - lb2), MOI.LessThan(relaxation.tau))
    idc2 = MOI.add_constraint(model, x1 * (x2 - ub2), MOI.LessThan(relaxation.tau))
    return [idc1, idc2]
end

"""
    FischerBurmeisterRelaxation <: AbstractComplementarityRelaxation

Implement the Fischer-Burmeister relaxation for complementarity constraints.
For `epsilon ≥ 0`, the complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
0 ≤ a
0 ≤ b
a + b - sqrt((a + b)^2 + epsilon) ≤ 0
```

"""
struct FischerBurmeisterRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _min_eps(a, b, eps)
    return MOI.ScalarNonlinearFunction(
        :-,
        Any[
            1.0*a+1.0*b,
            MOI.ScalarNonlinearFunction(:sqrt, Any[(1.0*a)^2+(1.0*b)^2+eps^2]),
        ],
    )
end

function _max_eps(a, b, eps)
    return MOI.ScalarNonlinearFunction(
        :+,
        Any[
            1.0*a+1.0*b,
            MOI.ScalarNonlinearFunction(:sqrt, Any[(1.0*a)^2+(1.0*b)^2+eps^2]),
        ],
    )
end

# x1 ⟂ (lb <= x2)   ≡  0 <= x1 ; lb <= x2 ; min(x1, x2 - lb) <= 0
function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    relaxation::FischerBurmeisterRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(lb1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert lb1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if ub1 is finite?
    end
    idc = MOI.add_constraint(
        model,
        _min_eps(x1, x2 - lb2, relaxation.epsilon),
        MOI.LessThan(0.0),
    )
    return [idc]
end

# x1 ⟂ (x2 <= ub)   ≡  x1 <= 0 ; lb <= x2 ; max(x1, x2 - ub) >= 0
function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    relaxation::FischerBurmeisterRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(ub1)
        MOI.add_constraint(model, x1, MOI.LessThan(0.0))
    else
        @assert ub1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if lb1 is finite?
    end
    idc = MOI.add_constraint(
        model,
        _max_eps(x1, x2 - ub2, relaxation.epsilon),
        MOI.GreaterThan(0.0),
    )
    return [idc]
end

# x1 ⟂ (lb <= x2 <= ub) ≡  lb <= x2 <= ub ; min(x1, x2 - lb) <= 0 ; max(x1, x2 - ub) >= 0
function _relax_complementarity_range!(
    model::MOI.ModelLike,
    relaxation::FischerBurmeisterRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    idc1 = MOI.add_constraint(
        model,
        _min_eps(x1, x2 - lb2, relaxation.epsilon),
        MOI.LessThan(0.0),
    )
    idc2 = MOI.add_constraint(
        model,
        _max_eps(x1, x2 - ub2, relaxation.epsilon),
        MOI.GreaterThan(0.0),
    )
    return [idc1, idc2]
end


"""
    LiuFukushimaRelaxation <: AbstractComplementarityRelaxation

Implement the Liu-Fukushima relaxation for complementarity constraints.
For `epsilon ≥ 0`, the complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
a . b ≤ epsilon^2
(a + epsilon) . (b + epsilon) ≥ epsilon^2

```
"""
struct LiuFukushimaRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    relaxation::LiuFukushimaRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    # Ensure we respect MOI's specs
    lb1, _ = MOIU.get_bounds(model, Float64, x1)
    @assert isinf(lb1) || iszero(lb1)

    idc1 = MOI.add_constraint(model, x1 * (x2 - lb2), MOI.LessThan(relaxation.epsilon^2))
    idc2 = MOI.add_constraint(
        model,
        (x1 + relaxation.epsilon) * (x2 - lb2 + relaxation.epsilon),
        MOI.GreaterThan(relaxation.epsilon^2),
    )

    # Remove bounds
    _remove_bounds!(model, x1)
    cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{Float64}}(x2.value)
    MOI.delete(model, cidx)

    return [idc1, idc2]
end

function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    relaxation::LiuFukushimaRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    # Ensure we respect MOI's specs
    _, ub1 = MOIU.get_bounds(model, Float64, x1)
    @assert isinf(ub1) || iszero(ub1)

    idc1 = MOI.add_constraint(model, x1 * (x2 - ub2), MOI.LessThan(relaxation.epsilon^2))
    idc2 = MOI.add_constraint(
        model,
        (x1 - relaxation.epsilon) * (x2 - ub2 - relaxation.epsilon),
        MOI.GreaterThan(relaxation.epsilon^2),
    )

    # Remove bounds
    _remove_bounds!(model, x1)
    cidx = MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{Float64}}(x2.value)
    MOI.delete(model, cidx)

    return [idc1, idc2]
end


"""
    KanzowSchwarzRelaxation <: AbstractComplementarityRelaxation

Implement the Kanzow-Schwarz relaxation for complementarity constraints.
For `epsilon ≥ 0`, the complementarity constraint `0 ≤ a ⟂ b ≥ 0` is reformulated as
```
0 ≤ a
0 ≤ b
ϕ(a, b) ≤ 0

```
with the function `ϕ`:
```
ϕ(a, b) = (a - epsilon) . (b - epsilon)                 if a + b > 2 epsilon
          -0.5 ((a -epsilon)^2 + (b - epsilon)^2)       otherwise

```

"""
struct KanzowSchwarzRelaxation{T} <: AbstractComplementarityRelaxation
    epsilon::T
end

function _kanzow_schwarz_relaxation(a, b, eps)
    return MOI.ScalarNonlinearFunction(
        :ifelse,
        Any[
            MOI.ScalarNonlinearFunction(:>, Any[1.0*a+1.0*b, 2*eps]),
            (a-eps)*(b-eps),
            - 0.5*((a-eps)^2+(b-eps)^2),
        ],
    )
end

function _relax_complementarity_lower_bound!(
    model::MOI.ModelLike,
    relaxation::KanzowSchwarzRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(lb1)
        MOI.add_constraint(model, x1, MOI.GreaterThan(0.0))
    else
        @assert lb1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if ub1 is finite?
    end
    idc = MOI.add_constraint(
        model,
        _kanzow_schwarz_relaxation(x1, x2 - lb2, relaxation.epsilon),
        MOI.LessThan(0.0),
    )
    return [idc]
end

# x1 ⟂ (x2 <= ub)   ≡  x1 <= 0 ; lb <= x2 ; max(x1, x2 - ub) >= 0
function _relax_complementarity_upper_bound!(
    model::MOI.ModelLike,
    relaxation::KanzowSchwarzRelaxation,
    x1,
    x2,
    lb2,
    ub2,
)
    lb1, ub1 = MOIU.get_bounds(model, Float64, x1)
    if isinf(ub1)
        MOI.add_constraint(model, x1, MOI.LessThan(0.0))
    else
        @assert ub1 == 0.0 # ensure we follow MOI's convention
        # TODO: what should we do if lb1 is finite?
    end
    idc = MOI.add_constraint(
        model,
        _kanzow_schwarz_relaxation(-1.0*x1, ub2 - x2, relaxation.epsilon),
        MOI.LessThan(0.0),
    )
    return [idc]
end
