"""
    SplitIntervalBridge{T,G} <: MOI.Bridges.Constraint.AbstractBridge

`SplitIntervalBridge` implements the following reformulation:

  * `[x, y]` in `ComplementsWithSetType{Interval{T}}` into
    `[xp, y]` in `ComplementsWithSetType{GreaterThan{T}}` and
    `[xn, y]` in `ComplementsWithSetType{LessThan{T}}`

with the equality constraint `x == xp + xn`, where `xp` and `xn` are new
variables representing the positive and negative parts of the activity.

The input function can be any `MOI.AbstractVectorFunction` (the first
component `x` may be affine or quadratic); only the second component `y`
must be a variable.

## Source node

`SplitIntervalBridge` supports:

  * [`MathOptInterface.AbstractVectorFunction`](@extref) in
    [`ComplementsWithSetType{MOI.Interval{T}}`](@ref)

## Target nodes

`SplitIntervalBridge` creates:

  * [`MathOptInterface.VectorOfVariables`](@extref) in
    [`ComplementsWithSetType{MOI.GreaterThan{T}}`](@ref)
  * [`MathOptInterface.VectorOfVariables`](@extref) in
    [`ComplementsWithSetType{MOI.LessThan{T}}`](@ref)
  * `G` in [`MathOptInterface.EqualTo`](@extref) (the splitting equality)

where `G` is the scalar function type of the first component.

"""
struct SplitIntervalBridge{T,G<:MOI.AbstractScalarFunction,F<:MOI.AbstractVectorFunction} <:
       MOI.Bridges.Constraint.AbstractBridge
    lower::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.GreaterThan{T}},
    }
    upper::MOI.ConstraintIndex{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.LessThan{T}},
    }
    equality::MOI.ConstraintIndex{G,MOI.EqualTo{T}}
    xp::MOI.VariableIndex
    xn::MOI.VariableIndex
    func::F
    set::ComplementsWithSetType{MOI.Interval{T}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{SplitIntervalBridge{T,G,F}},
    model::MOI.ModelLike,
    func::F,
    set::ComplementsWithSetType{MOI.Interval{T}},
) where {T,G,F<:MOI.AbstractVectorFunction}
    @assert set.dimension == 2
    scalars = MOIU.scalarize(func)
    x_func = scalars[1]  # activity (may be an expression)
    y = scalars[2]       # slack (must be a variable)
    # y must be a single variable
    y_var = if y isa MOI.VariableIndex
        y
    else
        # Extract the variable from a ScalarAffineFunction wrapping a single variable
        @assert length(y.terms) == 1 && isone(y.terms[1].coefficient) && iszero(y.constant)
        y.terms[1].variable
    end
    # Create xp >= 0 and xn <= 0
    xp, _ = MOI.add_constrained_variable(model, MOI.GreaterThan(zero(T)))
    xn, _ = MOI.add_constrained_variable(model, MOI.LessThan(zero(T)))
    # x == xp + xn
    eq_func = MOIU.operate(-, T, x_func, xp)
    eq_func = MOIU.operate!(-, T, eq_func, xn)
    equality = MOIU.normalize_and_add_constraint(model, eq_func, MOI.EqualTo(zero(T)))
    # [xp, y] in ComplementsWithSetType{GreaterThan{T}}
    lower = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([xp, y_var]),
        ComplementsWithSetType{MOI.GreaterThan{T}}(2),
    )
    # [xn, y] in ComplementsWithSetType{LessThan{T}}
    upper = MOI.add_constraint(
        model,
        MOI.VectorOfVariables([xn, y_var]),
        ComplementsWithSetType{MOI.LessThan{T}}(2),
    )
    return SplitIntervalBridge{T,G,typeof(func)}(lower, upper, equality, xp, xn, func, set)
end

function MOI.supports_constraint(
    ::Type{<:SplitIntervalBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{MOI.Interval{T}}},
) where {T}
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:SplitIntervalBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{MOI.Interval{T}}},
) where {T}
    G = MOIU.scalar_type(F)
    # After `operate(-, T, G, VariableIndex)`, the type may promote
    H = MOIU.promote_operation(-, T, G, MOI.VariableIndex)
    return SplitIntervalBridge{T,H,F}
end

MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:SplitIntervalBridge},
) = true

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::SplitIntervalBridge,
    value::AbstractComplementarityRelaxation,
)
    MOI.set(model, attr, bridge.lower, value)
    MOI.set(model, attr, bridge.upper, value)
    return
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:SplitIntervalBridge{T}},
) where {T}
    return Tuple{Type}[(MOI.GreaterThan{T},), (MOI.LessThan{T},)]
end

function MOI.Bridges.added_constraint_types(::Type{<:SplitIntervalBridge{T,G}}) where {T,G}
    return Tuple{Type,Type}[
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.GreaterThan{T}}),
        (MOI.VectorOfVariables, ComplementsWithSetType{MOI.LessThan{T}}),
        (G, MOI.EqualTo{T}),
    ]
end

function MOI.get(::SplitIntervalBridge, ::MOI.NumberOfVariables)::Int64
    return 2
end

function MOI.get(bridge::SplitIntervalBridge, ::MOI.ListOfVariableIndices)
    return [bridge.xp, bridge.xn]
end

# The constrained variables create VariableIndex-in-GreaterThan/LessThan
# constraints that must be reported as part of the bridge.

function MOI.get(
    ::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{MOI.VariableIndex,MOI.GreaterThan{T}},
)::Int64 where {T}
    return 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.GreaterThan{T}},
) where {T}
    return [MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{T}}(bridge.xp.value)]
end

function MOI.get(
    ::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{MOI.VariableIndex,MOI.LessThan{T}},
)::Int64 where {T}
    return 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.LessThan{T}},
) where {T}
    return [MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{T}}(bridge.xn.value)]
end

function MOI.get(
    ::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.GreaterThan{T}},
    },
)::Int64 where {T}
    return 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.GreaterThan{T}},
    },
) where {T}
    return [bridge.lower]
end

function MOI.get(
    ::SplitIntervalBridge{T},
    ::MOI.NumberOfConstraints{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.LessThan{T}},
    },
)::Int64 where {T}
    return 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T},
    ::MOI.ListOfConstraintIndices{
        MOI.VectorOfVariables,
        ComplementsWithSetType{MOI.LessThan{T}},
    },
) where {T}
    return [bridge.upper]
end

function MOI.get(
    ::SplitIntervalBridge{T,G},
    ::MOI.NumberOfConstraints{G,MOI.EqualTo{T}},
)::Int64 where {T,G}
    return 1
end

function MOI.get(
    bridge::SplitIntervalBridge{T,G},
    ::MOI.ListOfConstraintIndices{G,MOI.EqualTo{T}},
) where {T,G}
    return [bridge.equality]
end

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintFunction, bridge::SplitIntervalBridge)
    return bridge.func
end

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet, bridge::SplitIntervalBridge)
    return bridge.set
end

function MOI.delete(model::MOI.ModelLike, bridge::SplitIntervalBridge)
    MOI.delete(model, bridge.lower)
    MOI.delete(model, bridge.upper)
    MOI.delete(model, bridge.equality)
    MOI.delete(model, bridge.xp)
    MOI.delete(model, bridge.xn)
    return
end
