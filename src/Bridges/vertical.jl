"""
    VerticalBridge{T,S} <: MOI.Bridges.Constraint.AbstractBridge

`VerticalBridge` implements the following reformulation:

  * `f(x)` in `S` into `[x₁, x₂]` in `S`

where expression-based complementarity constraints are converted to vertical
form by introducing slack variables. If the left-hand side is an expression,
a slack variable `x₁` is created with an equality `lhs = x₁`. If the
right-hand side variable `x₂` is unbounded, the left-hand side is converted
to an equality constraint instead.

## Source node

`VerticalBridge` supports:

  * [`MathOptInterface.AbstractVectorFunction`](@extref) in [`MathOptInterface.Complements`](@extref)
  * [`MathOptInterface.AbstractVectorFunction`](@extref) in [`ComplementsWithSetType{S}`](@ref)

## Target nodes

`VerticalBridge` creates:

  * [`MathOptInterface.VectorOfVariables`](@extref) in `S`
  * [`MathOptInterface.ScalarAffineFunction`](@extref) in [`MathOptInterface.EqualTo`](@extref)

"""
struct VerticalBridge{T,S<:MOI.AbstractVectorSet} <: MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex{MOI.VectorOfVariables,S}
    equalities::Vector{MOI.ConstraintIndex}
    slacks::Vector{MOI.VariableIndex}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{VerticalBridge{T,MOI.Complements}},
    model::MOI.ModelLike,
    func::MOI.AbstractVectorFunction,
    set::MOI.Complements,
) where {T}
    ci, equalities, slacks = reformulate_to_vertical!(model, T, func, set)
    return VerticalBridge{T,MOI.Complements}(ci, equalities, slacks)
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{VerticalBridge{T,ComplementsWithSetType{S}}},
    model::MOI.ModelLike,
    func::MOI.AbstractVectorFunction,
    set::ComplementsWithSetType{S},
) where {T,S}
    ci, equalities, slacks = reformulate_to_vertical!(model, T, func, set)
    return VerticalBridge{T,ComplementsWithSetType{S}}(ci, equalities, slacks)
end

function MOI.supports_constraint(
    ::Type{<:VerticalBridge},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
)
    return true
end

function MOI.supports_constraint(
    ::Type{<:VerticalBridge},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ComplementsWithSetType},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:VerticalBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{MOI.Complements},
) where {T}
    return VerticalBridge{T,MOI.Complements}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:VerticalBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{ComplementsWithSetType{S}},
) where {T,S}
    return VerticalBridge{T,ComplementsWithSetType{S}}
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:VerticalBridge})
    return Tuple{Type}[(MOI.Reals,)]
end

function MOI.Bridges.added_constraint_types(::Type{VerticalBridge{T,S}}) where {T,S}
    return Tuple{Type,Type}[
        (MOI.VectorOfVariables, S),
        (MOI.ScalarAffineFunction{T}, MOI.EqualTo{T}),
    ]
end

function MOI.get(bridge::VerticalBridge, ::MOI.NumberOfVariables)::Int64
    return length(bridge.slacks)
end

function MOI.get(bridge::VerticalBridge, ::MOI.ListOfVariableIndices)
    return copy(bridge.slacks)
end

function MOI.get(
    ::VerticalBridge{T,S},
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,S},
)::Int64 where {T,S}
    return 1
end

function MOI.get(
    bridge::VerticalBridge{T,S},
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,S},
) where {T,S}
    return [bridge.constraint]
end

function MOI.get(
    bridge::VerticalBridge{T},
    ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}},
)::Int64 where {T}
    return count(
        ci -> ci isa MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}},
        bridge.equalities,
    )
end

function MOI.get(
    bridge::VerticalBridge{T},
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}},
) where {T}
    return MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}}[
        ci for ci in bridge.equalities if
        ci isa MOI.ConstraintIndex{MOI.ScalarAffineFunction{T},MOI.EqualTo{T}}
    ]
end

function MOI.delete(model::MOI.ModelLike, bridge::VerticalBridge)
    MOI.delete(model, bridge.constraint)
    for ci in bridge.equalities
        MOI.delete(model, ci)
    end
    for vi in bridge.slacks
        MOI.delete(model, vi)
    end
    return
end

MOI.supports(::MOI.ModelLike, ::ComplementarityReformulation, ::Type{<:VerticalBridge}) =
    true

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::VerticalBridge,
    value::AbstractComplementarityRelaxation,
)
    MOI.set(model, attr, bridge.constraint, value)
    return
end

#=
    Parser for JuMP problems with complementarity constraints.
=#

function _is_single_variable(func::MOI.ScalarAffineFunction)
    return length(func.terms) == 1 &&
           func.terms[1].coefficient == 1.0 &&
           iszero(func.constant)
end
function _is_single_variable(func::MOI.ScalarQuadraticFunction)
    return (
        length(func.quadratic_terms) == 0 &&
        length(func.affine_terms) == 1 &&
        func.affine_terms[1].coefficient == 1.0 &&
        iszero(func.constant)
    )
end
function _is_single_variable(func::MOI.ScalarNonlinearFunction)
    return func.head == :+ && length(func.args) == 1 && isa(func.args[1], MOI.VariableIndex)
end
_get_variable(func::MOI.ScalarAffineFunction) = func.terms[1].variable
_get_variable(func::MOI.ScalarQuadraticFunction) = func.affine_terms[1].variable
_get_variable(func::MOI.ScalarNonlinearFunction) = func.args[1]


# TODO: add support for ScalarNonlinearTerm
function _parse_complementarity_constraint(fun::MOI.AbstractVectorFunction, n_comp)
    exprs = MOIU.scalarize(fun)
    @assert length(exprs) == 2*n_comp

    cc_lhs = MOI.AbstractScalarFunction[]
    cc_rhs = MOI.VariableIndex[]

    for i in 1:n_comp
        # Parse LHS
        t1 = exprs[i]
        t2 = exprs[i+n_comp]
        if _is_single_variable(t1)
            push!(cc_lhs, _get_variable(t1))
        else
            push!(cc_lhs, t1)
        end

        # Parse RHS
        isvar2 = _is_single_variable(t2)
        if !isvar2
            # The RHS should be a variable if we follow MOI's specs
            # TODO: we should decide if we should add support complementarity
            # between expressions (see Issue #2)
            error(
                "Right-hand-side should be a single variable in complementarity constraints.",
            )
        end
        push!(cc_rhs, _get_variable(t2))
    end

    return cc_lhs, cc_rhs
end

"""
    reformulate_to_vertical!(model::MOI.ModelLike, ::Type{T}, fun, set)

Factorize all the complementarity constraints in `model` and formulate
an equivalent model in vertical form. The complementarity constraints involving
expressions are rewritten with a slack. `T` is the coefficient type used for
the generated equality constraints.

Once reformulated, the complementarity constraints involve only single variables.

"""
function reformulate_to_vertical!(model::MOI.ModelLike, ::Type{T}, fun, set) where {T}
    equalities = MOI.ConstraintIndex[]
    slacks = MOI.VariableIndex[]
    ind_cc1, ind_cc2 = MOI.VariableIndex[], MOI.VariableIndex[]
    n_comp = div(set.dimension, 2)
    @assert !(fun isa MOI.VectorOfVariables)
    # Read each complementarity constraint and get corresponding indices
    cc_lhs, cc_rhs = _parse_complementarity_constraint(fun, n_comp)
    for (lhs, x2) in zip(cc_lhs, cc_rhs)
        if set isa MOI.Complements
            # Check if x2 is bounded.
            lb, ub = MOIU.get_bounds(model, T, x2)
            if isinf(lb) && isinf(ub)
                # If x2 is unbounded, the LHS is directly converted to an equality constraint.
                push!(
                    equalities,
                    MOIU.normalize_and_add_constraint(model, lhs, MOI.EqualTo{T}(zero(T))),
                )
                continue
            end
        end
        if isa(lhs, MOI.VariableIndex)
            # If lhs is a variable, no need to reformulate the
            # complementarity constraint using a slack.
            # TODO: we should check if the variable lhs is bounded.
            push!(ind_cc1, lhs)
            push!(ind_cc2, x2)
        else
            # Else, reformulate LHS using vertical form
            x1 = MOI.add_variable(model)
            push!(slacks, x1)
            new_lhs = MOIU.operate!(-, T, lhs, x1)
            push!(
                equalities,
                MOIU.normalize_and_add_constraint(model, new_lhs, MOI.EqualTo{T}(zero(T))),
            )
            push!(ind_cc1, x1)
            push!(ind_cc2, x2)
        end
    end
    n_cc = length(ind_cc1)
    comp = MOI.VectorOfVariables([ind_cc1; ind_cc2])
    S = typeof(set)
    if set isa MOI.Complements
        ci = MOI.add_constraint(model, comp, MOI.Complements(2*n_cc))
    else
        ci = MOI.add_constraint(model, comp, S(2*n_cc))
    end
    return ci, equalities, slacks
end
