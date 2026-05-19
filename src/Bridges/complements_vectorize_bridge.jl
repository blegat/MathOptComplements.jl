"""
    ComplementsVectorizeBridge{T,F,S,SV} <: MOI.Bridges.Constraint.AbstractBridge

`ComplementsVectorizeBridge` implements the following reformulations:

  * `[x₁, x₂]` in `ComplementsWithSetType{GreaterThan{T}}` into
    `[x₁, x₂ - lb]` in `ComplementsWithSetType{Nonnegatives}`
  * `[x₁, x₂]` in `ComplementsWithSetType{LessThan{T}}` into
    `[x₁, x₂ - ub]` in `ComplementsWithSetType{Nonpositives}`
  * `[x₁, x₂]` in `ComplementsWithSetType{EqualTo{T}}` into
    `[x₁, x₂ - c]` in `ComplementsWithSetType{Zeros}`

where `T` is the coefficient type.

## Source node

`ComplementsVectorizeBridge` supports:

  * [`MathOptInterface.VectorOfVariables`](@extref) in [`ComplementsWithSetType{S}`](@ref)
    where `S` is [`MathOptInterface.GreaterThan`](@extref), [`MathOptInterface.LessThan`](@extref), or
    [`MathOptInterface.EqualTo`](@extref)

## Target nodes

`ComplementsVectorizeBridge` creates:

  * `F` in [`ComplementsWithSetType{SV}`](@ref), where `SV` is
    [`MathOptInterface.Nonnegatives`](@extref), [`MathOptInterface.Nonpositives`](@extref), or
    [`MathOptInterface.Zeros`](@extref) depending on the input set type

"""
struct ComplementsVectorizeBridge{T,F,S,SV} <: MOI.Bridges.Constraint.AbstractBridge
    constraint::MOI.ConstraintIndex{F,ComplementsWithSetType{SV}}
    set_constant::T
end

function _vector_set_type(::Type{<:MOI.GreaterThan})
    return MOI.Nonnegatives
end
function _vector_set_type(::Type{<:MOI.LessThan})
    return MOI.Nonpositives
end
function _vector_set_type(::Type{<:MOI.EqualTo})
    return MOI.Zeros
end

function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.GreaterThan},
    x2,
) where {T}
    return MOIU.get_bounds(model, T, x2)[1]
end
function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.LessThan},
    x2,
) where {T}
    return MOIU.get_bounds(model, T, x2)[2]
end
function _set_constant(
    ::Type{T},
    model,
    ::ComplementsWithSetType{<:MOI.EqualTo},
    x2,
) where {T}
    return MOIU.get_bounds(model, T, x2)[1]
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ComplementsVectorizeBridge{T,F,S,SV}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::ComplementsWithSetType{S},
) where {T,F,S,SV}
    @assert set.dimension == 2
    x1 = func.variables[1]
    x2 = func.variables[2]
    c = _set_constant(T, model, set, x2)
    # Shift only x2: [x1, x2] → [x1, x2 - c]
    shifted = MOIU.operate(vcat, T, one(T) * x1, one(T) * x2 - c)
    ci = MOI.add_constraint(model, shifted, ComplementsWithSetType{SV}(2))
    return ComplementsVectorizeBridge{T,F,S,SV}(ci, c)
end

function MOI.supports_constraint(
    ::Type{<:ComplementsVectorizeBridge},
    ::Type{MOI.VectorOfVariables},
    ::Type{<:ComplementsWithSetType{<:Union{MOI.GreaterThan,MOI.LessThan,MOI.EqualTo}}},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ComplementsVectorizeBridge{T}},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{S}},
) where {T,S<:Union{MOI.GreaterThan,MOI.LessThan,MOI.EqualTo}}
    G = MOI.Utilities.promote_operation(-, T, MOI.VariableIndex, T)
    F = MOI.Utilities.promote_operation(vcat, T, G, G)
    SV = _vector_set_type(S)
    return ComplementsVectorizeBridge{T,F,S,SV}
end

MOI.supports(
    ::MOI.ModelLike,
    ::ComplementarityReformulation,
    ::Type{<:ComplementsVectorizeBridge},
) = true

function MOI.set(
    model::MOI.ModelLike,
    attr::ComplementarityReformulation,
    bridge::ComplementsVectorizeBridge,
    value::AbstractComplementarityRelaxation,
)
    MOI.set(model, attr, bridge.constraint, value)
    return
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:ComplementsVectorizeBridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{ComplementsVectorizeBridge{T,F,S,SV}},
) where {T,F,S,SV}
    return Tuple{Type,Type}[(F, ComplementsWithSetType{SV})]
end

function MOI.get(
    ::ComplementsVectorizeBridge{T,F,S,SV},
    ::MOI.NumberOfConstraints{F,ComplementsWithSetType{SV}},
)::Int64 where {T,F,S,SV}
    return 1
end

function MOI.get(
    bridge::ComplementsVectorizeBridge{T,F,S,SV},
    ::MOI.ListOfConstraintIndices{F,ComplementsWithSetType{SV}},
) where {T,F,S,SV}
    return [bridge.constraint]
end

function MOI.get(
    model::MOI.ModelLike,
    ::MOI.ConstraintFunction,
    bridge::ComplementsVectorizeBridge,
)
    f = MOI.get(model, MOI.ConstraintFunction(), bridge.constraint)
    terms = f.terms
    vars = [t.scalar_term.variable for t in terms]
    return MOI.VectorOfVariables(vars)
end

function MOI.get(
    ::MOI.ModelLike,
    ::MOI.ConstraintSet,
    ::ComplementsVectorizeBridge{T,F,S},
) where {T,F,S}
    return ComplementsWithSetType{S}(2)
end

function MOI.delete(model::MOI.ModelLike, bridge::ComplementsVectorizeBridge)
    MOI.delete(model, bridge.constraint)
    return
end
