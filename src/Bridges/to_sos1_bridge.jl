"""
    ToSOS1Bridge{T} <: MOI.Bridges.Constraint.AbstractBridge

`ToSOS1Bridge` implements the following reformulation:

  * `[x₁, x₂]` in `ComplementsWithSetType{Nonnegatives}` into
    `[x₁, x₂]` in `MOI.SOS1`

Since both `x₁ ≥ 0` and `x₂ ≥ 0`, the SOS1 constraint enforces that at most
one of them is nonzero, which is equivalent to the complementarity condition.

## Source node

`ToSOS1Bridge` supports:

  * [`MathOptInterface.VectorOfVariables`](@extref) in
    [`ComplementsWithSetType{MOI.Nonnegatives}`](@ref)

## Target nodes

`ToSOS1Bridge` creates:

  * [`MathOptInterface.VectorOfVariables`](@extref) in [`MathOptInterface.SOS1`](@extref)

"""
struct ToSOS1Bridge{T} <: MOI.Bridges.Constraint.AbstractBridge
    sos1::MOI.ConstraintIndex{MOI.VectorOfVariables,MOI.SOS1{T}}
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ToSOS1Bridge{T}},
    model::MOI.ModelLike,
    func::MOI.VectorOfVariables,
    set::ComplementsWithSetType{MOI.Nonnegatives},
) where {T}
    @assert set.dimension == 2
    ci = MOI.add_constraint(model, func, MOI.SOS1(T[1.0, 2.0]))
    return ToSOS1Bridge{T}(ci)
end

function MOI.supports_constraint(
    ::Type{<:ToSOS1Bridge},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{MOI.Nonnegatives}},
)
    return true
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ToSOS1Bridge{T}},
    ::Type{MOI.VectorOfVariables},
    ::Type{ComplementsWithSetType{MOI.Nonnegatives}},
) where {T}
    return ToSOS1Bridge{T}
end

# Bridge metadata

function MOI.Bridges.added_constrained_variable_types(::Type{<:ToSOS1Bridge})
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(::Type{ToSOS1Bridge{T}}) where {T}
    return Tuple{Type,Type}[(MOI.VectorOfVariables, MOI.SOS1{T})]
end

function MOI.get(
    ::ToSOS1Bridge{T},
    ::MOI.NumberOfConstraints{MOI.VectorOfVariables,MOI.SOS1{T}},
)::Int64 where {T}
    return 1
end

function MOI.get(
    bridge::ToSOS1Bridge{T},
    ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables,MOI.SOS1{T}},
) where {T}
    return [bridge.sos1]
end

function MOI.get(model::MOI.ModelLike, ::MOI.ConstraintFunction, bridge::ToSOS1Bridge)
    return MOI.get(model, MOI.ConstraintFunction(), bridge.sos1)
end

function MOI.get(::MOI.ModelLike, ::MOI.ConstraintSet, ::ToSOS1Bridge)
    return ComplementsWithSetType{MOI.Nonnegatives}(2)
end

function MOI.delete(model::MOI.ModelLike, bridge::ToSOS1Bridge)
    MOI.delete(model, bridge.sos1)
    return
end
