module Bridges

import MathOptInterface as MOI
const MOIU = MOI.Utilities
using ..MathOptComplements:
    ComplementsWithSetType,
    AbstractComplementarityRelaxation,
    ComplementarityReformulation,
    _remove_bounds!

include("vertical.jl")
include("specify_set_type_bridge.jl")
include("complements_vectorize_bridge.jl")
include("split_interval_bridge.jl")
include("flip_sign_bridge.jl")
include("nonlinear.jl")
include("to_sos1_bridge.jl")

"""
    add_all_bridges(model::MOI.ModelLike, ::Type{T} = Float64)

Add all `MathOptComplements` bridges to `model`. The model is typically a
[`MathOptInterface.Bridges.LazyBridgeOptimizer`](@extref) so that the bridge graph is
extended with the bridges needed to reformulate
[`MathOptComplements.ComplementsWithSetType`](@ref) and [`MathOptInterface.Complements`](@extref)
constraints.

When used with a `LazyBridgeOptimizer`, the [`NonlinearBridge`](@ref) uses
the default [`ScholtesRelaxation`](@ref) because the
[`MathOptComplements.DefaultComplementarityReformulation`](@ref) optimizer
attribute is only supported by `MathOptComplements.Optimizer`.
"""
function add_all_bridges(model::MOI.ModelLike, ::Type{T} = Float64) where {T}
    MOI.Bridges.add_bridge(model, SpecifySetTypeBridge{T})
    MOI.Bridges.add_bridge(model, ComplementsVectorizeBridge{T})
    MOI.Bridges.add_bridge(model, SplitIntervalBridge{T})
    MOI.Bridges.add_bridge(model, FlipSignBridge{T})
    MOI.Bridges.add_bridge(model, ToSOS1Bridge{T})
    MOI.Bridges.add_bridge(model, VerticalBridge{T})
    MOI.Bridges.add_bridge(model, NonlinearBridge{T})
    return
end

end # module Bridges
