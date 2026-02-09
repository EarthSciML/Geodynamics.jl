module Geodynamics

using DocStringExtensions
using ModelingToolkit: t, D, System, @variables, @parameters, @named,
    @constants, @component
using DynamicQuantities: @u_str
using EarthSciMLBase

include("mogi_mctigue.jl")

end # module
