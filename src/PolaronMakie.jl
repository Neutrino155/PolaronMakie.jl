module PolaronMakie

export viewfunction
export ℜD
export ℑD
export ℜS
export ℑS

using AbstractPlotting.MakieLayout
using AbstractPlotting
using GLMakie
using QuadGK

include("../src/function.jl")
include("../src/Functions/D(u)_function.jl")
include("../src/Functions/S(u)_function.jl")

end # PolaronMakie
