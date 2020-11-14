module PolaronMakie

export viewfunction
export ℑS

using AbstractPlotting.MakieLayout
using AbstractPlotting
using GLMakie
using QuadGK
using Quadmath
using Plots

include("../src/function.jl")
include("../src/Functions/D(u)_function.jl")
include("../src/Functions/S(u)_function.jl")

end # PolaronMakie
