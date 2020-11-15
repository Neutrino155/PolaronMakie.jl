module PolaronMakie

export viewfunction, ℜD, ℑD, D, ℜS, ℑS, S, ℜχ, ℑχ, χ

using QuadGK
using AbstractPlotting.MakieLayout
using AbstractPlotting
using GLMakie
using Plots

include("../src/function.jl")
include("../src/Functions/D(u)_function.jl")
include("../src/Functions/S(u)_function.jl")
include("../src/Functions/X(u)_function.jl")

end # PolaronMakie
