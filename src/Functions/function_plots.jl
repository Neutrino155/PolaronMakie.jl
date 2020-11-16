# function_plots.jl

"""
Script for creating either 3D Makie plots, or 2D Plotly plots, of functions from the Functions directory.
"""

using ..PolaronMakie
using Plots
using ProgressBars
using Printf
plotly()

"""
D(u) Function
"""

# p = (:v, :w, :β)
# i = (1.0, 1.0, 5.0)
# pl = ((0.0, 10.0), (0.0, 10.0), (0.0, 10.0))
# al = ((-25, 25), (-25, 25), (-25, 25))
# a = Dict(
#     "ℜ_title" => "Real D(u)",
#     "ℜ_x" => "ℜ(u)",
#     "ℜ_y" => "ℑ(u)",
#     "ℜ_z" => "ℜ[D(u)]",
#     "ℑ_title" => "Imaginary D(u)",
#     "ℑ_x" => "ℜ(u)",
#     "ℑ_y" => "ℑ(u)",
#     "ℑ_z" => "ℑ[D(u)]",
# )

# vf = viewfunction(ℜD, ℑD, p, i, pl, a, axis_limits = al, len = 200)

"""
S(u) Function
"""

# p = (:v, :w, :β, :α)
# i = (2.0, 1.0, 10, 6.0)
# pl = ((0.0, 10.0), (0.0, 10.0), (0.0, 10.0), (0.0, 10.0))
# al = ((-25, 25), (-25, 25), (-25, 25))
# a = Dict(
#     "ℜ_title" => "Real S(u)",
#     "ℜ_x" => "ℜ(u)",
#     "ℜ_y" => "ℑ(u)",
#     "ℜ_z" => "ℜ[S(u)]",
#     "ℑ_title" => "Imaginary S(u)",
#     "ℑ_x" => "ℜ(u)",
#     "ℑ_y" => "ℑ(u)",
#     "ℑ_z" => "ℑ[S(u)]",
# )

# vf = viewfunction(ℜS, ℑS, p, i, pl, a, axis_limits = al, len = 100)

"""
X(u) Function
"""

# p = (:α, :v, :w)
# i = (5, 4.0, 2.1)
# pl = ((0.01, 10), (0.01, 10), (0.01, 10))
# al = ((0.01, 10), (0.01, 10), (-50, 100))
# a = Dict(
#     "ℜ_title" => "Real Χ(u)",
#     "ℜ_x" => "ℜ(u)",
#     "ℜ_y" => "ℑ(u)",
#     "ℜ_z" => "ℜ[X(u)]",
#     "ℑ_title" => "Imaginary Χ(u)",
#     "ℑ_x" => "ℜ(u)",
#     "ℑ_y" => "ℑ(u)",
#     "ℑ_z" => "ℑ[X(u)]",
# )
#
# vf = viewfunction(ℜχ, ℑχ, p, i, pl, a, axis_limits = al, len = 200)

Ω = range(0.0, stop = 10.0, length = 100)
params = (3.0, 7.0, 5.8, 1.6)
ReX = Array{Float64}(undef, length(Ω))
ImX = Array{Float64}(undef, length(Ω))
iter = ProgressBar(range(1, stop = length(Ω), step = 1))
Threads.@threads for i in iter
    ReX[i] = ℜχ(Ω[i], params...)
    ImX[i] = ℑχ(Ω[i], params...)
    set_description(iter, string(@sprintf("Ω: %.2f | ", Ω[i])) * string(@sprintf("ℜχ: %.2f | ", ℜχ(Ω[i], params...))) * string(@sprintf("ℑχ: %.2f |", ℑχ(Ω[i], params...))))
end
ReX_plot = Plots.plot(Ω, ReX, label = "ℜχ")
ImX_plot = Plots.plot!(Ω, ImX, label = "ℑχ")
Plots.xlabel!("Ω")
Plots.ylabel!("χ(Ω)")
