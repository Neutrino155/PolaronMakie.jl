import ..PolaronMakie
using ..PolaronMakie
using QuadGK

ℜΧ_integrand(x, y, v, w, β, α, ν, t) = ((1 - cos(ν * t / (1 - t)) * exp(-ν * y)) * ℑS(t / (1 - t), y, v, w, β, α) + sin(ν * x) * exp(-ν * t / (1 - t)) * ℑS(x, t / (1 - t), v, w, β, α)) * (1 - t)^(-2)

ℑΧ_integrand(x, y, v, w, β, α, ν, t) = ((1 - cos(ν * x) * exp(-ν * t / (1 - t))) * ℑS(x, t / (1 - t), v, w, β, α) - sin(ν * t / (1 - t)) * exp(-ν * y) * ℑS(t / (1 - t), y, v, w, β, α)) * (1 - t)^(-2)

ℜΧ(x, y, v, w, β, α, ν) = quadgk(t -> ℜΧ_integrand(x, y, v, w, β, α, ν, t), 0, 0.99)[1]

ℑΧ(x, y, v, w, β, α, ν) = quadgk(t -> ℑΧ_integrand(x, y, v, w, β, α, ν, t), 0, 0.9)[1]

p = (:v, :w, :β, :α, :ν)
i = (5.8, 1.6, 10, 7.0, 1.0)
pl = (10, 10, 100, 10, 10)
al = (20, 20, 200)
a = Dict(
    "ℜ_title" => "Real Χ(u)",
    "ℜ_x" => "ℜ(u)",
    "ℜ_y" => "ℑ(u)",
    "ℜ_z" => "ℜ[D(u)]",
    "ℑ_title" => "Imaginary Χ(u)",
    "ℑ_x" => "ℜ(u)",
    "ℑ_y" => "ℑ(u)",
    "ℑ_z" => "ℑ[D(u)]",
)

vf = viewfunction(ℜΧ, ℑΧ, p, i, pl, a, axis_limits = al, len = 100)
