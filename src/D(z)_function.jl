import .PolaronMakie
using .PolaronMakie

function ℜD(x, y, v, w, β, s)
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    R * cos(v * x) * (exp(-v * s*y) + 2 * P * cosh(v * y)) - (x^2 - y^2) / β - y -
    R * (2 * P + 1)
end

function ℑD(x, y, v, w, β, s)
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    R * sin(v * s * x) * (exp(-v * y) - 2 * P * sinh(v * y)) - x * (2 * y / β - 1)
end

p = (:v, :w, :β, :s)
pl = (10, 10, 20, 15.7)
al = (25, 25, 25)
a = Dict(
    "ℜ_title" => "Real D(u)",
    "ℜ_x" => "ℜ(u)",
    "ℜ_y" => "ℑ(u)",
    "ℜ_z" => "ℜ[D(u)]",
    "ℑ_title" => "Imaginary D(u)",
    "ℑ_x" => "ℜ(u)",
    "ℑ_y" => "ℑ(u)",
    "ℑ_z" => "ℑ[D(u)]",
)

vf = viewfunction(ℜD, ℑD, p, pl, a, axis_limits = al, len = 200)
