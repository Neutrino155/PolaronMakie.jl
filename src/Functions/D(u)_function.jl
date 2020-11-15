using ..PolaronMakie

function ℜD(x, y, v, w, β)
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    w^2 / v^2 * (
        R * cos(v * x) * (exp(-v * y) + 2 * P * cosh(v * y)) - (x^2 - y^2) / β - y -
        R * (2 * P + 1)
    )
end

function ℑD(x, y, v, w, β)
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    w^2 / v^2 * (R * sin(v * x) * (exp(-v * y) - 2 * P * sinh(v * y)) - x * (2 * y / β - 1))
end

D(x, y, v, w, β) = ℜD(x, y, v, w, β) + 1im * ℑD(x, y, v, w, β)

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
