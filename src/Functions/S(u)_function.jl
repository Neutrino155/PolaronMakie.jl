import ..PolaronMakie
using ..PolaronMakie

mod_squared_D(x, y, v, w, β) = (ℜD(x, y, v, w, β))^2 + (ℑD(x, y, v, w, β))^2

arg_D(x, y, v, w, β) = atan(ℑD(x, y, v, w, β), ℜD(x, y, v, w, β))

function ℜS(x, y, v, w, β, α)
    P = 1 / (exp(β) - 1)
    (2 * α / (3 * sqrt(π))) * (exp(-y) * (1 + P) * cos(x - 3 * arg_D(x, y, v, w, β) / 2) + exp(y) * P * cos(x + 3 * arg_D(x, y, v, w, β) / 2)) / (mod_squared_D(x, y, v, w, β))^(3 / 4)
end

function ℑS(x, y, v, w, β, α)
    P = 1 / (exp(β) - 1)
    (2 * α / (3 * sqrt(π))) * (exp(-y) * (1 + P) * sin(x - 3 * arg_D(x, y, v, w, β) / 2) - exp(y) * P * sin(x + 3 * arg_D(x, y, v, w, β) / 2)) / (mod_squared_D(x, y, v, w, β))^(3 / 4)
end

S(x, y, v, w, β, α) = ℜS(x, y, v, w, β, α) + 1im * ℑS(x, y, v, w, β, α)

# p = (:v, :w, :β, :α)
# i = (2.0, 1.0, 10, 6.0)
# pl = ((0.0, 10.0), (0.0, 10.0), (0.0, 10.0), (0.0, 10.0))
# al = ((-25, 25), (-25, 25), (-25, 25))
# a = Dict(
#     "ℜ_title" => "Real S(u)",
#     "ℜ_x" => "ℜ(u)",
#     "ℜ_y" => "ℑ(u)",
#     "ℜ_z" => "ℜ[D(u)]",
#     "ℑ_title" => "Imaginary S(u)",
#     "ℑ_x" => "ℜ(u)",
#     "ℑ_y" => "ℑ(u)",
#     "ℑ_z" => "ℑ[D(u)]",
# )

# vf = viewfunction(ℜS, ℑS, p, i, pl, a, axis_limits = al, len = 100)
