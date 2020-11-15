using ..PolaronMakie

function ℜχ(Ω, β, α, v, w)
    integrand(x) = (1 - cos(Ω * x)) * ℑS(x, 0, v, w, β, α)
    return QuadGK.quadgk(x -> integrand(x), 0, Inf)[1]
end

function ℑχ(Ω, β, α, v, w)
    integrand(x) = -sin(Ω * x) * ℑS(x, 0, v, w, β, α)
    return QuadGK.quadgk(x -> integrand(x), 0, Inf)[1]
end

χ(Ω, β, α, v, w) = ℜχ(Ω, β, α, v, w) + 1im * ℑχ(Ω, β, α, v, w)

# p = (:α, :v, :w)
# i = (5, 4.0, 2.1)
# pl = ((0.01, 10), (0.01, 10), (0.01, 10))
# al = ((0.01, 10), (0.01, 10), (-50, 100))
# a = Dict(
#     "ℜ_title" => "Real Χ(u)",
#     "ℜ_x" => "ℜ(u)",
#     "ℜ_y" => "ℑ(u)",
#     "ℜ_z" => "ℜ[D(u)]",
#     "ℑ_title" => "Imaginary Χ(u)",
#     "ℑ_x" => "ℜ(u)",
#     "ℑ_y" => "ℑ(u)",
#     "ℑ_z" => "ℑ[D(u)]",
# )
#
# vf = viewfunction(ℜχ, ℑχ, p, i, pl, a, axis_limits = al, len = 200)

v = range(0.0, stop = 100.0, length = 1000)
params = (100, 3.0, 3.4, 2.5)
ReX = [ℜχ(t, params...) for t in v]
ImX = [ℑχ(t, params...) for t in v]
plot(v, ReX, label = "ℜχ")
plot(v, ImX, label = "ℑχ")
xlabel!("Ω")
ylabel!("χ(Ω)")
