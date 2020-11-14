import.. PolaronMakie
using.. PolaronMakie
using QuadGK
using Plots

function ℜχ(Ω, β, α, v, w)
    @show(Ω)
    integrand(x) = (1 - cos(Ω * x)) * ℑS(x, 0, v, w, β, α)
    return QuadGK.quadgk(x -> integrand(x), 0, Inf)[1]
end

function ℑχ(Ω, β, α, v, w)
    @show(Ω)
    integrand(x) = -sin(Ω * x) * ℑS(x, 0, v, w, β, α)
    return QuadGK.quadgk(x -> integrand(x), 0, Inf)[1]
end

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

v = range(0.0, stop = 100.0, length = 1000)
params = (100, 3.0, 3.4, 2.5)
# ReX = [ℜχ(t, params...) for t in v]
ImX = [ℑχ(t, params...) for t in v]
# plot(v, ReX, label = "ℜχ")
plot(v, ImX, label = "ℑχ")
xlabel!("Ω")
ylabel!("χ(Ω)")
