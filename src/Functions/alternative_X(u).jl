# import ..PolaronMakie
# using ..PolaronMakie

"""
This is a playground for testing implementations of the X(u) function. Still need to find a detailed expansions of ReX(u) that is not horrific numerically, if such a thing exists!
"""

using BigCombinatorics
using Combinatorics
using SpecialFunctions: gamma, beta
using ArbNumerics, Readables
using Plots
using Struve
using MathLink
using QuadGK

function L(ν, z)
    ν = BigFloat(ν)
    z = BigFloat(z)
    coefficient = (z / 2)^(ν + 1)
    total_sum = 0.0
    k = 0
    prev_sum = 0.0
    while true
        next_sum = (z / 2)^(2 * k) / (gamma(3/2 + k) * gamma(3/2 + k + ν))
        if next_sum == prev_sum
            break
        end
        total_sum += next_sum
        k += 1
        prev_sum = next_sum
    end
    coefficient * total_sum
end

function I(ν, z)
    ν = BigInt(ν)
    z = BigFloat(z)
    total_sum = 0.0
    k = 0
    prev_sum = 0.0
    while true
        next_sum = (z / 2)^(2 * k + ν) / (gamma(k + 1) * gamma(k + ν + 1))
        if next_sum == prev_sum
            break
        end
        total_sum += next_sum
        k += 1
        prev_sum = next_sum
    end
    total_sum
end

function I_minus_L(n, x, a)
    setprecision(BigFloat, digits)
    n = BigFloat(n)
    x = BigFloat(x)
    a = BigFloat(a)

    d = weval(W"N"(W"Abs"(W"Subtract"(W"BesselI"(n+1, sign(x)^(n + 1) * x * a), W"StruveL"(-n-1, sign(x)^(n) * x * a))), digits)).value
    d_split = split(d, r"[`^]")
    if length(d_split) == 2
        return parse(BigFloat, d_split[1])
    else
        return parse(BigFloat, d_split[1] * "e" * d_split[3])
    end
end

function ℑχ_0(Ω, α, v, w, N = 10)

    R = (v^2 - w^2) / (w^2 * v)
    coefficient = 2 / 3 * α * (v / w)^3

    total_sum = 0.0
    for n in 0:Int(N)
        total_sum += - coefficient * sqrt(π) / (gamma(-n - 1/2) * gamma(n + 1)) * (-2 * R)^n / BigCombinatorics.doublefactorial(2 * n + 1) * abs(Ω - 1 - n * v)^(n + 1/2) * exp(-R * abs(Ω - 1 - n * v)) * (1 + sign(Ω - 1 - n * v))
    end
    return total_sum
end

function ℑχ(Ω, β, α, v, w, N = 10)

    Ω = BigFloat(Ω)
    β = BigFloat(β)
    α = BigFloat(α)
    v = BigFloat(v)
    w = BigFloat(w)

    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = 2 * α / (3 * sqrt(π)) * (v / w)^3 * β^(3/2) * sinh(β * Ω / 2) / sinh(β / 2)

    @show(Ω)

    total_sum = BigFloat(0.0)

    for n in 0:Int(N)

        n_coefficient = -2 * beta(-1/2 - n, n + 1)^(-1) * (-1)^n * b^n * (4 * a)^(-n - 1) * sqrt(π) / gamma(n + 3/2)

        if isodd(n)
            for k in -1:Int((n - 1) / 2)
                k_coefficient_inverse = ((n + 1) * beta(n - k + 1, k + 1))

                y = [(Ω + 1 + v * (n - 2 * k)), (Ω + 1 - v * (n - 2 * k)), (Ω - 1 + v * (n - 2 * k)), (Ω - 1 - v * (n - 2 * k))]

                for x in y
                    total_sum +=  coefficient * n_coefficient * besselk(Int(n + 1), abs(Float64(x)) * Float64(a)) * abs(x)^(n + 1) / k_coefficient_inverse
                    # total_sum += coefficient * n_coefficient * k_coefficient * besselk(ArbReal(n + 1), ArbReal(abs(x) * a)) * abs(x)^(n + 1)
                end
            end
        end

        if iseven(n)
            for k in -1:Int(n / 2 - 1)

                k_coefficient_inverse = ((n + 1) * beta(n - k + 1, k + 1))
                y = [(Ω + 1 + v * (n - 2 * k)), (Ω + 1 - v * (n - 2 * k)), (Ω - 1 + v * (n - 2 * k)), (Ω - 1 - v * (n - 2 * k))]

                for x in y
                    total_sum +=  coefficient * n_coefficient * besselk(Int(n + 1), abs(Float64(x)) * Float64(a)) * abs(x)^(n + 1) / k_coefficient_inverse
                    # total_sum += coefficient * n_coefficient * k_coefficient * besselk(ArbReal(n + 1), ArbReal(abs(x) * a)) * abs(x)^(n + 1)
                end

                k_even_coefficient_inverse = ((n + 1) * beta(n/2 + 1, n/2 + 1))
                y_even = [(Ω + 1), (Ω - 1)]

                for x in y_even
                    total_sum +=  coefficient * n_coefficient * besselk(Int(n + 1), abs(Float64(x)) * Float64(a)) * abs(x)^(n + 1) / k_even_coefficient_inverse
                    # total_sum += coefficient * n_coefficient * k_even_coefficient * besselk(ArbReal(n + 1), ArbReal(abs(x) * a)) * abs(x)^(n + 1)
                end
            end
        end
    end
    println("/n")
    return total_sum
end

function ℜχ(Ω, β, α, v, w, N = 10)
    Ω = BigFloat(Ω)
    β = BigFloat(β)
    α = BigFloat(α)
    v = BigFloat(v)
    w = BigFloat(w)

    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    @show(Ω)

    coefficient = 2 * α / (3 * sqrt(π)) * (v / w)^3 * β^(3/2) / sinh(β / 2)

    # total_sum = BigFloat(0.0)
    #
    # for n in 0:Int(N)
    #
    #     n_coefficient = sinh(β * Ω / 2) * (1 / (-1/2 * beta(n + 1, -n - 1/2))) * π *  b^n * (1/4) * (1/a^(n + 1)) * Factorial(n + 1) / Factorial(2 * (n + 1))
    #
    #     if isodd(n)
    #         for k in -1:Int((n - 1) / 2)
    #
    #             k_coefficient = Binomial(n, k) * (1 - (-1)^n)
    #
    #             y = [(Ω + 1 + v * (n - 2 * k)), (Ω + 1 - v * (n - 2 * k)), (Ω - 1 + v * (n - 2 * k)), (Ω - 1 - v * (n - 2 * k))]
    #
    #             for x in y
    #                 total_sum += coefficient * n_coefficient * k_coefficient * I_minus_L(n, x, a) * x * abs(x)^n
    #             end
    #         end
    #     end
    #
    #     if iseven(n)
    #         for k in -1:Int(n / 2 - 1)
    #
    #             k_coefficient = Binomial(n, k) * (1 - (-1)^n)
    #
    #             y = [(Ω + 1 + v * (n - 2 * k)), (Ω + 1 - v * (n - 2 * k)), (Ω - 1 + v * (n - 2 * k)), (Ω - 1 - v * (n - 2 * k))]
    #
    #             for x in y
    #                 total_sum += coefficient * n_coefficient * k_coefficient * I_minus_L(n, x, a) * x * abs(x)^n
    #             end
    #
    #             k_even_coefficient = -Binomial(n, k) * (1 - (-1)^n) / ((n + 1) * beta(n/2 + 1, n/2 + 1))
    #
    #             y_even = [(Ω + 1), (Ω - 1)]
    #
    #             for x in y_even
    #                 total_sum += coefficient * n_coefficient * k_coefficient * I_minus_L(n, x, a) * x * abs(x)^n
    #             end
    #         end
    #     end
    # end

    first_integrand(t) = sinh(Ω * β / 2) * sin(Ω * t) * cos(t) / (a^2 + t^2 - b*cos(v * t))^(3 / 2)

    first_term = coefficient * quadgk(t -> first_integrand(t), 0, Inf)[1]

    second_integrand(t) = (1 - cosh(Ω * (t - β / 2))) * cosh(t) / (a^2 - t^2 - b * cosh(v * t))^(3 / 2)

    second_term = coefficient * quadgk(t -> second_integrand(t), 0, β / 2)[1]

    return first_term + second_term
end

function realX(Ω, β, α, v, w, N=10)

    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    Q = 1 / (exp(β) - 1)

    ℜD(x) = w^2 / v^2 * (R * (1 - cos(v * x)) + 4 * R * P * (sin(v * x / 2))^2 + x^2 / β)

    ℑD(x) = -w^2 / v^2 * (R * sin(v * x) - x)

    RD(x) = sqrt(ℜD(x)^2 + ℑD(x)^2)
    θD(x) = atan(ℑD(x), ℜD(x))

    ℑS(x) = 2 * α / (3 * sqrt(π)) * ((1 + Q) * sin(x - 3 * θD(x) / 2) - Q * sin(x + 3 * θD(x) / 2)) / RD(x)^(3 / 2)

    integrand(x) = -(1 - cos(Ω * x)) * ℑS(x) / 150
    @show(Ω)
    return quadgk(x -> integrand(x), 0, Inf)[1]
end

function ℜχ_0(Ω, α, v, w, N = 10)
    R = (v^2 - w^2) / (w^2 * v)
    integrand(x, n) = ((n + 1/2) * x^(n - 1/2) * exp(-R * x) - R * x^(n + 1/2) * exp(-R * x)) * log(abs((1 + n * v + x)^2 / (Ω^2 - (1 + n * v + x)^2))^(1/2))
    total_sum = 0.0
    for n in 0:Int(N)
        total_sum += -(1 / gamma(n + 3/2)) * quadgk(x -> integrand(x, n), 0, Inf)[1]
    end
    return total_sum
end



# ImX(x, y, α, v, w) = ℑχ(x, y, α, v, w, 50)
#
# p = (:α, :v, :w)
# i = (5, 4.0, 2.1)
# pl = ((0.01, 10), (0.01, 10), (0.01, 10))
# al = ((0.01, 100), (0.0001, 100), (0.0, 200))
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
# vf = viewfunction(ImX, ImX, p, i, pl, a, axis_limits = al, len = 500)

β = 10
v = range(3.0, stop = 21.99, length = 100)
params = (5.0, 4.0, 2.1, 60)
# ImX = [ℑχ(t, β, params...) for t in v]
# ImX_0 = [ℑχ_0(t, params...) for t in v]
# plot(v, ImX_0, label = "β = ∞")
# plot!(v, ImX, label = "β = $(β)")
# xlabel!("Ω")
ReX_0 = [ℜχ_0(t, params...) for t in v]
# ReX = [ℜχ(t, β, params...) for t in v]
plot(v, ReX_0)
# plot!(v, ReX)
xlabel!("Ω")
ylabel!("Reχ(Ω)")
