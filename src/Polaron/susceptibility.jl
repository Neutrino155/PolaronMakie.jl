using Combinatorics
using SpecialFunctions
using Plots
using QuadGK
using Quadmath
# using PolaronMobility
# MAPIe=polaronmobility(300, 4.5, 24.1, 2.25E12, 0.12)

function ℑχ(Ω, β, α, v, w, N = 10)

    R = (v^2 - w^2) / (w^2 * v)
    a = sqrt(β^2 / 4 + R * β * coth(β * v / 2))
    b = R * β / sinh(β * v / 2)

    coefficient = Float128(2 * α / (3 * sqrt(π)) * (v / w)^3 * β^(3/2) * sinh(β * Ω / 2) / sinh(β / 2))

    if coefficient == Inf

        coefficient = 2 / 3 * α * (v / w)^3

        total_sum = 0.0
        for n in 0:N
            total_sum += - coefficient * sqrt(π) / (gamma(-n - 1/2) * gamma(n + 1)) * (-2 * R)^n / doublefactorial(2 * n + 1) * abs(Ω - 1 - n * v)^(n + 1/2) * exp(-R * abs(Ω - 1 - n * v)) * (1 + sign(Ω - 1 - n * v))
        end
        return total_sum

    else

        total_sum = Float128(0.0)
        for n in 0:N

            n_coefficient = Float128(-2 * beta(-1/2 - n, n + 1)^(-1) * (-1)^n * b^n * (4 * a)^(-n - 1) * sqrt(π) / gamma(n + 3/2))

            if isodd(n)
                for k in 0:Int((n - 1) / 2)

                    k_coefficient = Float128(((n + 1) * beta(n - k + 1, k + 1))^(-1))

                    total_sum += Float128(coefficient * n_coefficient * k_coefficient * (
                    besselk(n + 1, abs(Ω + 1 + v * (n - 2 * k)) * a) * abs(Ω + 1 + v * (n - 2 * k))^(n + 1) +
                    besselk(n + 1, abs(Ω - 1 + v * (n - 2 * k)) * a) * abs(Ω - 1 + v * (n - 2 * k))^(n + 1) +
                    besselk(n + 1, abs(Ω + 1 - v * (n - 2 * k)) * a) * abs(Ω + 1 - v * (n - 2 * k))^(n + 1) +
                    besselk(n + 1, abs(Ω - 1 - v * (n - 2 * k)) * a) * abs(Ω - 1 - v * (n - 2 * k))^(n + 1)
                    ))
                end
            end

            if iseven(n)
                for k in 0:Int(n / 2 - 1)

                    k_coefficient = Float128(((n + 1) * beta(n - k + 1, k + 1))^(-1))

                    total_sum += Float128(coefficient * n_coefficient * k_coefficient * (
                    besselk(n + 1, abs(Ω + 1) * a) * abs(Ω + 1)^(n + 1) +
                    besselk(n + 1, abs(Ω - 1) * a) * abs(Ω - 1)^(n + 1)
                    ))

                    total_sum += Float128(coefficient * n_coefficient * k_coefficient * (
                    besselk(n + 1, abs(Ω + 1 + v * (n - 2 * k)) * a) * abs(Ω + 1 + v * (n - 2 * k))^(n + 1) +
                    besselk(n + 1, abs(Ω - 1 + v * (n - 2 * k)) * a) * abs(Ω - 1 + v * (n - 2 * k))^(n + 1) +
                    besselk(n + 1, abs(Ω + 1 - v * (n - 2 * k)) * a) * abs(Ω + 1 - v * (n - 2 * k))^(n + 1) +
                    besselk(n + 1, abs(Ω - 1 - v * (n - 2 * k)) * a) * abs(Ω - 1 - v * (n - 2 * k))^(n + 1)
                    ))
                end
            end
        end
        return total_sum
    end
end

function ℜχ(Ω, α, v, w, N = 10)

    R = (v^2 - w^2) / (w^2 * v)
    coefficient = 2 / 3 * α * (v / w)^3 * sqrt(π)

    ℜχ_integrand(Ω, v, n, x) = (n + 1/2) * x^(n - 1/2) * exp(- R * x) -  R * x^(n + 1/2) * exp(- R * x) * (1/2) * log(abs((1 + n * v + x)^3 / (Ω^2 - (1 + n * v + x)^2)))

    total_sum = 0.0

    for n in 0:N
        total_sum += (gamma(n + 3/2))^(-1) * QuadGK.quadgk(x -> ℜχ_integrand(Ω, v, n, x), 0, Inf)[1]
    end
    return coefficient/100 * total_sum
end

v = range(0, stop = 800, length = 1500)
params = (7, 2.8, 1.6, 2)
X = [ℑχ(t, 3, 7, 16.4, 5.6, 50) for t in v]
# Re = [ℜχ(t, params...) for t in v]
# Im = [ℑχ(t, params...) for t in v]
# G = [Γ_z(t, params...) for t in v]
plot(v, X)
# plot(v, Re)
# plot(v, Im, legend=false)
xlabel!("Ω")
ylabel!("Imχ(Ω)")
