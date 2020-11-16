# D(u)_function.jl

"""
Implementation of FHIP's susceptibility formula D(u) for the Optical Absorption of Polarons.

See FHIP 1962, equation (35c):
https://link.aps.org/doi/10.1103/PhysRev.127.1004.
"""

"""
ℜD(x::Float64, y::Float64, v::Float64, w::Float64, β::Float64)

    Calculate the real part of D(u) (equation (35c) in FHIP) for a complex argument u = x + iy. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. β is the reduced thermodynamical beta = ħω/kT with k being Boltzman's constant.
"""
function ℜD(x, y, v, w, β)
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    w^2 / v^2 * (
        R * cos(v * x) * (exp(-v * y) + 2 * P * cosh(v * y)) - (x^2 - y^2) / β - y -
        R * (2 * P + 1)
    )
end

"""
ℑD(x::Float64, y::Float64, v::Float64, w::Float64, β::Float64)

    Calculate the imaginary part of D(u) (equation (35c) in FHIP) for a complex argument u = x + iy. v and w are the variational Polaron parameters that minimise the free energy, for the supplied α Frohlich coupling. β is the reduced thermodynamical beta = ħω/kT with k being Boltzman's constant.
"""
function ℑD(x, y, v, w, β)
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    w^2 / v^2 * (R * sin(v * x) * (exp(-v * y) - 2 * P * sinh(v * y)) - x * (2 * y / β - 1))
end

# Create a complex number from the real and imaginary parts of D(u).

D(x, y, v, w, β) = ℜD(x, y, v, w, β) + 1im * ℑD(x, y, v, w, β)
