module PolaronMakie

export viewfunction, frohlich_α, feynman_variation, osaka_free_energy, ℜD, ℑD, D, ℜS, ℑS, S, ℜχ, ℑχ, χ

using QuadGK
using ProgressBars
using AbstractPlotting.MakieLayout
using AbstractPlotting
# using GLMakie
using Optim

include("../src/function.jl")
include("../src/Functions/D(u)_function.jl")
include("../src/Functions/S(u)_function.jl")
include("../src/Functions/X(u)_function.jl")
include("../src/Polaron/frohlich.jl")
include("../src/Polaron/feynmantheory.jl")
include("../src/Polaron/osaka.jl")
include("../src/Polaron/hellwarththeory.jl")

# Physical constants
const ħ = 1.05457162825e-34; # kg m^2 s^{-1
const eV = 1.05457162825e-34; # kg m^2 s^{-2}
const m_e = 9.10938188e-31; # kg
const k_B = 1.3806504e-23; # kg m^2 K^{-1} s^2
const ϵ_0 = 8.854e-12 # C^2 N^{-1} m^{-2}

struct Polaron
    α::Float64      # Frohlich alpha
    v::Float64      # Variational parameter
    w::Float64      # Variational parameter
    F::Float64      # Free energy
    E::Float64      # Ground-state energy
    m_b::Float64    # Effective mass
    μ::Float64      # Mobility
end

function Polaron(α, v, w, F, E, m_b, μ)
    new(α, v, w, F, E, m_b, μ)
end

function Base.show(io::IO, quantities::Polaron)
    print(io, "$α, $v, $w, $F, $E, $m_b, $μ")
end

end # PolaronMakie
