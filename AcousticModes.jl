"""
    Modal acoustic resonators

    Acoustic reosnators defined by their resonant modes
"""
module AcousticModes

using OrdinaryDiffEq


abstract type AcousticMode end

struct SingleAcousticMode <: AcousticMode
    ω::Real
    Q::Real
    A::Real
end

struct SingleComplexAcMode <: AcousticMode
    s::Complex
    c::Complex
end

abstract type AbstractAcousticResonator end

struct ModalAcousticResonator{T} <: AbstractAcousticResonator where {T<:AcousticMode}
    modes::Vector{T}
end

struct DelayLineStraightTube <: AbstractAcousticResonator
    length::Real
    radius::Real
    reflectionCoeff::Real
end


function SingleComplexAcMode(ω::Real, Q::Real, A::Real)

    inv2q = 1 / (2 * Q)
    inv2q2 = inv2q^2
    s = ω * (-inv2q + 1im * sqrt(Complex(1 - inv2q2)))
    c = A * ω / Q / 2 * (1 - 1im / (2Q * sqrt(Complex(1 - inv2q2))))
    SingleComplexAcMode(s, c)
end

function acousticModeDE(du, u, par, t; ω₀=1.0, Q=1.0)
    ω₀, Q, A = par
    du[1] = u[2]
    du[2] = -ω₀^2 * u[1] - ω₀ / Q * u[2]
end

function generateAcousticProblem(sam::SingleAcousticMode)
    ω₀ = sam.ω
    Q = sam.Q
    A = sam.A
    f!(du, u, par, t) =
        let
            dv = par[2](t)
            du[1] = u[2]
            du[2] = -ω₀^2 * u[1] - ω₀ / Q * u[2] + A * dv
        end
    return f!, [1.0 0.0]
end

function generateAcousticProblem(sam::SingleComplexAcMode)
    s = sam.s
    c = sam.c
    f!(du, u, par, t) =
        let
            v = par[1](t)
            p_r = u[1]
            p_i = u[2]
            p_c = Complex(p_r, p_i)
            dp = c * v + s * p_c
            du[1] = real(dp)
            du[2] = imag(dp)
        end
    transformer = [1.0 0.0]
    return f!, transformer
end

function impulse_function(sr=1.0)
    impulse_t = 4.0 / sr
    x -> (x < impulse_t ? 1.0 / impulse_t : 0.0)
end

function impResp(sam::T; dt=1.0, tmax=1.0) where {T<:AcousticMode}
    f!, transformer = generateAcousticProblem(sam)
    outer_u₀ = [0.0]
    u₀ = outer_u₀ * transformer
    tspan = [0.0, tmax]
    impulse_t = dt * 4.0
    impulse_fun = x -> (x < impulse_t ? 1.0 / impulse_t : 0.0)
    impulse_dfun(x) =
        let
            if x < impulse_t / 2
                1.0 / impulse_t / dt / 2
            elseif x < impulse_t
                -1.0 / impulse_t / dt / 2
            else
                0.0
            end
        end
    pars = [impulse_fun, impulse_dfun]
    prob = ODEProblem(f!, u₀, tspan, pars, dt=dt)
    sol = solve(prob, Tsit5(), saveat=0:dt:tmax)#Rosenbrock23())
    return sol
end
end
