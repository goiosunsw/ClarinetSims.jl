
using OrdinaryDiffEq

"""
	modalPoleRes(ω::Real, q::Real; a::Real=1)

Return pole and residue for a second order system
"""
function modalPoleRes(ω::Real, q::Real; a::Real=1)
        inv2q = 1 / (2 * q)
        inv2q2 = inv2q^2
        s = ω * (-inv2q + 1im * sqrt(Complex(1 - inv2q2)))
        c = a * ω / q / 2 * (1 - 1im / (2q * sqrt(Complex(1 - inv2q2))))
        s, c
end


"""
	nlstiff(x::Real; x_st::Real=0.6, x_ev::Real=1.2)

Non-linear additional stiffness function from Chatziioanou
This stiffness should be added to a linear term as such:

``F = k_{lin} x + k_{nl}(x)``

where ``k_nl`` is this function

Returns the stiffness at a position `x`

Parameters:
`x_st`: Starting position for non-linear stiffness
`x_ev`: non-linearity coefficient
"""
function nlstiff(x::Real; x_st::Real=0.6, x_ev::Real=1.2)
        xdiff = -x - x_st
        if xdiff > 0.0
                # parabolic
                return (sqrt(xdiff^2 + 1) - 1) * 0.8
                # chatziioanou
                #return log((6.4*xdiff/(x_ev-x_st))+1)
        else
                return 0.0
        end
end


function oneModeClarinet(du, u, par, t)
        p_r = u[1]
        p_i = u[2]
        p = 2 * p_r

        # γ, ζ = par

        x = (p - γ) / (1 + nlstiff(p - γ))

        a = √(abs(γ - p) + α^2)
        if x < -10.0
                v = 0.0
        else
                v = ζ * (1.0 + x) * sign(γ - p) * (a - α)
        end

        pc = Complex(p_r, p_i)
        dpc = s * pc + c * v

        du[1] = real(dpc)
        du[2] = imag(dpc)
end


γ = 0.4
ζ = 0.3

α = 0.001
# resonator

ω₁ = 3141
q₁ = 30
f₁ = ω₁ * q₁
# reed
fᵣ = 1
ωᵣ = 12000
qᵣ = 2


s, c = modalPoleRes(ω₁, q₁, a=q₁)

# initial conditions (complex pressure)
u₀ = [0.001, 0]

sr = 48000.0
duration = 0.1
initialValues = [0.0, 0.0]
tspan = (0, duration)

# diffEqs, initialValues = simToDiffEqs(sim)
dt = 1.0 / sr
prob = ODEProblem(oneModeClarinet, u₀, tspan, dt=dt)
sol = solve(prob, Rosenbrock23())

