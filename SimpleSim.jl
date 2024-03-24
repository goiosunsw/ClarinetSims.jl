
using OrdinaryDiffEq, Plots, FFTW

function modalPoleRes(ω, q; a=1)
    inv2q = 1 / (2 * q)
    inv2q2 = inv2q^2
    s = ω * (-inv2q + 1im * sqrt(Complex(1 - inv2q2)))
    c = a * ω / q / 2 * (1 - 1im / (2q * sqrt(Complex(1 - inv2q2))))
    s, c
end

function nlstiff(x; x_st=0.6, x_ev=1.2)
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

γ = 0.4
ζ = 0.3

α = 0.001

ω₁ = 3141
q₁ = 30
f₁ = ω₁ * q₁
# reed
fᵣ = 1
ωᵣ = 12000
qᵣ = 2

# initial conditions (complex pressure)
u₀ = [0.001, 0]
tspan = (0, 1.5)

s, c = modalPoleRes(ω₁, q₁, a=q₁)

function resonator(du, u, par, t)
    du[1] = c * par[1](t) + s * u[1]
end

function resonator2(du, u, par, t)
    dv = u[3]
    du[1] = u[2]
    du[2] = f₁ * dv - ω₁ / q₁ * u[2] - ω₁^2 * u[1]
end

prob_res = ODEProblem(resonator, [complex(1.0, 0)], tspan, [x -> (x < 0.001 ? 1 : 0)], dt=0.0001)
sol_res = solve(prob_res, Tsit5())#,Rosenbrock23())
#soln = solve(prob,Tsit5())

plot(sol_res.t, 2 .* [real(x[1]) for x in sol_res.u], linewidth=2, title="Resonator response", xaxis="Time", yaxis="Solution", idxs=[1], label=["p"])


function clarinet(du, u, par, t)

    pr = u[1]
    pi = u[2]
    p = 2 * pr

    x = (p - γ) / (1 + nlstiff(p - γ))

    a = √(abs(γ - p) + α^2)
    if x < -10.0
        v = 0.0
    else
        v = ζ * (1 + x) * sign(γ - p) * (a - α)
    end

    pc = Complex(pr, pi)
    #   dpc = c*v + s*pc
    dpc = [Complex(0.0, 0.0)]
    resonator(dpc, [pc], [x -> v], t)

    du[1] = real(dpc[1])
    du[2] = imag(dpc[1])

end

prob = ODEProblem(clarinet, u₀, tspan, dt=0.0001)
sol = solve(prob, Rosenbrock23(autodiff=false))
#soln = solve(prob,Tsit5())

function clarinet2(du, u, par, t)

    p = u[1]
    dp = u[2]
    dv = u[3]

    x = (p - γ) / (1 + nlstiff(p - γ))

    a = √(abs(γ - p) + α^2)
    if x < -10.0
        v = 0.0
    else
        v = ζ * (1 + x) * sign(γ - p) * (a - α)
    end

    resonator(du, u, [], t)


end

prob = ODEProblem(clarinet, u₀, tspan, dt=0.0001)
sol = solve(prob, Rosenbrock23(autodiff=false))
#soln = solve(prob,Tsit5())
dtf = 0.0001
nf = 1024
tf = 0:dtf:dtf*(nf-1)
ff = 1/dtf/nf:1/dtf/nf:1/dtf
fval = abs.(fft(2 * real(sol_res.(tf, idxs=1))))
print(ff[argmax(fval)], " ", ω₁ / 2pi)
plot(ff, fval)


# ╔═╡ 8592389a-524d-4677-ac90-c28e17ff7894
plot(sol, linewidth=2, title="Single-mode clarinet", xaxis="Time", yaxis="Solution", idxs=[1], label=["p"], xlim=(0, 0.06))


N = length(sol.u)
J = length(sol.u[1])
u = hcat(sol.u...)
p = 2 * u[1, :]
xx = (p .- γ) ./ (1 .+ nlstiff.(p .- γ))
flow = ζ * (1 .+ xx) .* sqrt.(abs.(γ .- p)) .* sign.(γ .- p)
plot(p, flow)
vline!([γ])

# ╔═╡ dbb459c9-eeaa-4462-9fc0-355850103b86
plot(p, xx)

# ╔═╡ 716eb020-a1b7-481e-88c4-29be3f310b5c

begin
    dt = 0.0001
    psol = hcat((sol.((maximum(sol.t)/3):dt:maximum(sol.t)))...)[1, :]
    plot(psol)
end
