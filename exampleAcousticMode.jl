include("AcousticModes.jl")
using .AcousticModes
using CairoMakie

ω = 10.0
q = 5.0
a = 1.0
dt = 0.001
tmax = 10.0

fig = Figure()
ax = Axis(fig[1,1])

# Real modes
sam = AcousticModes.SingleAcousticMode(ω, q, a)
sol = AcousticModes.impResp(sam, dt=dt, tmax=tmax)
vals = convert(Array, sol)
lines!(ax, sol.t, vals[1,1,:])

# Complex modes
sam = AcousticModes.SingleComplexAcMode(ω, q, a)
sol = AcousticModes.impResp(sam, dt=dt, tmax=tmax)
vals = convert(Array, sol)
lines!(ax, sol.t, vals[1,1,:])
#xlims!(0,0.5)

display(fig)