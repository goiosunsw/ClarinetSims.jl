module ModalClarinet
	using OrdinaryDiffEq

	function resonatorDE end

	abstract type AcousticMode end

	struct ComplexAcousticMode
		pole::Complex
		residue::Complex
	end

	# Acoustic resonators can be implemented in multiple ways
	abstract type AcousticResonator end

	abstract type ModalResonator <: AcousticResonator end


	struct ComplexModesResonator <: ModalResonator
		modes::Vector{ComplexAcousticMode}
	end

	# Reeds, lips and vocal folds can also be implemented in different ways
	abstract type FlowValve end

	struct MasslessReed <: FlowValve
		stiffness::Real
	end

	# instruments may also have different implementations?

	abstract type AbstractReedInstrument end

	struct ReedInstrument{T} <: AbstractReedInstrument where {T <: AcousticResonator}
		resonator::T
		exciter::FlowValve
		couplingFunction::Function
	end

	struct SimulationParameters
		rate::Real
		duration::Real
	end

	struct ReedInstrumentSimulation
		instrument::AbstractReedInstrument
		simulationParameters::SimulationParameters
	end

	"""
		modalPoleRes(ω::Real, q::Real; a::Real=1)

	Return pole and residue for a second order system
	"""
	function modalPoleRes(ω::Real, q::Real; a::Real=1)
		inv2q = 1 / (2*q)
		inv2q2 = inv2q^2
		s = ω * (-inv2q + 1im* sqrt( Complex(1 - inv2q2)))
		c = a * ω / q / 2 * (1 - 1im/( 2q * sqrt( Complex(1 - inv2q2))))
		ComplexAcousticMode(s,c)
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
		xdiff = -x-x_st
		if xdiff > 0.0
			# parabolic
			return (sqrt(xdiff^2+1)-1)*.8
			# chatziioanou
			#return log((6.4*xdiff/(x_ev-x_st))+1)
		else
			return 0.0
		end
	end


	function oneModeClarinet(du,u,par,t)
		p_r = u[1]
		p_i = u[2]
		p = 2*p_r

		γ, ζ = par

		x=(p-γ)/(1+nlstiff(p-γ))

		a = √(abs(γ-p)+α^2)
		if x < -10.0
			v = 0
		else
			v = ζ * (1+x) * sign(γ - p) * (a-α) 
		end

		pc = Complex(p_r, p_i)
		dpc = resonatorDE(pc, v)
		
		du[1] = real(dpc)
		du[2] = imag(dpc)
	end

	function simToDiffEqs(sim)
		inst = sim.instrument
		res = inst.resonator
		setResonator(res)
		exc = inst.exciter
		par = (exc.γ, exc.ζ)
	end

	function setResonator(res::ComplexModesResonator)
		resonatorDE(p_complex::Complex, v::Real) = res.c*v + res.s*p_complex
	end


	function runSimulation(sim::ReedInstrumentSimulation)
		diffEqs, initialValues = simToDiffEqs(sim)
		dt = 1/sim.rate
		prob = ODEProblem(diffEqs, initialValues, sim.duration, dt=dt)
		sol = solve(prob, Rosenbrock23())
	end
end #module