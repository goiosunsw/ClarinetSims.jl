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

    function SingleComplexAcMode(ω::Real, Q::Real, A::Real)

        inv2q = 1 / (2*Q)
        inv2q2 = inv2q^2
        s = ω * (-inv2q + 1im* sqrt( Complex(1 - inv2q2)))
        c = A * ω / Q / 2 * (1 - 1im/( 2Q * sqrt( Complex(1 - inv2q2))))
        SingleComplexAcMode( s, c )		
    end

    function acousticModeDE(du, u, par, t; ω₀=1.0, Q=1.0)
        ω₀, Q, A = par
        du[1] = u[2]
        du[2] = - ω₀^2 * u[1] - ω₀/Q * u[2]
    end

    function generateAcousticProblem(sam::SingleAcousticMode)
        ω₀ = sam.ω
        Q = sam.Q
        A = sam.A
        f!(du,u,par,t) = let 
            dv = par[2](t)
            du[1] = u[2]
            du[2] = - ω₀^2 * u[1] - ω₀/Q * u[2] + A*dv
        end 
        return f!
    end

    function generateAcousticProblem(sam::SingleComplexAcMode)
        s = sam.s
        c = sam.c
        f!(du,u,par,t) = let 
            v = par[1](t)
            p_r = u[1]
            p_i = u[2]
            p_c = Complex(p_r, p_i)
            dp = c*v + s*p_c
            du[1] = real(dp)
            du[2] = imag(dp)
        end 
        return f!
    end
    
    function ImpResp(sam::T) where T <: AcousticMode
        f! = generateAcousticProblem(sam)
        u₀ = [1.0,0.0]
        tspan = [0.0,1.0]
        pars = [x -> (x < 0.001 ? 1.0 : 0.0),x -> (x < 0.001 ? 1.0 : 0.0)]
        
        prob = ODEProblem(f!, u₀, tspan, pars, dt=0.0001)
        sol = solve(prob,Tsit5())#Rosenbrock23())
        return sol
    end
end