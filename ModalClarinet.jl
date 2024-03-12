using OrdinaryDiffEq

struct ComplexAcousticMode
	pole::Complex
	residue::Complex
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



function clarinet(du,u,par,t)
    pr = u[1]
	pi = u[2]
	p = 2*pr

    x=(p-γ)/(1+nlstiff(p-γ))

    a = √(abs(γ-p)+α^2)
    if x < -10.0
        v = 0
    else
        v = ζ * (1+x) * sign(γ - p) * (a-α) 
    end

	pc = Complex(pr,pi)
    dpc = c*v + s*pc
    
    du[1] = real(dpc)
	du[2] = imag(dpc)

end
