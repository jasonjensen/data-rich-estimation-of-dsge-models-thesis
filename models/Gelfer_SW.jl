# Linearized Smets & Wouters 2007 model as written in Gelfer (2019)

@model SW begin
	K[0] = (1-τ) * K[-1] + τ * I[0] + τ * (1 + β)*Φ * εI[0]

	L[0] = -w[0] + (1 + 1/Ψ)*rk[0] + K[-1]

	Y[0] = Cy_bar * C[0] + Iy * I[0] + rk_ss*Ky_bar/Ψ * rk[0] + εG[0]

	Y[0] = φ*εa[0] + φ*α*K[-1] + φ*α/Ψ*rk[0] + φ*(1-α)*L[0]

	R[0] = ρ*R[-1] + (1-ρ)*(rπ1*π[0] + ry1*Y[0] + rπ2*π[-1] + ry2*Y[-1]) + εr[0]

	C[0] = h/(1+h)*C[-1] + 1/(1+h)*C[1] - (1-h)/((1+h)*σC)*(R[0] - π[1]) + εb[0]

	I[0] = 1/(1 + β)*I[-1] + β/(1+β)*I[1] + 1/((1+β)*Φ)*q[0] + εI[0]

	q[0] = -1*(R[0] - π[1]) +(1-τ)/(1-τ+rk_ss)*q[1] + rk_ss/(1-τ+rk_ss)*rk[1] + εq[0]

	π[0] = β/(1+β*ιp)*π[1] + ιp/(1 + β*ιp)*π[-1] + (1-β*ξp)*(1-ξp)/((1+β*ιp)*ξp)*(α*rk[0] + (1-α)*w[0] - εa[0]) + εp[0]

	w[0] = β/(1+β)*w[1] + 1/(1+β)*w[-1] + β/(1+β)*π[+1] - (1+β*ιw)/(1+β)*π[0] + ιw/(1+β)*π[-1] - (1-β*ξw)*(1-ξw)/((1+β)*(1+νl*(1+λw)/λw)*ξw)*(w[0] - νl*L[0] - σC/(1-h)*(C[0]-h*C[-1])) + εw[0]

	# shocks
	εI[0] = ρI*εI[-1] + σI*ϵI[x]
	εG[0] = ρG*εG[-1] + σG*ϵG[x]
	εb[0] = ρb*εb[-1] + σb*ϵb[x]
	εa[0] = ρa*εa[-1] + σa*ϵa[x]
	εw[0] = ρw*εw[-1] + σw*ϵw[x]
	εp[0] = ρp*εp[-1] + σp*ϵp[x]
	εq[0] = σq*ϵq[x]
	εr[0] = σr*ϵr[x]
end


@parameters SW begin
	Ψ  = 0.345
	ιp = 0.261
	ιw = 0.223
	ξp = 0.838
	ξw = 0.853
	νl = 2.009
	σC = 1.678
	h  = 0.688
	φ  = 0.445
	Φ  = 5.348

	rπ1 = 2.161
	ry1 = 0.345
	rπ2 = -0.222
	ry2 = -0.084
	ρ   = 0.867

	ρa = 0.911
	ρb = 0.772
	ρG = 0.974
	ρI = 0.710
	ρp = 0.827
	ρw = 0.524

	σa = 0.500
	σb = 0.085
	σG = 0.322
	σr = 0.125
	σI = 0.737
	σq = 0.104
	σp = 0.061
	σw = 0.048
	
	
	β = 0.99
	α = 0.3
	τ = 0.025
	Iy = 0.18
	Gy = 0.18
	λw = 0.3

	# derived parameters
	rk_ss = 1/β - (1-τ)
	Ky_bar = Iy * 1/τ 
	Cy_bar = 1 - Gy - Iy								# steady state federal funds rate ($\bar r$)
	
	# steadystates
	rk[ss] -> 0
end

# make sure it solves
get_solution(SW)

