# Linearized Smets & Wouters 2007 model with financial frictions as written in Gelfer (2019)
# Uses reciprocal function approximations to the COEF table columns from the Gelfer (2019) replication codes

# add some functions for δ coefficients
@eval MacroModelling begin
	δ_model(x,p) = p[1] .+ p[2] ./ (x .+ p[3])
	get_δRk(x) = δ_model(x, [0.9990447162600182, 0.009215599922294652, -0.020403301033697047])
	get_δR(x) = δ_model(x, [-0.004330298425182292, 0.009189115808462114, -0.020403200451239976])
	get_δqK(x) = δ_model(x, [0.006749813715453756, 5.2977534638673e-5, -0.020437235231941287])
	get_δn(x) = δ_model(x, [1.0033749068690714, 2.6488767192755553e-5, -0.020437235269791302])
	get_δσ(x) = δ_model(x, [0.0004199898881014732, 3.704968154136453e-5, -0.019684072435417508]) / δ_model(x, [0.027877821655470147, -0.0005156736022788782, 0.03192774545311349])
end
@model SWFF begin
	K[0] = (1-τ) * K[-1] + τ * I[0] + τ * (1 + β)*Φ * εI[0]

	L[0] = -w[0] + (1 + 1/Ψ)*rk[0] + K[-1]

	Y[0] = Cy_bar * C[0] + Iy * I[0] + rk_ss*Ky_bar/Ψ * rk[0] + εG[0]

	Y[0] = φ*εa[0] + φ*α*K[-1] + φ*α/Ψ*rk[0] + φ*(1-α)*L[0]

	R[0] = ρ*R[-1] + (1-ρ)*(rπ1*π[0] + ry1*Y[0] + rπ2*π[-1] + ry2*Y[-1]) + εr[0]

	C[0] = h/(1+h)*C[-1] + 1/(1+h)*C[1] - (1-h)/((1+h)*σC)*(R[0] - π[1]) + εb[0]

	I[0] = 1/(1 + β)*I[-1] + β/(1+β)*I[1] + 1/((1+β)*Φ)*q[0] + εI[0]

	Rk[0] - π[0] = (1-τ)/(1-τ+rk_ss)*q[0] + rk_ss/(1 - τ +rk_ss)*rk[0] - q[-1]

	π[0] = β/(1+β*ιp)*π[1] + ιp/(1 + β*ιp)*π[-1] + (1-β*ξp)*(1-ξp)/((1+β*ιp)*ξp)*(α*rk[0] + (1-α)*w[0] - εa[0]) + εp[0]

	w[0] = β/(1+β)*w[1] + 1/(1+β)*w[-1] + β/(1+β)*π[+1] - (1+β*ιw)/(1+β)*π[0] + ιw/(1+β)*π[-1] - (1-β*ξw)*(1-ξw)/((1+β)*(1+νl*(1+λw)/λw)*ξw)*(w[0] - νl*L[0] - σC/(1-h)*(C[0]-h*C[-1])) + εw[0]

	S[0] = χ*(q[0] + K[0] - n[0]) + εF[0] 

	Rk[1] = S[0] + R[0]
	
	n[0] = δRk*(Rk[0] - π[0]) - δR*(R[-1] - π[0]) + δqK*(q[-1] + K[-1]) + δn*n[-1] - δσ*εF[-1]
	
	# shocks
	εI[0] = ρI*εI[-1] + σI*ϵI[x]
	εG[0] = ρG*εG[-1] + σG*ϵG[x]
	εb[0] = ρb*εb[-1] + σb*ϵb[x]
	εa[0] = ρa*εa[-1] + σa*ϵa[x]
	εw[0] = ρw*εw[-1] + σw*ϵw[x]
	εp[0] = ρp*εp[-1] + σp*ϵp[x]
	εF[0] = ρF*εF[-1] + σF*ϵF[x]
	εr[0] = σr*ϵr[x]

end


@parameters SWFF begin
	Ψ  = 0.491
	ιp = 0.261
	ιw = 0.250
	ξp = 0.837
	ξw = 0.833
	νl = 1.782
	σC = 1.6624
	h  = 0.672
	φ  = 0.467
	Φ  = 2.716
	χ  = 0.051

	rπ1 = 2.196
	ry1 = 0.336
	rπ2 = -0.216
	ry2 = -0.103
	ρ   = 0.853

	ρa = 0.910
	ρb = 0.755
	ρG = 0.971
	ρI = 0.664
	ρF = 0.964
	ρp = 0.826
	ρw = 0.600

	σa = 0.487
	σb = 0.094
	σG = 0.327
	σr = 0.127
	σI = 0.955
	σF = 0.063
	σp = 0.061
	σw = 0.045

	β = 0.99
	α = 0.3
	τ = 0.025
	Iy = 0.18
	Gy = 0.19
	λw = 0.3
	γ = 0.99
	F = 0.0075
	S_ss = 1.00343 # from Gelfer code
	
	# derived parameters
	δRk = get_δRk(χ)
	δR = get_δR(χ)
	δqK = get_δqK(χ)
	δn = get_δn(χ)
	δσ = get_δσ(χ)
	
	rk_ss = S_ss * 1/β - (1-τ)
	Ky_bar = Iy * 1/τ
	Cy_bar = 1 - Gy - Iy
	
	# steadystates
	rk[ss] -> 0
	Y[ss] -> 0
	L[ss] -> 0
	K[ss] -> 0
	C[ss] -> 0
	I[ss] -> 0
	π[ss] -> 0
	R[ss] -> 0
	w[ss] -> 0
end

# check that it solves
get_solution(SWFF)
