mata
mata clear
/* Set a structure including two matrices */
struct two_matrices {
    real Matrix1
    real Matrix2
}


function compute_SE(IF, q){
// 	struct two_matrices scalar res
	
	M = rows(IF)
	Avar = q*IF'*IF/(M-1)
	SE = diagonal(sqrt(Avar/M))
	
// 	res.Matrix1 = Avar/M
// 	res.Matrix2 = SE
	
	return(Avar/M)
}

real matrix compute_IF(X,Xr){
    //struct two_matrices scalar res
	
	N = rows(X)			/* Number of observation */
	K = cols(X)			/* Number of regressros */
	M = rows(Xr)		/* Number of clusters */
	q = ((N-1)/(N-K)) 	/* Finite Sample Adjustment */		
	IF = (invsym(X'*X/M)*Xr')'
	
	return(IF)
}

/* Set mata functions to compute empirical influence and standard error */
// function compute_IF(X,Xr){
//     struct two_matrices scalar res
//	
// 	N = rows(X)			/* Number of observation */
// 	K = cols(X)			/* Number of regressros */
// 	M = rows(Xr)		/* Number of clusters */
// 	q = ((N-1)/(N-K)) 	/* Finite Sample Adjustment */	
//	
// 	IF = (invsym(X'*X/M)*Xr')'	/* Influence function of LS estimator */
// 	Avar = q*IF'*IF/(M-1)		/* Compute asymptotic variance matrix */
// 	SE = diagonal(sqrt(Avar/M))	/* Compute plug-in clustered standard error */
//	
// 	res.Matrix1 = IF
// 	res.Matrix2 = SE
//     return(res)
// }

/* Recover the decomposition */
real matrix est_pi(beta, zeta, xi, NN, IF_b, IF_z, IF_x,q){
	pi 		= compute_decomp(beta, zeta, xi, NN)
	IF 		= compute_IF_decomp(beta,zeta,xi,NN,IF_b,IF_z,IF_x)
	Vhat 	= compute_SE(IF,q)	 	
	return(pi, Vhat)
}

/* Compute influence functions of the decomposition */
real matrix compute_decomp(b,z,x,N){
	pi1 = b[2]
	pi2 = N*b[4]*x[2]
	pi3 = (b[3]-b[4])*z[1]
	pi4 = b[3]*x[3]
	pi = (pi1\pi2\pi3\pi4)
	return(pi)
}
real matrix compute_IF_decomp(b,z,x,N,IF_b,IF_z,IF_x){
	IF_pi1 = IF_b[,2]
	IF_pi2 = N*(b[4]*IF_x[,2] + IF_b[,4]*x[2])
	IF_pi3 = (IF_b[,3]-IF_b[,4])*z[1] + (b[3]-b[4])*IF_z[,1]
	IF_pi4 = IF_b[,3]*x[3] + b[3]*IF_x[,3]
	IF_pi = (IF_pi1, IF_pi2, IF_pi3, IF_pi4)
	return(IF_pi)
}

end



cap program drop reg_dyadic
program reg_dyadic, eclass
/* Program for dyadic regression */
/* Input
	- varlist: (link) (Di) (Dj) 
	- cov: additional covariates if any
	- idx: specify indices of group, i, j
	- cl_g: clustering within groups
	- cl_i: clustering within individuals */
qui{
syntax varlist (min = 3) [,cov(varlist) idx(varlist) cl_g cl_i]
tokenize `varlist'

/* Generate Dyadic Regressor W = (1, Di, Dj, Di x Dj) */
for any iota Di_tmp Dj_tmp Dij_tmp: cap drop X
gen iota		= 1
gen Di_tmp 		= `2'
gen Dj_tmp 		= `3'
gen Dij_tmp		= Di_tmp*Dj_tmp

local W iota Di_tmp Dj_tmp Dij_tmp

/* Store Index Variables (g,i,j) */
if "`idx'" != "" {
	local i 1
	foreach var of varlist `idx' {
		local idx_var`i' `var'
		local ++i
	}
}

/* Set clustering group and run dyadic regression */	
local cluster_g = "`cl_g'" != ""
local cluster_i = "`cl_i'" != ""

if `cluster_g' == 1{
	local cl `idx_var1'
	noi: reg `1' `W' `cov', nocon cluster( `cl' )
}
else if `cluster_i' == 1{
	local cl `idx_var2'
	noi: reg `1' `W' `cov', nocon cluster( `cl' )
}
else{
	local cl 
	noi: reg `1' `W' `cov', nocon r
}

/* Store residual and predicted links */
for any r_tmp `1'_hat: cap drop X
predict `1'_hat
predict r_tmp, resid

/* Compute Influence Functions */
mata{
	W 	= st_data(.,"`W'")
	r 	= st_data(.,"r_tmp")
	Wr 	= W:*r
}

/* Convert mata matrix to variables */
cap drop 	Wr_tmp*
getmata 	Wr_tmp_* = Wr

if "`cl'" != "" {
*	di "clustered"
	preserve
		collapse (sum) Wr_tmp_*, by( "`cl'" )
		mata: Wr_cl = st_data(.,("Wr_tmp_*"))
		mata: IF_`1' = compute_IF(W, Wr_cl)
	restore
}
else{	
*	di "not clustered"
	mata: IF_`1' = compute_IF(W, Wr)

}

/* Remove temporal variables */
for any iota Di_tmp Dj_tmp Dij_tmp r_tmp Wr_tmp*: cap drop X
}


end





cap program drop compute_QRS
program compute_QRS
syntax varlist (min = 3 max = 3) [,savenm(string)]
qui{
	tokenize `varlist'
	local A_hat `1'
	local Di  	`2'
	local Dj  	`3'
	
	for any iota Di_tmp Dj_tmp Dij_tmp Q_tmp R_tmp S_tmp A_hat_tmp DQ*_tmp DR*_tmp DS*_tmp: cap drop X
	gen iota		= 1
	gen Di_tmp 		= `2'
	gen Dj_tmp 		= `3'
	gen Dij_tmp		= Di_tmp*Dj_tmp
	gen A_hat_tmp  	= `1'
	gen Q_tmp 		= A_hat_tmp * Dj_tmp
	gen R_tmp 		= A_hat_tmp * (1-Dj_tmp)
	gen S_tmp 		= A_hat_tmp
	
	local W iota Di_tmp Dj_tmp Dij_tmp

	local vn = 0
	foreach v in `W'{
		local ++vn
		gen DQ`vn'_tmp = `v'*Dj_tmp
		gen DR`vn'_tmp = `v'*(1-Dj_tmp)
		gen DS`vn'_tmp = `v'
	}
	preserve
		collapse (sum) Q_tmp R_tmp S_tmp DQ*_tmp DR*_tmp DS*_tmp, by(g i)
		save `savenm', replace
	restore
	for any iota Di_tmp Dj_tmp Dij_tmp Q_tmp R_tmp S_tmp A_hat_tmp DQ*_tmp DR*_tmp DS*_tmp: cap drop X
}
end



cap program drop est_display
program est_display
qui{
	syntax namelist(min=2 max=2) [,names(string)]

	gettoken mat_b  namelist : namelist
	gettoken mat_V namelist : namelist

	matrix colnames `mat_b' = `names'
	matrix colnames `mat_V' = `names'
	matrix rownames `mat_V' = `names'
	ereturn post `mat_b' `mat_V'
	}
	noi ereturn display
end



cap program drop estim_RE
program estim_RE, eclass

syntax varlist (min = 4)
tokenize `varlist'
local Y  `1'
local D  `2'
local Q  `3'
local R  `4'
macro shift 4
local regressors `*'


cap gen iota = 1
local Z iota `D' `Q' `R'
qui reg `Y' `Z' `regressors', noc
mat beta_RE = e(b)
mat beta_RE = beta_RE[1,1..4]

cap drop 	r_tmp
predict 	r_tmp, resid

/* Compute Influence Functions */
mata{
	Z 		= st_data(.,"`Z'")
	r 		= st_data(.,"r_tmp")
	Zr 		= Z:*r
}
cap drop 	Zr_tmp* 
getmata 	Zr_tmp_* = Zr
preserve
	collapse (sum) Zr_tmp_*, by( g )
	mata: Zr_cl = st_data(.,("Zr_tmp_*"))
restore


mata{
N = rows(Z)
K = cols(Z)
q = ((N-1)/(N-K)) 
M = rows(Zr_cl)

/* Influence of beta without adjustment */
IF_Y1 		= compute_IF(Z, Zr_cl)

/* Standard error adjustment for beta*/
Dq 			= st_data(.,"DQ*_tmp")
Dr 			= st_data(.,"DR*_tmp")
beta_RE 	= st_matrix("beta_RE")'
zeta1 		= st_matrix("zeta1")

Diff_Zb  	= beta_RE[3]:*Dq + beta_RE[4]:*Dr
IF_adj  	= -(invsym(Z'*Z/M)*((Z'*Diff_Zb )/M)*IF_A1'/M)'
IF_beta_RE 	= IF_Y1 + IF_adj
Vhat_RE 	= compute_SE(IF_beta_RE,q)

/* Compute SE of the decomposition */
res_pi = est_pi(beta_RE, zeta1, zeta1, NN, IF_beta_RE,IF_A1,IF_A1,q)


st_matrix("V_beta_RE", Vhat_RE)
st_matrix("pi_RE", res_pi[,1]')
st_matrix("V_pi_RE", res_pi[,2..5]')
}

est_display beta_RE V_beta_RE, names(`Z')
est_display pi_RE V_pi_RE, names(Direct_Treatment Direct_Network Indirect_Treatment Indirect_Network)
end




cap program drop estim_PT
program estim_PT, eclass

syntax varlist (min = 4)
tokenize `varlist'
local DY  `1'
local D  `2'
local Q  `3'
local RS  `4'
macro shift 4
local regressors `*'

cap gen iota = 1
local X iota `D' `Q' `RS' 
qui reg `DY' `X' `regressors', noc r 
mat beta_PT = e(b)
mat beta_PT = beta_PT[1,1..4]

cap drop e_tmp
predict e_tmp, resid

/* Compute Influence Functions */
mata{
	X		= st_data(.,"`X'")
	e 		= st_data(.,"e_tmp")
	Xe 		= X:*e
}

cap drop Xe_tmp*
getmata Xe_tmp_* = Xe

preserve
	collapse (sum) Xe_tmp_*, by( g )
	mata: Xe_cl = st_data(.,("Xe_tmp_*"))
restore


mata{
N = rows(X)
K = cols(X)
q = ((N-1)/(N-K)) 
M = rows(Xe_cl)
IF_DY = compute_IF(X, Xe_cl)

/* Standard error adjustment */
Dq 		= st_data(.,"DQ*_tmp")
Dr 		= st_data(.,"DR*_tmp")
Ds 		= st_data(.,"DS*_tmp")
beta_PT = st_matrix("beta_PT")'
zeta1 	= st_matrix("zeta1")
zeta0 	= st_matrix("zeta0")
xi 		= st_matrix("xi") 

Diff_Zb1 = beta_PT[3]:*Dq + beta_PT[4]:*Dr
Diff_Zb0 = beta_PT[4]:*Ds

IF_adj1 = -(invsym(X'*X/M)*((X'*Diff_Zb1)/M)*IF_A1'/M)'
IF_adj0 = -(invsym(X'*X/M)*((X'*Diff_Zb0)/M)*IF_A0'/M)'

IF_beta_PT 	= IF_DY + IF_adj1 - IF_adj0
Vhat_PT 	= compute_SE(IF_beta_PT,q)/M

res_pi = est_pi(beta_PT, zeta1, xi, NN, IF_beta_PT,IF_A1,IF_DA,q)

st_matrix("V_beta_PT", Vhat_PT)
st_matrix("pi_PT", res_pi[,1]')
st_matrix("V_pi_PT", res_pi[,2..5]')
}

est_display beta_PT V_beta_PT, names(`X')
est_display pi_PT V_pi_PT, names(Direct_Treatment Direct_Network Indirect_Treatment Indirect_Network)

end

