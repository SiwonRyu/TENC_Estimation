clear all

cd "`c(pwd)'"
qui do dofiles/def_estim.do

qui{
	/* Dyadic Data */
	use dyadic_dataset_ComolaPrina.dta,clear
	
	gen g = village_code

	/* Set Indices */
	gen 	i 		= HHNO1
	gen 	j 		= HHNO2

	/* Set Regressors */
	gen 	iota 	= 1
	gen 	Di 		= ITT_H11
	gen 	Dj 		= ITT_H21
	gen 	Dij 	= Di*Dj
	gen 	DipDj 	= Di+Dj
	egen 	DimDj 	= rowmax(Di Dj) 

	/* Drop self link */
	drop if i == j

	/* Set Links */
	foreach t in 0 1{
		gen A`t' 	= link_bin`t'
		gen A`t'_rn	= A`t'
		bys i: egen deg`t' = total(A`t')
		replace A`t'_rn = A`t'/deg`t' if deg`t' > 0
	}
	gen DA = A1-A0
	order g i j Di Dj A0 A1 DA

	foreach var in abs_diff_HHmembers abs_diff_chl16 both_married ///
		some_shock_livestock some_shock_death {
		gen D_`var' = `var'1-`var'0
	}
}

reg_dyadic A0 Di Dj, cov(${controls}) idx(g i j) cl_g
mat zeta0 = e(b)'
reg_dyadic A1 Di Dj, idx(g i j) cl_g
mat zeta1 = e(b)'
reg_dyadic DA Di Dj, idx(g i j) cl_g
mat xi = e(b)'

preserve
	collapse (sum) A0 A1 (mean) deg0 deg1 ,by (g i)
	noi: di "Degree Distribution"
	noi: su A0 A1 deg0 deg1
	gen nn = 1
	collapse (sum) nn , by(g)
	noi: di "Average N"
	noi: su nn
restore

compute_QRS A1_hat Di Dj, savenm(Z1_tmp)
compute_QRS A0_hat Di Dj, savenm(Z0_tmp)


/* Outcome Regression */
qui{
	mata: NN = 20 /* counterfactual number of neighbors */
	use individual_dataset_ComolaPrina.dta, clear

	/* Set Indices */
	gen g = village_code
	tab g, generate(gdum)
	gen i = HHNO
	gen D = ITT
	gen iota = 1
	bys g: egen Ng = count(i)

	merge 1:1 i using Z1_tmp
	drop if _merge != 3
	drop _merge

	gen YD = y1 > 0
	gen Y1 = y1
	gen lY1e = log(y1)
	gen lY1 = log(y1+1)
	gen lY0 = log(y0+1)
	ren Q_tmp Q
	ren R_tmp R
}
estim_RE Y1 D Q R gdum*
estim_RE YD D Q R gdum*
estim_RE lY1e D Q R gdum*

qui{
	drop S*
	merge 1:1 i using Z0_tmp
	drop if _merge != 3
	drop _merge

	gen DY 	= y1-y0
	gen DlY = log(y1+1)-log(y0+1)
	ren S_tmp S
	gen RS = R-S
}
estim_PT DY D Q RS gdum*
estim_PT DlY D Q RS gdum*

cap erase Z0_tmp.dta
cap erase Z1_tmp.dta




