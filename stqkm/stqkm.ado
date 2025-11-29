*! v.1.0.0 - 23sept2014
*! N.Orsini, A.Discacciati 

capture program drop stqkm
program stqkm, eclass
version 10	

local cmdline : copy local 0

syntax [varlist(default=none)] [if] [in]  [ ,  ///
   Reps(integer -1)   ///
   MINReps(integer -1) /// undocumented
   Level(integer $S_level) ///
   SEED(string) ///
   Quantiles(numlist) ///
   WGTv(varlist min=1 max=1) /// undocumented
   NODOTS ///
   ]
	
	marksample touse, strok novarlist

	if "`log'"!="" | "`dots'"!="" {
			local log "*"
	}
	 
	// Check if any covariate. If none 
	
		if "`varlist'" == "" local overall 1
		else local overall 0
	
	// Check data st
		
		st_is 2 analysis
		
	// Check if data is already weighted

		if "`_dta[st_wv]'" != "" {
			di in red "weighted data not allowed"
			exit 198
		}

	// Check bootstrap
	
		if `reps' != -1 & `minreps' != -1 {
			di in red "reps() and minreps() can not be used at the same time"
			exit 198
		}
		
		if `reps' == -1 & `minreps' == -1 local reps 20

		if (`reps' < 2 & `minreps' == -1) | (`reps' == -1 & `minreps' < 2) { 
			di in red "insufficient bootstrap replications"
			exit 198
		}	
			
	// Check number of observations 

		quietly count if `touse'
		if r(N) < 4 { 
			di in red "insufficient observations"
			exit 2001
		}

	// Get exposure variable name
 
		local firstname : word 1 of `varlist'
	    local xname "`varlist'"
		local ixlist "`varlist'"

		
		if regexm("`firstname'","_I") == 1 {
				local varlab: variable label `firstname'
				local poseq =  strpos("`varlab'","=")
				local xname = substr("`varlab'",1,`poseq'-1)
				local ixlist "`varlist'"
		}
		else {		
			if `: word count `varlist'' > 1 {
					di as err "specify just one variable name"
					exit 198
			}
			if `: word count `varlist'' == 1 { 
				qui tab `varlist' 
				if r(r) > 2 {
						di as err "`varlist' has more than two levels, use the xi: notation"
						exit 198
				}
				local xname "`varlist'"
			}
			
		} 
	
		if `overall' == 1  {
						tempvar const1
						gen `const1' =  1
						local xname "`const1'"
						local ixlist "`xname'"
		}
	
	// Recode exposure as 1,2,3...
	
		tempvar expr
		qui gen `expr' = . 
		local j = 0
		qui levelsof `xname', local(exposurel)
		capture confirm numeric variable `xname'
		if !_rc {
			foreach name of local exposurel  {
				qui replace `expr' = `j' if `xname' == `name'
				local j = `j' + 1
			}
		}
		else {
			foreach name of local exposurel  {
				qui replace `expr' = `j' if `xname' == "`name'"
				local j = `j' + 1
			}
		}
		
	// Identify the reference level
	
		qui su `xname' 
		local bl = r(min)

	// More than 2 levels 
		
		if "`ixlist'" != "" {
			if "`: char `xname'[omit]'" != "" local bl : char `xname'[omit]
		}	
	 	
	// Check specified quantiles 
	
		tempname quantlist
		SetQ `quantiles'
		local listquant "`r(quants)'"
		tokenize "`r(quants)'"
		local nq 1
		while "``nq''" != "" {
			local q`nq' ``nq''
			local nq = `nq' + 1
		}
		local nq = `nq' - 1
	
	// Check quantiles and get a list 

		local pow = 10
		local done = 0
		while !`done' {
			local pow = `pow'*10
			local done = 1
			forvalues k=1/`nq' {
				local q = (`q`k'')*`pow'
				local done = `done' & /* 
					*/ (reldif(`q',trunc(`q'))<c(epsfloat))
			}
		}
		local f = trunc(log10(`pow'))
		forvalues k=1/`nq' {
			local myeq = "q" + /*
				*/ string(trunc((`q`k'')*`pow'),"%0`f'.0f")
			local eqnames `eqnames' `myeq'
		}
		local eqnames : list uniq eqnames
		local k : word count `eqnames'
		if `k' < `nq' {
			di as err "only `k' out of `nq' quantiles are " /*
			 */ "unique within a relative precision of c(epsfloat)"
			exit 498
		}
		
	// Define labels for the regression coefficients for each quantile 
	
		tokenize  "`varlist'"

		local eclist ""
		local k 1 
		while `k' <= `nq' {
			local myeq : word `k' of `eqnames'
			local eclist "`eclist' `myeq'"
			local i 1
			while "``i''" != "" {
				local eqnams "`eqnams' `myeq'"
				local conams "`conams' ``i''"
				tempvar v
				local vl "`vl' `v'"
				local i = `i' + 1
			}
 			local eqnams "`eqnams' `myeq'"
			local conams "`conams' _cons"
			 
			tempvar v
			local vl "`vl' `v'"
			local k = `k' + 1
		}
		
	
	tempname coefb
	
	local rs : word count `eqnams'
	mat `coefb' = J(1, `rs', .)
	
	forv s = 1/`rs' {
		local result "`result' `coefb'[1,`s']"
	}
	
	// Get a clean list of percentiles (25, 50, ...) rather than quantiles (.25, .50, ...) 
	
	local listquant2 : subinstr local eqnames "q" "", all  
	  
	// Get name of surival time and failure 
	
	local failure : char _dta[st_bd]
	local depname : char _dta[st_bt] 
    local idname : char _dta[st_id]
	preserve

	quietly keep if `touse'
	quietly keep if _st != 0 
	quietly keep _* `xname'  `wgtv' `depname'  `failure' `expr' `idname'
 
	local nobs = `c(N)'
	local tdf = `c(N)' - (`:word count `varlist'' + 1)
			
	tempfile BOOTRES BOOTJ
	tempname coefs VCE  handle 
		
	qui save `BOOTJ', replace // To be used in the bootstrap loop

	// Point estimates on original sample
	stqkm_f , expvar(`expr' `xname') quantilesf(`listquant2') wgtvv(`wgtv') base(`bl')  nexp(`ixlist') boot(0)
 	mat `coefs' = e(_beta)
	mat colnames `coefs' = `conams'
	mat coleq `coefs' = `eqnams'
	
	// Confidence intervals on bootstrap samples 
 
	*`log' di in gr "(bootstrapping " _c // Dots by Nicola
	 qui postfile `handle' `vl' using "`BOOTRES'", double
 	
	// set the seed
	if "`seed'" != "" {
		set seed `seed'
	}
	local seed `c(seed)' 	
	
	local j 1
	
	// Reps
	if `reps' != -1 {
	
	di ""
	if "`nodots'" == "" nois _dots 0, title(Bootstrap replications (`reps')) // reps(`reps')
	
		while `j'<=`reps' { 
				
			qui use `BOOTJ', clear // Necessary to sample every time from the original dataset
		
			bsample 

			local listquj : subinstr local eqnames "q" "", all  
					
			stqkm_f , expvar(`expr' `xname') quantilesf(`listquj') wgtvv(`wgtv') base(`bl') nexp(`ixlist')
			mat `coefb' = e(_beta)
			version 6: post `handle' `result'
			*if e(mp) == 0 `log' di in gr "." _c // Dots by Nicola
			*else `log' di in red "x" _c // Dots by Nicola
			 
			if "`nodots'" == "" { 
				if e(mp) == 0 nois _dots `j' 0
				else nois _dots `j' 1
			}
				 
			local j = `j' + 1

		}
		
	if "`nodots'" == "" di _n
	}
	
	// MINReps
	local nsuccess 0
	
	if `minreps' != -1 {
	
	di ""
	nois _dots 0, title(Bootstrapping...)
	
		while `nsuccess' < `minreps' {
		
			qui use `BOOTJ', clear
		
			bsample 

			local listquj : subinstr local eqnames "q" "", all  
					
			stqkm_f , expvar(`expr' `xname') quantilesf(`listquj') wgtvv(`wgtv') base(`bl') nexp(`ixlist')
			mat `coefb' = e(_beta)
			version 6: post `handle' `result'
			 
			if e(mp) == 0 local fail 0
			else local fail 1
			
			local nsuccess = `nsuccess' + (`fail' == 0)
				 
			nois _dots `j' `fail'
			
			local j = `j' + 1
		}
		
	local reps = `j' - 1
	
	}
	
	postclose `handle'

	qui use "`BOOTRES'", clear
	*qui bstat `vl' , stat(`coefs') n(`nobs')	  
	*estat bootstrap,  all
	*mat list e(ci_normal)
	*mat list e(ci_percentile)
	
	tempvar missingpercentiles
	qui egen `missingpercentiles' = rowmiss(`vl')
	qui count if `missingpercentiles' != 0
	local missperc = `r(N)'
	local reps = `reps' - `missperc'
	quietly mat accum `VCE' = `vl', dev nocons
	mat rownames `VCE' = `conams'
	mat roweq `VCE' = `eqnams'
	mat colnames `VCE' = `conams'
	mat coleq `VCE' = `eqnams'
	mat `VCE'=`VCE'*(1/(`reps'-1))	 
	version 6: est post `coefs' `VCE', obs(`c(N)')  /*dof(`tdf')*/  depn(_t)
	
	*`log' noi di in gr ")" // Dots by Nicola
	restore
	version 6: est repost, esample(`touse')
	capture erase "`BOOTRES'"
	version 11	
	
	ereturn scalar reps = `reps'
	ereturn scalar N = `nobs'
 
	// Display header 
	
	if `nq' == 1 di in gr "Quantiles from Kaplan-Meier" _col(54) _c
	else di in gr "Simultaneous Quantiles from Kaplan-Meier" _col(54) _c
	di in gr "Number of obs =" in ye %10.0g e(N)

	
	// Save results
	
	ereturn scalar reps = `reps'
	ereturn scalar N = `nobs'
	ereturn scalar n_q = `nq'
	ereturn local eqnames "`eqnames'"
	ereturn local vce  "bootstrap"
	ereturn local vcetype "Bootstrap"
	ereturn local predict "sqreg_p"
	ereturn local cmdline `"stqkm `cmdline'"'
	ereturn local cmd "stqkm"
		
	// Display results
	
	ereturn display, level(`level')
	
	// Warning message
	
	if `missperc' != 0 {
	di as text "Note: one or more quantile estimate could not be estimated in `missperc' bootstrap"
	di as text "replicates; standard errors estimates include only complete replications."
	}  
	
	// Restore original stset
	qui streset [pw = 1]
	
end



*******************
**				 **
**	SUBPROGRAMS  **
**				 **
*******************

capture program drop stqkm_f
program stqkm_f, eclass
version 11
syntax , expvar(string) quantilesf(numlist) base(string) nexp(string) [ wgtvv(string) boot(string) ]

tokenize `expvar'
tempvar exposure origexp 
gen `exposure' = `1'
gen `origexp'  = `2'
 
	qui tab `exposure'
	local levx = r(r)

    
	if "`wgtvv'" != "" {
		qui streset  [pw = `wgtvv']  	
	}
 
	* Get the percentiles of the Kaplan-Meier 

	tempvar atrisk
	qui gen double `atrisk' = _t - _t0  
	
	tempvar flag 
	qui bysort `exposure':  gen `flag' = _n == 1
	 
	tempvar grp subuse  
	qui bysort  `exposure': gen long `grp'=1 if _n==1  
	qui replace `grp'=sum(`grp')  
	local ng = `grp'[_N]

    foreach k of local quantilesf {
			tempvar out`k'
			qui gen `out`k'' = . 
			label var `out`k'' "q`k'" 
			char `out`k'' [varname] "q`k'" 
			local listvarquant "`listvarquant' `out`k''"
	}
 

	local i 1
	while `i' <= `ng' {			
			    capture drop `subuse'
				qui gen byte `subuse'= `grp'==`i' 
				
				DoStats `subuse' `atrisk' `exposure' , percentiles(`quantilesf')
				
				foreach k of local quantilesf {

					if "`boot'" == "0" { // Check if all quantiles are estimable in the original dataset; Andrea 12sep2011
							if r(p`k') == . {
								di as err "Quantile `k'% does not exist in at least one of the exposure groups"
								exit 198								
							}
						
					
					}
				
					qui replace `out`k'' = r(p`k') if `grp'==`i'

				}
				
				local i = `i' + 1
	}
	
	* list `origexp'  `listvarquant' if `flag' == 1, clean noobs subvarname

	// Save differences among percentiles as vectors 

	tempname _bacc
	char `origexp'[omit] `base'
	local s = 1
	foreach v of local listvarquant {
			if `levx' > 1  qui reg `v' `nexp' // basta passare i nomi e non usare i.
			// Check missing values 			
			local names : colfullnames e(b)
			 
			foreach x of local names {
				if substr("`x'",1,2)=="o." local missingp = substr("`x'",3,length("`x'")-2)
			}
			else qui reg `v'  
			
			tempname matclean
			mat `matclean' = e(b)

			if "`missingp'" != "" mat `matclean'[1, colnumb(`matclean',"`missingp'")]  = .  
			
			if `s' == 1 {
					mat `_bacc' = `matclean'
					}
			else  {
					mat `_bacc' = `_bacc' , `matclean'
			}
			
		local `++s'
	}	
	
	if "`missingp'" != "" ereturn scalar mp = 1
	else ereturn scalar  mp   = 0
	ereturn matrix _beta = `_bacc'
	
end

capture program drop DoStats
program define DoStats, rclass /* touse atrisk by */
	syntax varlist  , percentiles(string)
	
	tempvar touse atrisk by
	tokenize `varlist'
    qui gen `touse' = `1'
	qui gen `atrisk' = `2'
	qui gen `by' = `3'
	 
	local id : char _dta[st_id]
	local wv  : char _dta[st_wv]
 
	tempname tatr m s j0
	quietly {
		if `"`wv'"'==`""' {
			summ `atrisk' if `touse'
			scalar `tatr' = r(sum)
			summ _d if `touse'
			local fail = r(sum)
			if `"`id'"'==`""' {
				count if `touse'
				local nsubj = r(N)
			}
			else {
				sort `touse' `id'
				by `touse' `id': /*
				*/ gen byte `m'=1 if _n==1 & `touse'
				summ `m'
				local nsubj = r(sum)
			}
		}
		else {
			tempvar z
			gen double `z' = `wv'*`atrisk'
			summ `z' if `touse'
			scalar `tatr' = r(sum)
			replace `z' = `wv'*_d
			summ `z' if `touse'
			local fail = r(sum)
			drop `z'
			if `"`id'"'==`""' {
				summ `wv' if `touse'
				local nsubj = r(sum)
			}
			else {
				sort `touse' `id'
				by `touse' `id': /*
				*/ gen double `z'=`wv' if _n==1 & `touse'
				summ `z' if `touse'
				local nsubj = r(sum)
			}
		}

		local ttl `"$S_1"'
		capture GetS `touse' `s'
		if _rc==0 { 
			replace `touse' = 0 if `s'>=.
			sort `touse' _t 
			gen long `j0' = _n if `touse'==0
			summ `j0'
			local j0 = cond(r(max)>=.,1,r(max)+1)
			
			foreach k of local percentiles { 
				  Findptl `=(100-`k')/100' `s' `j0'
				  local p`k' $S_1
			}
		}
		else {
		
				foreach k of local percentiles { 			  
				local p`k' .
			}
		}
	}

	
    noisily foreach k of local percentiles {
	 
	ret scalar p`k' = `p`k''
	}
	
	
end
 
 
capture program drop SetQ
program define SetQ /* <nothing> | # [,] # ... */ , rclass
	version 6.0
		if "`*'"=="" {
		ret local quants ".5"
		exit
	}
	local orig "`*'"
	tokenize "`*'", parse(" ,")

	while "`1'" != "" {
		FixNumb "`orig'" `1'
		ret local quants "`return(quants)' `r(q)'"
		mac shift 
		if "`1'"=="," {
			mac shift
		}
	}
end

capture program drop FixNumb
program define FixNumb /* # */ , rclass
	version 6.0
	local orig "`1'"
	mac shift
	capture confirm number `1'
	if _rc {
		Invalid "`orig'" "`1' not a number"
	}
	if `1' >= 1 {
		ret local q = `1'/100
	}
	else 	ret local q `1'
	if `return(q)'<=0 | `return(q)'>=1 {
		Invalid "`orig'" "`return(q)' out of range"
	}
end
		
capture program drop Invalid 
program define Invalid /* "<orig>" "<extra>" */
	version 6.0
	di in red "quantiles(`1') invalid"
	if "`2'" != "" {
		di in red "`2'"
	}
	exit 198
end

capture program drop Findptl
program define Findptl /* percentile s `j0' */
	args p s j0

	if `j0'>=. {
		global S_1 . 
		exit
	}

	tempvar j
	quietly {
		gen long  `j' = _n if float(`s')<=float(`p') in `j0'/l
		summ `j' in `j0'/l
		local j=r(min)
	}
	local t : char _dta[st_t]
	global S_1 = _t[`j']
end

capture program drop GetS
program define GetS /* touse s */
	args touse s 

	tempname prior
	capture {
		capture estimate hold `prior'
		stcox if `touse', estimate bases(`s')
	}
	local rc=_rc
	capture estimate unhold `prior'
	exit `rc'
end

exit


webuse catheter, clear
xi: stqkm i.female, q(25 50 75) reps(20)
 
