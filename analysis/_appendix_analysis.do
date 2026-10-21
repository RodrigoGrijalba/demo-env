/**********************/
/*Clean and set memory*/
/**********************/
clear all
set more off
set scheme s1mono
set maxvar 10000

/******************/
/*Install programs*/
/******************/
/*cap ssc install xtabond2 
cap ssc install xtivreg2 
cap ssc install spmat 
cap ssc install spmack
*/
/*********************/
/*Sets base directory*/
/*********************/
global project "~/Dropbox/demo_env/democracy_outcomes/ours/replication_package"

/*************************/
/*Create auxiliary locals*/
/*************************/
local coef_renames ""
forvalues i = 1/20 {
	local coef_renames `coef_renames' Placebo_`i' = -`i'
	local coef_renames `coef_renames' Effect_`i' = +`i'
}

local placebo_lc
local tree_placebo_lc
local effect_lc
local pixel_effect_lc

forval i =  1/15 {
	local actuali = 16 - `i'
	
	if `i' != 1 {
		local placebo_lc `placebo_lc'  +
	}
	
	local placebo_lc `placebo_lc' Placebo_`actuali' / 15
}

forval i =  1/13 {
	local actuali = 14 - `i'
	
	if `i' != 1 {
		local tree_placebo_lc `tree_placebo_lc'  +
	}
	
	local tree_placebo_lc `tree_placebo_lc' Placebo_`actuali' / 13
}

forval i = 1/20 {
	
	if `i' != 1 {
		local effect_lc `effect_lc'  +
	}
	
	local effect_lc `effect_lc' Effect_`i' / 20
}

forval i = 1/18 {
	
	if `i' != 1 {
		local pixel_effect_lc `pixel_effect_lc'  +
	}
	
	local pixel_effect_lc `pixel_effect_lc' Effect_`i' / 18
}


colorpalette cblind, select(3 4) nograph
local pre_color `r(p1)'
local post_color `r(p2)'

/****************************************/
/*Program for presenting the DCH results*/
/****************************************/
program define make_dch_table, eclass
syntax, name(name) [outcome(varname) pre(integer 15) post(integer 20)]
scalar pre = `pre'
scalar post = `post'

if "`outcome'" == "" {
	local varlist outcome outcome
}

local N = e(N_avg_total_effect)
local placebo_pval = e(p_jointplacebo)
qui sum `outcome' if _did_sample == 1 & _effect == .
local mdv = r(mean)

mata {
	pre = st_numscalar("pre")
	post = st_numscalar("post")
	avg_effect_index = post + pre + 1
	st_matrix("b", st_matrix("e(b)")[avg_effect_index])
	st_matrix("V", st_matrix("e(V)")[avg_effect_index, avg_effect_index])
}

local cnames "avg_total"

matrix colnames b = `cnames'
matrix colnames V = `cnames'
matrix rownames V = `cnames'

ereturn post b V, obs(`N')
ereturn scalar placebo_pval = `placebo_pval'
ereturn local ctry_fe Y
ereturn local year_fe Y
ereturn local regtrends Y
ereturn scalar mdv = `mdv'
ereturn local cmd "eventstudy_table"
estimates store `name'
end

program define add_dch_omitted, eclass
syntax, [pre(integer 15) post(integer 20)]
scalar pre = `pre'
scalar post = `post'
local N = e(N_avg_total_effect)

mata {
	pre = st_numscalar("pre")
	post = st_numscalar("post")
	first_pre_index = post + pre
	index = first_pre_index..21, 1..20

	b0 = st_matrix("e(b)")[index]
	V0 = st_matrix("e(V)")[index, index]

	num_coefs = cols(b0)
	ommited_coef_pos = pre + 1

	b = J(1, num_coefs+1, 0)
	b[1, 1..ommited_coef_pos-1] = b0[1..ommited_coef_pos-1]
	b[1, ommited_coef_pos+1..num_coefs+1] = b0[ommited_coef_pos..num_coefs]

	V = J(num_coefs+1, num_coefs+1, 0)
	V[1..ommited_coef_pos-1, 1..ommited_coef_pos-1] = V0[1..ommited_coef_pos-1, 1..ommited_coef_pos-1]
	V[ommited_coef_pos+1..num_coefs+1, ommited_coef_pos+1..num_coefs+1] = V0[ommited_coef_pos..num_coefs, ommited_coef_pos..num_coefs]
	V[ommited_coef_pos+1..num_coefs+1, 1..ommited_coef_pos-1] = V0[ommited_coef_pos..num_coefs, 1..ommited_coef_pos-1]
	V[1..ommited_coef_pos-1, ommited_coef_pos+1..num_coefs+1] = V0[1..ommited_coef_pos-1, ommited_coef_pos..num_coefs]

	st_matrix("b", b)
	st_matrix("V", V)

}

local cnames

forval i =  1/`pre' {
	local actuali = `pre' + 1 - `i'
	local cnames `cnames' "Placebo_`actuali'"
}

local cnames `cnames' "0"

forval i = 1/`post' {
	local cnames `cnames' "Effect_`i'"
}


matrix colnames b = `cnames'
matrix colnames V = `cnames'
matrix rownames V = `cnames'

ereturn post b V, obs(`N') findomitted
ereturn local cmd "eventstudy_table"
end


/**********************************/
/*Program for rotating DCH results*/
/**********************************/

program define make_rotated_coefs, eclass
syntax, [pre(integer 15) post(integer 20)]
scalar pre = `pre'
scalar post = `post'
local N = e(N)

mata {
	pre = st_numscalar("pre")
	post = st_numscalar("post")
	b0 = st_matrix("e(b)")
	V0 = st_matrix("e(V)")
	t = (-pre:: 0) \ (1::1::post)
	pre_idx = 1..pre
	y = b0[pre_idx]'
	X = t[pre_idx]
	beta = invsym(X'X) * X'y
	trend = t * beta
	b = b0 - trend'

	st_matrix("b", b)
	st_matrix("V", V)

}

local cnames

forval i =  1/`pre' {
	local actuali = `pre' + 1 - `i'
	local cnames `cnames' "Placebo_`actuali'"
}

local cnames `cnames' "0"

forval i = 1/`post' {
	local cnames `cnames' "Effect_`i'"
}

matrix colnames b = `cnames'
matrix colnames V = `cnames'
matrix rownames V = `cnames'

ereturn post b V, obs(`N') findomitted
ereturn local cmd "eventstudy_table"
end

program define get_rotated_avgs, eclass
syntax , [outcome(varname) pre(integer 15) post(integer 20)]
scalar pre = `pre'
scalar post = `post'
local N = e(N)
if "`outcome'" == "" {
	local varlist outcome outcome
}
qui sum `outcome' if _did_sample == 1 & _effect == .
local mdv = r(mean)
local placebo_test
local effect_lc

forval i =  1/`pre' {	
	local placebo_test `placebo_test' Placebo_`i'
}

test `placebo_test'
local placebo_pval = r(p)

forval i = 1/`post' {
	
	if `i' != 1 {
		local effect_lc `effect_lc'  +
	}
	
	local effect_lc `effect_lc' Effect_`i' / `post'
}

lincom `effect_lc'
matrix b = r(estimate)
matrix V = r(se) ^ 2
matrix colnames b = "avg_total"
matrix colnames V = "avg_total"
matrix rownames V = "avg_total"

ereturn post b V, obs(`N')
ereturn local ctry_fe Y
ereturn local year_fe Y
ereturn local regtrends Y
ereturn scalar mdv = `mdv'
ereturn scalar placebo_pval = `placebo_pval'
ereturn local cmd "rot_avg_table"
end

/**************************************/
/*Program for parsing negative weights*/
/**************************************/
program define make_negweights_table, eclass
mata {
	st_local("N", strofreal(st_matrix("e(M)")[3, 1]))
	b = st_matrix("e(M)")[(1::2), 1]'
	V = J(2, 2, 0)
	V[1, 1] = st_matrix("e(M)")[1, 2] ^ 2
	V[2, 2] = st_matrix("e(M)")[2, 2] ^ 2
	st_matrix("b", b)
	st_matrix("V", V)
}

matrix colnames b = "numpos" "numneg"
matrix colnames V = "numpos" "numneg"
matrix rownames V = "numpos" "numneg"

local lb_se_te = `e(lb_se_te)'
local lb_se_te2 = `e(lb_se_te2)'
local beta = `e(beta)'

ereturn post b V, obs(`N')
ereturn scalar beta = `beta'
ereturn scalar lb_se_te = `lb_se_te'
ereturn scalar lb_se_te2 = `lb_se_te2'
ereturn local cmd "negweights"
end

/***********************************************/
/*Program for getting upper CI for drawing text*/
/***********************************************/
program define get_upper_ci, rclass
matrix b = e(b)
matrix V = e(V)
local z = 1.96
local ymax_pre = .
local ymax_post = .
local ymax = .

forvalues i = 10/13 {
	local bi = b[1, `i']
	local sei = sqrt(V[`i', `i'])
	local ci = `bi' + `z' * `sei'
	local ymax_pre = cond(missing(`ymax_pre'), `bi', max(`ymax_pre', `bi'))
	local ymax_pre = max(`ymax_pre', `ci')
}

forvalues i = 19/22 {
	local bi = b[1, `i']
	local sei = sqrt(V[`i', `i'])
	local ci = `bi' + `z' * `sei'
	local ymax_post = cond(missing(`ymax_post'), `bi', max(`ymax_post', `bi'))
	local ymax_post = max(`ymax_post', `ci')
}

forvalues i = 1/34 {
	local bi = b[1, `i']
	local sei = sqrt(V[`i', `i'])
	local ci = `bi' + `z' * `sei'
	local ymax = cond(missing(`ymax'), `bi', max(`ymax', `bi'))
	local ymax = max(`ymax', `ci')
}

local pad = 0.2 * (`ymax' - min(`ymax_pre', `ymax_post'))


local pre_pos = `ymax_pre' + `pad'
local post_pos = `ymax_post' + `pad'
return scalar pre_pos = `pre_pos'
return scalar post_pos = `post_pos'
end


/*****************************/
/*Load data and generate vars*/
/*****************************/
use ${project}/data/processed/pollution.dta, clear
merge 1:1 wbcode year using ${project}/data/processed/tree_cover.dta, nogen ///
	keepusing(treecover) keep(master matched)

*drop gdp1970
gen year_num = year - 1970 + 1
bysort wbcode: egen treated = max(dem)
xtile decile = gdp_pc, nq(10)
gen gdp=log(gdp_pc)
gen outcome=log(1000*co2/PopulationtotalSPPOPTOTL)
gen ch4_nolog = 1000*ch4/PopulationtotalSPPOPTOTL
replace ch4 = log(1000*ch4/PopulationtotalSPPOPTOTL)
gen n2o_nolog = 1000*n2o/PopulationtotalSPPOPTOTL
replace n2o = log(1000*n2o/PopulationtotalSPPOPTOTL)
gen co2_nolog = 1000 * co2 / PopulationtotalSPPOPTOTL
sum co2_nolog
gen co2_std = (co2_nolog - r(mean)) / r(sd)
tsset wbcode2 year
gen pop = PopulationtotalSPPOPTOTL / 1000000

gen gdp_nolog = exp(gdp)
summ gdp_nolog
gen gdpnonlog_std = (gdp_nolog - r(mean)) / r(sd)

summ treecover
gen tc_std = (treecover - r(mean)) / r(sd)
gen ltc = log(treecover)

gen outcomelag1 = l.outcome
gen outcomelag2 = l2.outcome
gen outcomelag3 = l3.outcome
gen outcomelag4 = l4.outcome

gen tclag1 = l.treecover
gen tclag2 = l2.treecover
gen tclag3 = l3.treecover
gen tclag4 = l4.treecover


gen gdplag1 = l.gdp
gen gdplag2 = l2.gdp
gen gdplag3 = l3.gdp
gen gdplag4 = l4.gdp

gen ynonlog = 1000*co2/PopulationtotalSPPOPTOTL

summ ynonlog
gen ynonlog_std = (ynonlog - r(mean)) / r(sd)

gen temp = outcome if year == 1970
bysort wbcode: egen outcome1970 = max(temp)
drop temp
/*
gen temp = gdp if year == 1970
bysort wbcode: egen gdp1970 = max(temp)
*/
xtile outcome_median = outcome1970, nq(2)
xtile gdp_median = gdp1970, nq(2)


/// Need to generate a variable that identifies treated or not treated observations
bys wbcode2: egen switcher = max(dem)

// Make panelview
gen temp = 1 if dem == . & year == 1970
bysort wbcode: egen startmissing = min(temp)
replace startmissing = 0 if startmissing == .
drop temp


panelview dem dem startmissing if startmissing == 1, i(wbcode) t(year) type(treat) ///
	mycolor(cblind) ytitle("World Bank country code") xtitle(Year) bytiming ///
	legend(off) xlabdist(5) leavegap
graph export "${project}/results/figures/_app_panelview_wostart.png", replace width(1600) height(2000)

panelview dem dem startmissing if startmissing == 0, i(wbcode) t(year) type(treat) ///
	ylabel(, labsize(*.45)) mycolor(cblind) ytitle("World Bank country code") ///
	xtitle(Year) bytiming legend(off) xlabdist(5)
graph export "${project}/results/figures/_app_panelview_wstart.png", replace width(1600) height(2000)

// Testing negative weights and heterogeneity
twowayfeweights outcome wbcode2 year dem, type(feTR) summary_measures
make_negweights_table
est sto e1

twowayfeweights treecover wbcode2 year dem, type(feTR) summary_measures
make_negweights_table
est sto e2

twowayfeweights gdp wbcode2 year dem, type(feTR) summary_measures
make_negweights_table
est sto e3

twowayfeweights outcome wbcode2 year dem, type(feTR) summary_measures ///
	controls(outcomelag1 outcomelag2 outcomelag3 outcomelag4)
make_negweights_table
est sto e4

twowayfeweights treecover wbcode2 year dem, type(feTR) summary_measures ///
	controls(tclag1 tclag2 tclag3 tclag4)
make_negweights_table
est sto e5

twowayfeweights gdp wbcode2 year dem, type(feTR) summary_measures ///
	controls(gdplag1 gdplag2 gdplag3 gdplag4)
make_negweights_table
est sto e6

#delimit ;
esttab e1 e2 e3 e4 e5 e6 using "${project}/results/tables/_app_negweigh.tex", replace booktabs
se varlabels(numpos "N. of ATT w/ positive weights" 
numneg "N. of ATT w/ negative weights") nostar
stats(N beta lb_se_te lb_se_te2, fmt(0 4 4 4)
label("Observations" "TWFE estimate" "Min. s.e. for opposite sign ATT" 
"Min. s.e. for opposite sign in all country-years"))
prehead( "\begin{tabular}{l*{6}{c}}" \hline 
"& \multicolumn{3}{c}{No Lags} & \multicolumn{3}{c}{Four Lags} \\" 
"\cmidrule(lr){2-4} \cmidrule(lr){5-7}") mtitles("CO\textsubscript{2}" "TC\%"
"GDP" "CO\textsubscript{2}" "TC\%" "GDP")
posthead(\hline) postfoot(\hline \end{tabular})
;
#delimit cr


// CH4 and N2O
did_multiplegt_dyn ch4 wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on CH4) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(-.085 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-.19 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_ch4.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on CH4) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.09 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-.175 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_ch4_rot.png", replace
est restore temp
make_dch_table, name(e1) outcome(ch4_nolog)

did_multiplegt_dyn n2o wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on N2O) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(-.115 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-.145 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_n2o.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on N2O) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(-.115 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-.145 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_n2o_rot.png", replace
est restore temp
make_dch_table, name(e2) outcome(n2o_nolog)

#delimit ;
esttab e1 e2 using "${project}/results/tables/_app_pollutants.tex", replace booktabs
se varlabels(avg_total "Avg. Cumm. Effect of Dem.") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends, fmt(0 4 0 0 4)
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"))
prehead( "\begin{tabular}{l*{2}{c}}" \hline ) 
mtitles("Log-Tonnes CH\textsubscript{4} p.c." "Log-Tonnes N\textsubscript{2}O p.c.")
posthead(\hline)
postfoot(\hline 
"\multicolumn{3}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr


// Population
did_multiplegt_dyn pop wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on Population) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(5 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(5 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_pop.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on Population) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(5 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(5 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_pop_rot.png", replace
get_rotated_avgs, outcome(pop)
est sto e1_rot
est restore temp
make_dch_table, name(e1) outcome(pop)
est sto e1_rot

#delimit ;
esttab e1 using "${project}/results/tables/_app_pop.tex", replace booktabs
se varlabels(avg_total "Avg. Cumm. Effect of Dem.") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends, fmt(0 4 0 0 4)
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"))
prehead( "\begin{tabular}{l*{1}{c}}" \hline ) 
mtitles("Population")
posthead(\hline)
postfoot(\hline 
"\multicolumn{2}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr

// Split by median

* CO2
did_multiplegt_dyn outcome wbcode2 year dem if gdp_median == 1, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on CO2) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.32 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.22 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_outcome50_low.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on CO2) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.32 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.22 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_outcome50_lowrot.png", replace
get_rotated_avgs, outcome(co2_nolog)
est sto e1_rot
est restore temp
make_dch_table, name(e1) outcome(co2_nolog)

did_multiplegt_dyn outcome wbcode2 year dem if gdp_median == 2, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on CO2) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.4 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.285 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_outcome50_high.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on CO2) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.4 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.285 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_outcome50_highrot.png", replace
get_rotated_avgs, outcome(co2_nolog)
est sto e2_rot
est restore temp
make_dch_table, name(e2) outcome(co2_nolog)

* GDP
did_multiplegt_dyn gdp wbcode2 year dem if gdp_median == 1, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on GDP) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.25 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.2 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_gdp50_low.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on GDP) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.25 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.2 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_gdp50_lowrot.png", replace
get_rotated_avgs, outcome(gdp_nolog)
est sto e3_rot
est restore temp
make_dch_table, name(e3) outcome(gdp_nolog)

did_multiplegt_dyn gdp wbcode2 year dem if gdp_median == 2, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on GDP) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.235 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.17 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_gdp50_high.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on GDP) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.235 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.17 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_gdp50_highrot.png", replace
get_rotated_avgs, outcome(gdp_nolog)
est sto e4_rot
est restore temp
make_dch_table, name(e4) outcome(gdp_nolog)

* Forest Cover
	
did_multiplegt_dyn treecover wbcode2 year dem if gdp_median == 1, effects(20) placebo(13) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(14, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on Tree Cover) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 14, lc("`pre_color'") ///
	text(5 9.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 14 `avg_post' 34.5, lc("`post_color'") text(5 18.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_tc50_low.png", replace
make_rotated_coefs, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(14, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on Tree Cover) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 14, lc("`pre_color'") ///
	text(5 9.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 14 `avg_post' 34.5, lc("`post_color'") text(5 18.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_tc50_lowrot.png", replace
get_rotated_avgs, outcome(treecover) pre(13)
est sto e5_rot
est restore temp
make_dch_table, name(e5) outcome(treecover) pre(13)

did_multiplegt_dyn treecover wbcode2 year dem if gdp_median == 2, effects(20) placebo(13) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(14, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on Tree Cover) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 14, lc("`pre_color'") ///
	text(6 9.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 14 `avg_post' 34.5, lc("`post_color'") text(6 18.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_tc50_high.png", replace
make_rotated_coefs, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(14, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on Tree Cover) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 14, lc("`pre_color'") ///
	text(6 9.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 14 `avg_post' 34.5, lc("`post_color'") text(6 18.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_tc50_highrot.png", replace
get_rotated_avgs, outcome(treecover) pre(13)
est sto e6_rot
est restore temp
make_dch_table, name(e6) outcome(treecover) pre(13)

#delimit ;
esttab e1 e2 e5 e6 e3 e4 using "${project}/results/tables/_app_split.tex", replace booktabs
se varlabels(avg_total "\textbf{Avg. Cumm. Effect of Dem.}") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends, fmt(0 4 0 0 4) 
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"))
prehead( "\begin{tabular}{l*{6}{c}}" \hline 
"& \multicolumn{2}{c}{CO\textsubscript{2}} & \multicolumn{2}{c}{TC\%} & \multicolumn{2}{c}{GDP} \\" 
"\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7}") 
mtitles( "Low GDP" "High GDP" "Low GDP" 
"High GDP" "Low GDP" "High GDP")
posthead(\hline)
postfoot(\hline 
"\multicolumn{7}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr

#delimit ;
esttab e1_rot e2_rot e5_rot e6_rot e3_rot e4_rot using "${project}/results/tables/_app_split_rot.tex", replace booktabs
se varlabels(avg_total "Avg. Rotated Post-switch Coef.") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends, fmt(0 4 0 0 4) 
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"))
prehead( "\begin{tabular}{l*{6}{c}}" \hline 
"& \multicolumn{2}{c}{CO\textsubscript{2}} & \multicolumn{2}{c}{TC\%} & \multicolumn{2}{c}{GDP} \\" 
"\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7}") 
mtitles( "Low GDP" "High GDP" "Low GDP" 
"High GDP" "Low GDP" "High GDP")
posthead(\hline)
postfoot(\hline 
"\multicolumn{7}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr

// Rotation Analysis
* CO2
did_multiplegt_dyn outcome wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
add_dch_omitted
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(.2 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.25 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_co2rot.png", replace
get_rotated_avgs, outcome(co2_nolog)
est sto e1

* GDP
did_multiplegt_dyn gdp wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
add_dch_omitted
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.175 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.27 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_gdprot.png", replace
get_rotated_avgs, outcome(gdp_nolog)
est sto e2

* Tree Cover
did_multiplegt_dyn treecover wbcode2 year dem, effects(20) placebo(13) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
add_dch_omitted, pre(13)
make_rotated_coefs, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(14, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 14, lc("`pre_color'")  ///
	text(1.75 9.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 14 `avg_post' 34.5, lc("`post_color'") text(2.7 18.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_tcrot.png", replace
get_rotated_avgs, pre(13) outcome(treecover)
est sto e3

#delimit ;
esttab e1 e3 e2 using "${project}/results/tables/_app_rot.tex", replace booktabs
se varlabels(avg_total "Avg. Rotated Post-switch Coef.") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends, fmt(0 4 0 0 4)
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"))
prehead( "\begin{tabular}{l*{3}{c}}" \hline ) 
mtitles("Log-Tonnes CO\textsubscript{2} p.c." "Tree Cover \%" "Log-GDP p.c.")
posthead(\hline)
postfoot(\hline 
"\multicolumn{4}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr

// Analysis for Non log outcome

** Tree Cover
did_multiplegt_dyn tc_std wbcode2 year dem, effects(20) placebo(13) ///
	cluster(wbcode2) switchers(in) trends_nonparam(wb_region2) save_sample ///
	graph_off
est sto temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(1.05 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-1.1 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_tcstd.png", replace
est restore temp
make_dch_table, name(e1) outcome(treecover) pre(13)

did_multiplegt_dyn ltc wbcode2 year dem, effects(20) placebo(13) ///
	cluster(wbcode2) switchers(in) trends_nonparam(wb_region2) save_sample ///
	graph_off
est sto temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(1.05 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-1.1 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_tcstd.png", replace
est restore temp
make_dch_table, name(e2) outcome(treecover) pre(13)


#delimit ;
esttab e1 e2 using "${project}/results/tables/_app_tc_nolog.tex", replace booktabs
se varlabels(avg_total "\textbf{Avg. Cumm. Effect of Dem.}") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends, fmt(0 4 0 0 4) 
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"))
prehead( "\begin{tabular}{l*{2}{c}}" \hline ) 
mtitles("Std. TC\%" "Log-TC\%")
posthead(\hline)
postfoot(\hline 
"\multicolumn{3}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr


** GDP
did_multiplegt_dyn gdp_nolog wbcode2 year dem, effects(20) placebo(15) ///
	cluster(wbcode2) switchers(in) trends_nonparam(wb_region2) save_sample ///
	graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(1.05 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-1.1 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_gdpnonlog.png", replace
est restore temp
make_dch_table, name(e1) outcome(gdp_nolog)

did_multiplegt_dyn gdpnonlog_std wbcode2 year dem, effects(20) placebo(15) ///
	cluster(wbcode2) switchers(in) trends_nonparam(wb_region2) save_sample ///
	graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.07 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-.08 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) // note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_stdgdpnonlog.png", replace
est restore temp
make_dch_table, name(e2) outcome(gdpnonlog_std)

#delimit ;
esttab e1 e2 using "${project}/results/tables/_app_gdp_nolog.tex", replace booktabs
se varlabels(avg_total "\textbf{Avg. Cumm. Effect of Dem.}") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends, fmt(0 4 0 0 4) 
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"))
prehead( "\begin{tabular}{l*{3}{c}}" \hline ) 
mtitles("GDP p.c." "Std. GDP p.c." )
posthead(\hline)
postfoot(\hline 
"\multicolumn{3}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr


** CO2
did_multiplegt_dyn ynonlog wbcode2 year dem, effects(20) placebo(15) ///
	cluster(wbcode2) switchers(in) trends_nonparam(wb_region2) save_sample ///
	graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(1.05 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-1.1 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_co2nonlog.png", replace
est restore temp
make_dch_table, name(e1) outcome(ynonlog)

did_multiplegt_dyn ynonlog_std wbcode2 year dem, effects(20) placebo(15) ///
	cluster(wbcode2) switchers(in) trends_nonparam(wb_region2) save_sample ///
	graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.07 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(-.08 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) // note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_app_figure_stdco2nonlog.png", replace
est restore temp
make_dch_table, name(e2) outcome(ynonlog_std)

#delimit ;
esttab e1 e2 using "${project}/results/tables/_app_1_table1.tex", replace booktabs
se varlabels(avg_total "\textbf{Avg. Cumm. Effect of Dem.}") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends, fmt(0 4 0 0 4) 
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"))
prehead( "\begin{tabular}{l*{3}{c}}" \hline ) 
mtitles("Tonnes CO\textsubscript{2} p.c." "Std. Tonnes CO\textsubscript{2} p.c." )
posthead(\hline)
postfoot(\hline 
"\multicolumn{3}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr

// Different specifications using DCH

replace outcome1970 = outcome1970 * year_num
replace gdp1970 = gdp1970 * year_num

* Tree Cover

did_multiplegt_dyn treecover wbcode2 year dem, effects(20) placebo(13) ///
	cluster(wbcode2) switchers(in) save_sample graph_off
est sto temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_noconts_gdp.png", replace
make_rotated_coefs, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_noconts_rot_gdp.png", replace
est restore temp
make_dch_table, name(e1) outcome(treecover) pre(13)
estadd local startco2 "N"
est restore e1
estadd local startgdp "N"
est sto e1
est restore e1
estadd local regtrends2 "N"
est sto e1


did_multiplegt_dyn treecover wbcode2 year dem, controls(outcome1970) effects(20) ///
	placebo(13) cluster(wbcode2) switchers(in) graph_off /// 
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_co2_tc.png", replace
make_rotated_coefs, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_co2_rot_tc.png", replace
est restore temp
make_dch_table, name(e2) outcome(treecover) pre(13)
estadd local startco2 "Y"
est restore e2
estadd local startgdp "N"
est sto e2
est restore e2
estadd local regtrends2 "Y"
est sto e2

did_multiplegt_dyn treecover wbcode2 year dem, controls(gdp1970) effects(20) /// 
	placebo(13) cluster(wbcode2) switchers(in) graph_off /// 
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_gdp_tc.png", replace
make_rotated_coefs, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_gdp_rot_tc.png", replace
est restore temp
make_dch_table, name(e3) outcome(treecover) pre(13)
estadd local startco2 "N"
est restore e3
estadd local startgdp "Y"
est sto e3
est restore e3
estadd local regtrends2 "Y"
est sto e3

did_multiplegt_dyn treecover wbcode2 year dem, controls(gdp1970 outcome1970) /// 
	effects(20) placebo(13) cluster(wbcode2) switchers(in) graph_off ///
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_both_tc.png", replace
make_rotated_coefs, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_both_rot_tc.png", replace
est restore temp
make_dch_table, name(e4) outcome(treecover) pre(13)
estadd local startco2 "Y"
est restore e4
estadd local startgdp "Y"
est sto e4
est restore e4
estadd local regtrends2 "Y"
est sto e4

#delimit ;
esttab e1 e2 e3 e4 using "${project}/results/tables/_app_table_specs_tc.tex", replace booktabs
se varlabels(avg_total "\textbf{Avg. Cumm. Effect of Dem.}") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends2 startco2 startgdp, fmt(0 4 0 0 4) 
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"
"Starting Emissions" "Starting GDP")) prehead( "\begin{tabular}{l*{4}{c}}" \hline ) 
mtitles("No Controls" "1970 CO\textsubscript{2}" "1970 GDP" "Both" )
posthead(\hline)
postfoot(\hline 
"\multicolumn{5}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr



* GDP

did_multiplegt_dyn gdp wbcode2 year dem, effects(20) placebo(15) ///
	cluster(wbcode2) switchers(in) save_sample graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_noconts_gdp.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_noconts_rot_gdp.png", replace
est restore temp
make_dch_table, name(e1) outcome(gdp_nolog)
estadd local startco2 "N"
est restore e1
estadd local startgdp "N"
est sto e1
est restore e1
estadd local regtrends2 "N"
est sto e1


did_multiplegt_dyn gdp wbcode2 year dem, controls(outcome1970) effects(20) ///
	placebo(15) cluster(wbcode2) switchers(in) graph_off /// 
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_co2_gdp.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_co2_rot_gdp.png", replace
est restore temp
make_dch_table, name(e2) outcome(gdp_nolog)
estadd local startco2 "Y"
est restore e2
estadd local startgdp "N"
est sto e2
est restore e2
estadd local regtrends2 "Y"
est sto e2

did_multiplegt_dyn gdp wbcode2 year dem, controls(gdp1970) effects(20) /// 
	placebo(15) cluster(wbcode2) switchers(in) graph_off /// 
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_gdp_gdp.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_gdp_rot_gdp.png", replace
est restore temp
make_dch_table, name(e3) outcome(gdp_nolog)
estadd local startco2 "N"
est restore e3
estadd local startgdp "Y"
est sto e3
est restore e3
estadd local regtrends2 "Y"
est sto e3

did_multiplegt_dyn gdp wbcode2 year dem, controls(gdp1970 outcome1970) /// 
	effects(20) placebo(15) cluster(wbcode2) switchers(in) graph_off ///
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_both_gdp.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_both_rot_gdp.png", replace
est restore temp
make_dch_table, name(e4) outcome(gdp_nolog)
estadd local startco2 "Y"
est restore e4
estadd local startgdp "Y"
est sto e4
est restore e4
estadd local regtrends2 "Y"
est sto e4

#delimit ;
esttab e1 e2 e3 e4 using "${project}/results/tables/_app_table_specs_gdp.tex", replace booktabs
se varlabels(avg_total "\textbf{Avg. Cumm. Effect of Dem.}") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends2 startco2 startgdp, fmt(0 4 0 0 4) 
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"
"Starting Emissions" "Starting GDP")) prehead( "\begin{tabular}{l*{4}{c}}" \hline ) 
mtitles("No Controls" "1970 CO\textsubscript{2}" "1970 GDP" "Both" )
posthead(\hline)
postfoot(\hline 
"\multicolumn{5}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr


* CO2


did_multiplegt_dyn outcome wbcode2 year dem, effects(20) placebo(15) ///
	cluster(wbcode2) switchers(in) save_sample graph_off
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_noconts.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_noconts_rot.png", replace
est restore temp
make_dch_table, name(e1) outcome(ynonlog)
estadd local startco2 "N"
est restore e1
estadd local startgdp "N"
est sto e1
est restore e1
estadd local regtrends2 "N"
est sto e1


did_multiplegt_dyn outcome wbcode2 year dem, controls(outcome1970) effects(20) ///
	placebo(15) cluster(wbcode2) switchers(in) graph_off /// 
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_co2.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_co2_rot.png", replace
est restore temp
make_dch_table, name(e2) outcome(ynonlog)
estadd local startco2 "Y"
est restore e2
estadd local startgdp "N"
est sto e2
est restore e2
estadd local regtrends2 "Y"
est sto e2

did_multiplegt_dyn outcome wbcode2 year dem, controls(gdp1970) effects(20) /// 
	placebo(15) cluster(wbcode2) switchers(in) graph_off /// 
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_gdp.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_gdp_rot.png", replace
est restore temp
make_dch_table, name(e3) outcome(ynonlog)
estadd local startco2 "N"
est restore e3
estadd local startgdp "Y"
est sto e3
est restore e3
estadd local regtrends2 "Y"
est sto e3

did_multiplegt_dyn outcome wbcode2 year dem, controls(gdp1970 outcome1970) /// 
	effects(20) placebo(15) cluster(wbcode2) switchers(in) graph_off ///
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_both.png", replace
make_rotated_coefs
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'")  ///
	text(.5 7 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.5 22 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_int_specs_both_rot.png", replace
est restore temp
make_dch_table, name(e4) outcome(ynonlog)
estadd local startco2 "Y"
est restore e4
estadd local startgdp "Y"
est sto e4
est restore e4
estadd local regtrends2 "Y"
est sto e4

preserve
use ${project}/data/processed/pixels_co2.dta, clear

did_multiplegt_dyn lco2_pc wbcode2 year dem, controls(lco2_pc_base_pixel) /// 
	effects(18) placebo(15) cluster(wbcode2) switchers(in) graph_off ///
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted, post(18)
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `pixel_effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
get_upper_ci
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(`r(pre_pos)' 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 34.5, lc("`post_color'") text(`r(post_pos)' 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_alt_specs_pixel.png", replace
make_rotated_coefs, post(18)
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `pixel_effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
get_upper_ci
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(`r(pre_pos)' 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 34.5, lc("`post_color'") text(`r(post_pos)' 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_alt_specs_pixel_rot.png", replace
est restore temp
make_dch_table, name(e5) outcome(co2_pc) post(18)
estadd local startco2 "N"
est restore e5
estadd local startco2pix "Y"
est sto e5
est restore e5
estadd local startgdp "N"
est sto e5
est restore e5
estadd local regtrends2 "Y"
est sto e5

did_multiplegt_dyn lco2_pc wbcode2 year dem, controls(gdp_base lco2_pc_base_pixel) /// 
	effects(18) placebo(15) cluster(wbcode2) switchers(in) graph_off ///
	trends_nonparam(wb_region2) save_sample
est sto temp
add_dch_omitted, post(18)
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `pixel_effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
get_upper_ci
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(`r(pre_pos)' 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 34.5, lc("`post_color'") text(`r(post_pos)' 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_alt_specs_both_pixel.png", replace
make_rotated_coefs, post(18)
lincom `placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(`r(se)', 0.0001)
lincom `pixel_effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
get_upper_ci
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(16, lc(gray) lp(dash)) xtitle("Years to Democratization") ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 16, lc("`pre_color'") ///
	text(`r(pre_pos)' 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 34.5, lc("`post_color'") text(`r(post_pos)' 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap)) //note("Avg. of pre-switch coef.: `avg_pre' (`avg_pre_se'). Avg. of post-switch coef.: `avg_post' (`avg_post_se')", ring(10) pos(11))
graph export "${project}/results/figures/_alt_specs_both_pixel_rot.png", replace
est restore temp
make_dch_table, name(e6) outcome(co2_pc) post(18)
estadd local startco2 "N"
est restore e6
estadd local startco2pix "Y"
est sto e6
est restore e6
estadd local startgdp "Y"
est sto e6
est restore e6
estadd local regtrends2 "Y"
est sto e6

restore

#delimit ;
esttab e1 e2 e3 e4 e5 e6 using "${project}/results/tables/_app_table_specs.tex", replace booktabs
se varlabels(avg_total "\textbf{Avg. Cumm. Effect of Dem.}") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends2 startco2 startgdp, fmt(0 4 0 0 4) 
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"
"Starting Emissions" "Starting GDP")) prehead( "\begin{tabular}{l*{6}{c}}" \hline ) 
mtitles("No Controls" "1970 CO\textsubscript{2}" "1970 GDP" "Both" 
"1970 CO\textsubscript{2} Pixels" "Both Pixels")
posthead(\hline)
postfoot(\hline 
"\multicolumn{7}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr



/*
// Within-method specifications
levelsof wb_region2, local(regions)

foreach region of local regions {
	gen dumreg`region' = cond(wb_region2 == `region', 1, 0)
	gen trend`region' = year_num * dumreg`region'
}

tsset wbcode2 year
local out_spec y co2_nolog co2_std
local i = 1
local ctry Y
local yr Y
local regtrends Y

foreach yout of local out_spec {
	xtreg `yout' dem yy* trend*, fe r cluster(wbcode2)
	est sto within`i'
	sum `yout' if treated == 1 & dem == 0
	local ymean = r(mean)
	estadd scalar ymean `ymean'
	est restore within`i'
	estadd local ctry `ctry'
	est sto within`i'
	est restore within`i'
	estadd local yr `yr'
	est sto within`i'
	est restore within`i'
	estadd local regtrends `regtrends'
	est sto within`i'
	local i = `i' + 1
}


#delimit ;
esttab within* using "${project}/results/tables/_app_table_ols.tex", replace
booktabs se keep(dem) varlabels(dem "Effect of Democracy") star(* 0.1 ** 0.05 *** 0.01)
stats(N ymean ctry yr regtrends, fmt(0 4 0 0 0)
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends")) 
prehead( "\begin{tabular}{l*{3}{c}}" \hline ) nomtitles 
postfoot(\hline 
"\multicolumn{4}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr
*/
// Honest DiD

mata index = 35..21, 1..20
did_multiplegt_dyn outcome wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) graph_off save_sample switchers(in)
mata st_matrix("b", st_matrix("e(b)")[index])
mata st_matrix("V", st_matrix("e(V)")[index, index])
honestdid, numpre(15) mvec(0.5(0.5)2) coefplot colorspec(black gray) /// 
	ytitle("95% Robust Conf. Int.") xtitle("Extent of Violation of Parallel Trends Assumption")
graph export ${project}/results/figures/_app_honest_co2.jpg, as(jpg) replace

did_multiplegt_dyn gdp wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) graph_off save_sample switchers(in)
mata st_matrix("b", st_matrix("e(b)")[index])
mata st_matrix("V", st_matrix("e(V)")[index, index])
honestdid, numpre(15) mvec(0.5(0.5)2) coefplot colorspec(black gray) /// 
	ytitle("95% Robust Conf. Int.") xtitle("Extent of Violation of Parallel Trends Assumption")
graph export ${project}/results/figures/_app_honest_gdp.jpg, as(jpg) replace

mata index = 33..21, 1..20
did_multiplegt_dyn treecover wbcode2 year dem, effects(20) placebo(13) cluster(wbcode2) ///
	trends_nonparam(wb_region2) graph_off save_sample switchers(in)
mata st_matrix("b", st_matrix("e(b)")[index])
mata st_matrix("V", st_matrix("e(V)")[index, index])
honestdid, numpre(15) mvec(0.5(0.5)2) coefplot colorspec(black gray) /// 
	ytitle("95% Robust Conf. Int.") xtitle("Extent of Violation of Parallel Trends Assumption")
graph export ${project}/results/figures/_app_honest_tc.jpg, as(jpg) replace

/*
//// Working on Stacked dataset

use "${project}/temp/stacked_flex_Full.dta", clear
// the original variable is this one
// replace y=100*log(1000*co2/PopulationtotalSPPOPTOTL)
// This means tons of CO2 per capita

levelsof wb_region2, local(regions)

foreach region of local regions {
	gen dumreg`region' = cond(wb_region2 == `region', 1, 0)
	gen trend`region' = year_num * dumreg`region'
}

char event_time[omit] -1
xi i.event_time
sum event_time

forvalues i=26/71 {
if `i'!=40 {
	gen treat_ev_`i' = _Ievent_tim_`i' * treated
	local lab_value = `i' - 41
	label var treat_ev_`i' "`lab_value'"
}
if `i'==40 {
	gen treat_ev_`i' = 0
	label var treat_ev_40 "-1"
}
}

char bin_time[omit] -1
xi i.bin_time
sum bin_time

forvalues i=6/14 {
if `i'!=8 {
	gen treat_bin_`i' = _Ibin_time_`i' * treated
	local low_value = `i' * 5 - 45
	local high_value = `i' * 5 - 41
	label var treat_bin_`i' "`low_value' to `high_value'"
}
if `i'==8 {
	gen treat_bin_`i' = 0
	label var treat_bin_8 "-5 to -1"
}
}

gen post_event = (event_time >= 0)
gen post_treatment = post_event * treated

local out_spec y co2_nolog co2_std
local i = 1
local ctry Y
local yr Y
local regtrends Y
foreach yout of local out_spec{
	/// Effect of democracy on tons of co2 by thousand population
	reghdfe `yout' post_treatment trend*, abs(wbcode2#cohort year#cohort) cluster(wbcode2)
	est sto reg`i'
	sum `yout' if treated == 1 & post_event == 0
	local ymean = r(mean)
	estadd scalar ymean `ymean'
	est restore reg`i'
	estadd local ctry `ctry'
	est sto reg`i'
	est restore reg`i'
	estadd local yr `yr'
	est sto reg`i'
	est restore reg`i'
	estadd local regtrends `regtrends'
	est sto reg`i'
	local i = `i' + 1
}

#delimit ;
esttab reg* using "${project}/results/tables/_app_stacked_table1.tex", replace booktabs
se keep(post_treatment) varlabels(post_treatment "Effect of Democracy")
stats(N ymean ctry yr regtrends, fmt(0 4 0 0 0)
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends")) 
prehead( "\begin{tabular}{l*{3}{c}}" \hline ) star(* 0.1 ** 0.05 *** 0.01)
mtitles("Log-Tonnes CO\textsubscript{2} p.c." "Tonnes CO\textsubscript{2} p.c."
"Std. Tonnes CO\textsubscript{2} p.c.") posthead(\hline)
postfoot(\hline 
"\multicolumn{4}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr
*/
