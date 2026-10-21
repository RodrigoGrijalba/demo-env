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

colorpalette cblind, select(3 4) nograph
local pre_color `r(p1)'
local post_color `r(p2)'


/****************************************/
/*Program for reordering the DCH results*/
/****************************************/
program define make_event_table, eclass
syntax, name(name) [outcome(varname) pre(integer 15)]
scalar pre = `pre'
if "`outcome'" == "" {
	local varlist outcome outcome
}

local N = e(N_avg_total_effect)
local placebo_pval = e(p_jointplacebo)
qui sum `outcome' if _did_sample == 1 & _effect == .
local mdv = r(mean)

mata {
	pre = st_numscalar("pre")
	avg_effect_index = 20 + pre + 1
	index = 26..21, 1..5, avg_effect_index

	b0 = st_matrix("e(b)")[index]
	V0 = st_matrix("e(V)")[index, index]

	num_coefs = cols(b0)
	ommited_coef_pos = 7

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

forval i = -6/-1 {
	local cnames `cnames' "t=`i'"
}

local cnames `cnames' "t=0 (ommited)"

forval i = 1/5 {
	local cnames `cnames' "t=`i'"
}

local cnames `cnames' "avg_total"

matrix colnames b = `cnames'
matrix colnames V = `cnames'
matrix rownames V = `cnames'

ereturn post b V, obs(`N') findomitted
ereturn scalar placebo_pval = `placebo_pval'
ereturn local ctry_fe Y
ereturn local year_fe Y
ereturn local regtrends Y
ereturn scalar mdv = `mdv'
ereturn local cmd "eventstudy_table"
estimates store `name'
end

program define add_dch_omitted, eclass

local N = e(N_avg_total_effect)
syntax, [pre(integer 15)]
scalar pre = `pre'
mata {
	pre = st_numscalar("pre")
	first_pre_index = 20 + pre
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

forval i = 1/20 {
	local cnames `cnames' "Effect_`i'"
}


matrix colnames b = `cnames'
matrix colnames V = `cnames'
matrix rownames V = `cnames'

ereturn post b V, obs(`N') findomitted
ereturn local cmd "eventstudy_table"
end

/*****************************/
/*Load data and generate vars*/
/*****************************/
use ${project}/data/processed/pollution.dta, clear
gen gdp = log(gdp_pc)
gen gdp_nolog = gdp_pc
gen outcome=log(1000*co2/PopulationtotalSPPOPTOTL)
gen outcome_nolog = 1000*co2/PopulationtotalSPPOPTOTL
tsset wbcode2 year

* CO2
did_multiplegt_dyn outcome wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
make_event_table, name(e1) outcome(outcome_nolog)
est restore temp
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
	text(.225 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.275 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap))
graph export "${project}/results/figures/_main_co2log.png", replace

* GDP
did_multiplegt_dyn gdp wbcode2 year dem, effects(20) placebo(15) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
est sto temp
make_event_table, name(e2) outcome(gdp_nolog)
est restore temp
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
	text(.2 11.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 16 `avg_post' 36.5, lc("`post_color'") text(.22 20.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap))
graph export "${project}/results/figures/_main_gdp.png", replace

* Forest Cover
use ${project}/data/processed/tree_cover.dta, clear

did_multiplegt_dyn treecover wbcode2 year dem, effects(20) placebo(13) cluster(wbcode2) ///
	trends_nonparam(wb_region2) save_sample switchers(in) graph_off
di e(se)
est sto temp
make_event_table, name(e3) outcome(treecover) pre(13)
est restore temp
add_dch_omitted, pre(13)
lincom `tree_placebo_lc'
local avg_pre = round(`r(estimate)', 0.0001)
local avg_pre_se = round(round(`r(se)', 0.0001), 0.0001)
lincom `effect_lc'
local avg_post = round(`r(estimate)', 0.0001)
local avg_post_se = round(`r(se)', 0.0001)
coefplot, vertical omitted xlabel(, labsize(2) angle(90)) yline(0, lc(red)) ///
	xline(14, lc(gray) lp(dash)) xtitle("Years to Democratization") ytitle(Effect on Tree Cover) ///
	rename(`coef_renames') connect(l) addplot(pci `avg_pre' 0 `avg_pre' 14, lc("`pre_color'") ///
	text(2.25 9.5 "Avg. pre-switch:" "`avg_pre'" "(`avg_pre_se')", c(`pre_color')) || ///
	pci `avg_post' 14 `avg_post' 34.5, lc("`post_color'") text(2.75 18.5 "Avg. post-switch:" "`avg_post'" "(`avg_post_se')", c(`post_color'))) ///
	ciopts(recast(rcap))
graph export "${project}/results/figures/_main_forest.png", replace


#delimit ;
esttab e1 e3 e2 using "${project}/results/tables/_main_table.tex", replace booktabs
se varlabels(avg_total "\hline\textbf{Avg. Cumm. Effect of Dem.}") star(* 0.1 ** 0.05 *** 0.01)
stats(N mdv ctry_fe year_fe regtrends placebo_pval, fmt(0 4 0 0 4)
label("Observations" "Mean DV" "Country FE" "Year FE" "Regional Trends"
"Joint p-val of pre-periods")) prehead( "\begin{tabular}{l*{3}{c}}" \hline ) 
mtitles("Log-Tonnes CO\textsubscript{2} p.c." "Tree Cover \%" "GDP p.c.")
posthead(\hline)
postfoot(\hline 
"\multicolumn{4}{l}{\footnotesize \sym{*} \(p<0.10\), \sym{**} \(p<0.05\), \sym{***} \(p<0.01\)}\\"
\end{tabular})
;
#delimit cr



