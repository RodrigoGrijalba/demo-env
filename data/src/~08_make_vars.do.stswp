/**/
/*Process Aggregate Data*/
/**/

use ${project}/data/temp/agg_pop_gdp.dta, clear
gen gdp_pc = gdp / pop
gen lgdp_pc = log(gdp_pc)
gen co2_pc = co2 / pop
gen lco2_pc = log(co2_pc)

gen temp = lco2_pc if year == 0
bysort wbcode: egen lco2_pc1970 = max(temp)
drop temp
gen temp = lgdp_pc if year == 0
bysort wbcode: egen gdp_1970 = max(temp)
drop temp

gen lco2_pc_base = lco2_pc1970 * year
xtile gdp_median = gdp_1970, nq(2)
gen gdp_base = gdp_1970 * year

tsset wbcode2 year
gen gdp_lag1 = L1.lgdp_pc
gen gdp_lag2 = L2.lgdp_pc
gen gdp_lag3 = L3.lgdp_pc
gen gdp_lag4 = L4.lgdp_pc

compress

save ${project}/data/output/agg_pop_gdp.dta, replace

keep wbcode2 year gdp gdp_pc lgdp_pc gdp_median lco2_pc_base gdp_base gdp_lag*

/**/
/*Process Pixel Data*/
/**/

merge 1:m year wbcode2 using ${project}/data/temp/pixels_co2.dta, nogen

replace pop = 1 if pop == 0

gen co2_pc = co2 / pop
gen ch4_pc = ch4 / pop
gen n2o_pc = n2o / pop
gen lco2_pc = log(co2_pc)
gen lch4_pc = log(ch4_pc)
gen ln2o_pc = log(n2o_pc)

sum co2_pc
gen co2_pc_std = (co2_pc - r(mean)) / r(sd) * 1000000

gen temp = lco2_pc if year == 0
bysort pixel_id: egen lco2_pc1970_pixel = max(temp)
drop temp

gen lco2_pc_base_pixel = lco2_pc1970_pixel * year

drop co2 pop ch4 n2o lco2_pc1970_pixel

tsset pixel_id year
gen co2_lag1 = L1.lco2_pc
gen co2_lag2 = L2.lco2_pc
gen co2_lag3 = L3.lco2_pc
gen co2_lag4 = L4.lco2_pc

compress

save ${project}/data/output/pixels_co2.dta, replace
