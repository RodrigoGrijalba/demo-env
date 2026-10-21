# Replication Package (draft)

## Contents

- `data`: Folder. Contains all data, as well as the source code to process it for analysis.
    - `raw`: Folder. Contains all raw data, either in identical form to their original sources or extracted directly from the archive acquired from their original sources.
        - `geo_data`: Folder. Contains data for country shapes. Sourced from the replication package in [Acemoglu et al. (2019)](https://www.journals.uchicago.edu/doi/10.1086/700763).
        - `population`: Folder. Contains data for pixel-level population. Each file within is a geotiff raster at a resolution of 30 arcseconds, named with the format "GHS_POP_E\[year\]_GLOBE_R2023A_4326_30ss_V1_0", where "year" between brackets is in place of the corresponding year of the population data. These data are only available in 5 year increments starting from 1975. Sourced from the Global Human Settlement (GHS) population grid multitemporal dataset ([GHS-POP](https://human-settlement.emergency.copernicus.eu/ghs_pop2023.php), by [Schiavina et al., 2023](https://doi.org/10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE)).
        - `TOTAL_EMI_NC_[pollutant]`: Three folders. "pollutant" between brackets is replaced by the name of the corresponding pollutant. Each folder contains contains files named as "EDGAR_2024_GHG_\[pollutant\]_\[year\]_TOTALS_emi", where "year" between brackets is replaced by each file's corresponding year. These files are NetCDF 4 raster files containing data on tonnes of overall emissions of the corresponding pollutant in each corresponding year. Sourced and extractred from the archive files for the Emissions Database for Global Atmospheric Research (EDGAR) Yearly Emissions gridmaps in ton substance for Carbon Dioxide ([link](https://edgar.jrc.ec.europa.eu/gallery?release=v2024ghg&substance=CO2&sector=TOTALS)), Methane ([link](https://edgar.jrc.ec.europa.eu/gallery?release=v2024ghg&substance=CH4&sector=TOTALS)), and Nitrous Oxide ([link](https://edgar.jrc.ec.europa.eu/gallery?release=v2024ghg&substance=N2O&sector=TOTALS)) by [Crippa et al. (2024)](https://dx.doi.org/10.2760/4002897).
        - `[IEA_]EDGAR_[pollutant]_1970_2025`: Three folders. "pollutant" between brackets is replaced by the name of the corresponding pollutant. The prefix "IEA_" between brackets is only present for CO2. Each folder contains contains files named as "\[IEA_\]EDGAR_\[pollutant\]_1070_2025.xlsx". These files are Excel sheets containing data on tonnes of overall emissions of the corresponding pollutant in each corresponding year, by country. Sourced and extractred from the files for the Emissions Database for Global Atmospheric Research (EDGAR) Yearly Emissions gridmaps in ton substance for Carbon Dioxide ([link](https://edgar.jrc.ec.europa.eu/dataset_ghg2026)).
        - `DDCGdata_final.dta`: File. Contains country-year level democracy and GDP data, as prepared by [Acemoglu et al. (2019)](https://www.journals.uchicago.edu/doi/10.1086/700763).
        - `wb-regions.csv`: File. Contains mappings between countries and world geographical regions, as defined by the World Bank ([link](https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups)). Sourced from Our World in Data ([link](https://ourworldindata.org/grapher/world-regions-according-to-the-world-bank))
        - `vcf_ctry_year_treecover.csv`: File. Contains country year tree cover proportions. Tree cover is from the MODIS/Terra Vegetation Continuous Fields Yearly L3 Global 250 m SIN Grid product (MOD44B, version 061), distributed by the NASA Land Processes Distributed Active Archive Center

    - `clean`: Folder. Destination of clean data without processing or merging.
    - `temp`: Folder. Destination of temporary middle steps during the data preparation process.
    - `processed`: Folder. Destination of the final output of the data preparation process, i.e. the files that will be used for analysis.
    - `src`: Folder. Contains the source code for all steps in data processing.
        - `01_process_populations.py`: This file reprojects and aggregates the raw population data at the 30 arc-seconds resolution to be compatible with the emissions data at the 6 arc-minute resolution. It also linearly interpolates the population data for each missing year. Saves the resulting interpolated, coarser resolution data in `data/temp/pop_coarse.tif`.
        - `02_create_lonlat_table.py`: Creates a mapping between each pixel in the coarse population and pollutant datasets, and a corresponding pixel ID, latitude, and longitude. Saves the mapping in `data/temp/lon_lat.csv`
        - `03_merge_pixels_codes.py`: Filters for countries with available democracy indicators, and maps each pixel's ID with their corresponding country. Saves the mapping to `data/temp/pixel_countries.csv`
        - `04_gridmaps_to_table.py`: Assigns each pixel-year the corresponding emission and population measurements. The mapping is saved to `data/temp/emissions_pixels`, which is a collection of parquet part files.
        - `05_make_wb_regions.py`: Maps countries to geographical regions as defined by the World Bank, together with the corresponding yearly democracy and GDP data. This mapping is saved in `data/temp/pollution.dta`
        - `06_merge_pixels_democracy.py`: Merges the country-year level democracy data with the corresponding pixel-year level emission data. The resulting table is saved in `data/temp/dem_pixels`, a collection of partial parquet files.
        - `07_merge_populations.py`: Merges the pixel-year level data on emissions with pixel-year population. It also creates a separate table with the pixel-level population data aggregated at the country-year level, and merged with GDP and democracy. The pixel-year level data are saved in `data/temp/pixels_co2.dta`. The country-year level data are saved in `data/temp/agg_pop_gdp.dta`
        - `08_make_vars.do`: creates ansillary variables for the pixel-year and country-year level tables. Specifically, per-capita counts, logged per-capita counts, normally standardized, and lagged versions of GDP and emissions. The country-year level data are saved in `data/output/agg_pop_gdp.dta`, and the pixel-year level data are saved in `data/output/pixels_co2.dta`
- `analysis`: Folder. Contains the .do files for analysis of the data.
    - `_main.do`: Contains the script for the main analysis in the article.
    - `_appendix.do`: Contains the script for generating the results within the appendix
- `results`: Folder. Contains two subfolders: `figures` and `tables` which are the destinations of the corresponding analysis results, generated by the scripts in `analysis`.

## Running
The code within this replication package has been tested for Python 3.14.3 and Stata 18 MP. The libraries necessary to run the Python script are listed in the `requirements.txt` file. The stata packages necessary to run the Do files are the following:

- panelview
- coefplot
- twowaynegweights
- honestdid
- did_multiplegt_dyn

In a computer with an AMD Ryzen 5 7533HS with 32GB of RAM, this code takes approximately 3 and a half hours to run.