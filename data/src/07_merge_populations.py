import numpy as np
import pandas as pd
import dask.dataframe as dd

POP_DATA_PATH: str = "data/input/population/GHS_POP_E1975_GLOBE_R2023A_4326_30ss_V1_0.tif"
LON_LAT_NC_PATH: str = "data/input/TOTALS_emi_nc_CO2/EDGAR_2024_GHG_CO2_1970_TOTALS_emi.nc"
POP_COARSE_PATH: str = "data/temp/pop_coarse.csv"
DEM_PIXELS_DATA_PATH: str = "data/temp/dem_pixels"
DEMOCRACY_DATA_PATH: str = "data/temp/pollution.dta"
DEMOCRACY_COLUMNS: list[str] = [
        "wbcode", "wbcode2", "year", "dem", "y", "PopulationtotalSPPOPTOTL",
        "wb_region"
] #, "gdp_median"
GDP_MERGE_COLUMNS: list[str] = ["wbcode2", "year", "gdp"]
MAIN_PIXEL_DATA_PATH: str = "data/temp/pixels_co2.dta"
SECONDARY_GDP_DATA_PATH: str = "data/temp/agg_pop_gdp.dta"

def prepare_gdp_data() -> pd.DataFrame:
        gdp_data: pd.DataFrame = pd.read_stata(DEMOCRACY_DATA_PATH)[DEMOCRACY_COLUMNS]
        gdp_data = gdp_data.rename(columns={"y": "gdp"})
        gdp_data["gdp"] = np.exp(gdp_data["gdp"] / 100) * gdp_data["PopulationtotalSPPOPTOTL"]
        gdp_data = gdp_data.drop(columns="PopulationtotalSPPOPTOTL")
        gdp_data["wb_region2"] = pd.factorize(gdp_data["wb_region"])[0]
        gdp_data = dd.from_pandas(gdp_data).drop(columns="wb_region")
        gdp_data["dem"] = gdp_data["dem"].astype("bool")
        gdp_data["wbcode2"] = gdp_data["wbcode2"].astype("int16")
        gdp_data["year"] = (gdp_data["year"] - 1970).astype("int8")
        return gdp_data

def merge_pixels_with_gdp() -> tuple[dd.DataFrame, ...]:
        gdp_data: pd.DataFrame = prepare_gdp_data()
        pixel_co2_data: dd.DataFrame = dd.read_parquet(DEM_PIXELS_DATA_PATH)
        aggregated_pixels: dd.DataFrame = pixel_co2_data.groupby(["wbcode2", "year"])[["pop", "co2"]].sum().reset_index()
        pixel_co2_data = pixel_co2_data.merge(
                gdp_data[GDP_MERGE_COLUMNS], how="left"
        )
        aggregated_pixels = gdp_data.merge(aggregated_pixels, on=["wbcode2", "year"], how="left")
        return pixel_co2_data, aggregated_pixels

def main() -> None:
        pixel_data, country_data = merge_pixels_with_gdp()
        pixel_data.compute().to_stata(MAIN_PIXEL_DATA_PATH, write_index=False)
        country_data.compute().to_stata(SECONDARY_GDP_DATA_PATH, write_index=False)
        return

if __name__ == "__main__": 
        main()