import dask.dataframe as dd
import pandas as pd

DEMOCRACY_DATA_PATH = "data/temp/pollution.dta"
PIXELS_PATH = "data/temp/emissions_pixels"
DEM_PIXEL_DATA_PATH = "data/temp/dem_pixels"
DEMOCRACY_COLUMNS = ["wbcode2", "year", "dem", "wb_region"] #, "gdp_median"
POLLUTANTS: list[str] = ["co2", "ch4", "n2o"]

def fetch_democracy_data() -> dd.DataFrame:
	democracy_data: pd.DataFrame = pd.read_stata(DEMOCRACY_DATA_PATH)
	democracy_data = democracy_data[DEMOCRACY_COLUMNS]
	democracy_data["wb_region2"] = pd.factorize(democracy_data["wb_region"])[0]
	democracy_data_delayed: dd.DataFrame = dd.from_pandas(democracy_data).drop(columns="wb_region")
	democracy_data_delayed["dem"] = democracy_data_delayed["dem"].astype("bool")
	democracy_data_delayed["wbcode2"] = democracy_data_delayed["wbcode2"].astype("uint8")
	democracy_data_delayed["year"] = (democracy_data_delayed["year"] - 1970).astype("uint8")
	democracy_data_delayed["wb_region2"] = dd.to_numeric(
		democracy_data_delayed["wb_region2"], downcast="unsigned"
	)
	return democracy_data_delayed

def transform_pollutant_data(pixels_data: dd.DataFrame, pollutant: str) -> dd.DataFrame:
	pixels_data[pollutant] = pixels_data[pollutant]
	return pixels_data

def transform_pixels_data(pixels_data: dd.DataFrame) -> dd.DataFrame:
	for pollutant in POLLUTANTS:
		pixels_data = transform_pollutant_data(
			pixels_data, pollutant
		)
	return pixels_data

def fetch_pixels() -> dd.DataFrame:
	pixels_data: dd.DataFrame = dd.read_parquet(PIXELS_PATH)
	pixels_data = transform_pixels_data(pixels_data)
	return pixels_data

def main() -> None:
	democracy_data: dd.DataFrame = fetch_democracy_data()
	pixels_data: dd.DataFrame = fetch_pixels()
	main_data = pixels_data.merge(democracy_data, on=["wbcode2", "year"])
	print(main_data)
	main_data.to_parquet(DEM_PIXEL_DATA_PATH)
	return

if __name__ == "__main__":
	main()
