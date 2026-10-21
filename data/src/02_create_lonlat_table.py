import netCDF4 as nc
import pandas as pd
import numpy as np
from numpy.typing import NDArray

LON_LAT_NC_PATH = "data/input/TOTALS_emi_nc_CO2/EDGAR_2024_GHG_CO2_1970_TOTALS_emi.nc"
LON_LAT_CSV_PATH = "data/temp/lon_lat.csv"

def make_lon_lat_array(lon_lat_nc_dataset: nc.Dataset) -> tuple[NDArray, ...]:
	lon_array = np.array(lon_lat_nc_dataset.variables["lon"][:])
	lat_array = np.array(lon_lat_nc_dataset.variables["lat"][:])
	lon_array, lat_array = np.meshgrid(lon_array, lat_array)
	return lon_array.ravel(), lat_array.ravel()

def main() -> None:
	lon_lat_nc_dataset = nc.Dataset(LON_LAT_NC_PATH)
	lon_array, lat_array = make_lon_lat_array(lon_lat_nc_dataset)
	lon_lat_table = pd.DataFrame({"lon": lon_array, "lat": lat_array})
	lon_lat_table.to_csv(LON_LAT_CSV_PATH, index_label="pixel_id")
	lon_lat_nc_dataset.close()
	return

if __name__ == "__main__":
	main()