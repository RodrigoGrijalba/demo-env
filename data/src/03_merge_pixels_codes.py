import numpy as np
import dask.dataframe as dd
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from joblib import Parallel, parallel_config, delayed
from typing import Any

COUNTRY_SHAPES_PATH: str = "data/input/geo_data/shapecountries.dbf"
DEMOCRACY_DATA_PATH: str = "data/input/DDCGdata_final.dta"
LON_LAT_DATA_PATH: str = "data/temp/lon_lat.csv"
PIXEL_COUNTRIES_PATH: str = "data/temp/pixel_countries.csv"

def fetch_country_shapes() -> gpd.GeoDataFrame:
        country_shapes: gpd.GeoDataFrame = gpd.read_file(COUNTRY_SHAPES_PATH)
        country_shapes = country_shapes[["geometry"]]
        country_shapes["_ID"] = np.linspace(1, 246, num=246).astype("uint8")
        return country_shapes

def fetch_valid_country_ids() -> pd.DataFrame:
        democracy_data: pd.DataFrame = pd.read_stata(DEMOCRACY_DATA_PATH)
        valid_country_ids: pd.DataFrame = democracy_data.query("year == 1970")
        valid_country_ids = valid_country_ids[["_ID", "wbcode2"]]
        return valid_country_ids

def filter_valid_shapes() -> gpd.GeoDataFrame:
        country_shapes: gpd.GeoDataFrame = fetch_country_shapes()
        valid_country_ids: pd.DataFrame = fetch_valid_country_ids()
        valid_country_shapes: gpd.GeoDataFrame = country_shapes.merge(
                valid_country_ids, how="right", on="_ID"
        )
        return valid_country_shapes

def get_pixel_lon_lats() -> gpd.GeoDataFrame:
        lon_lat_table: pd.DataFrame = pd.read_csv(LON_LAT_DATA_PATH)
        lon_lat: zip[tuple[Any, Any]] = zip(lon_lat_table["lon"], lon_lat_table["lat"])
        
        with parallel_config(n_jobs=-1, verbose=5):
                geometry = Parallel()(delayed(Point)(lon, lat) for lon, lat in lon_lat)

        pixel_lon_lats: gpd.GeoDataFrame = gpd.GeoDataFrame(
                lon_lat_table, geometry=geometry, crs="EPSG:4326"
        )
        return pixel_lon_lats

def perform_spatial_join(valid_country_shapes) -> gpd.GeoDataFrame:
        pixel_lon_lats: gpd.GeoDataFrame = get_pixel_lon_lats()
        valid_country_pixels: gpd.GeoDataFrame = gpd.sjoin(
                pixel_lon_lats, valid_country_shapes, how="left", predicate="within"
        ).fillna(0)[["wbcode2"]]
        return valid_country_pixels
        
def main() -> None:
        valid_country_shapes: gpd.GeoDataFrame = filter_valid_shapes()
        valid_country_pixels: gpd.GeoDataFrame = perform_spatial_join(
                valid_country_shapes
        )
        valid_country_pixels.to_csv(
                PIXEL_COUNTRIES_PATH,
                index_label = "pixel_id"
        )
        return

if __name__ == "__main__":
        main()