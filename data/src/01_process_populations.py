import numpy as np
import rasterio
from affine import Affine
import netCDF4 as nc
from rasterio.warp import reproject, Resampling
from rasterio.crs import CRS
from numpy.typing import NDArray
from typing import Any
from joblib import delayed, Parallel
from tqdm import tqdm

POP_DATA_PATH_TEMPLATE: str = "data/input/population/GHS_POP_E{year}_GLOBE_R2023A_4326_30ss_V1_0.tif"
LON_LAT_NC_PATH: str = "data/input/TOTALS_emi_nc_CO2/EDGAR_2024_GHG_CO2_1970_TOTALS_emi.nc"
POP_COARSE_PATH: str = "data/temp/pop_coarse.tif"
LON_LAT_NC_DATASET = nc.Dataset(LON_LAT_NC_PATH)
LON_COARSE = np.array(LON_LAT_NC_DATASET.variables["lon"][:])
LAT_COARSE = np.array(LON_LAT_NC_DATASET.variables["lat"][:])
RES_LON = LON_COARSE[1] - LON_COARSE[0]
RES_LAT = LAT_COARSE[1] - LAT_COARSE[0]
if RES_LAT > 0:
    LAT_MAX = LAT_COARSE[-1]
else:
    LAT_MAX = LAT_COARSE[0]
LON_MIN = LON_COARSE[0]
COARSE_TRANSFORM = Affine(
    RES_LON, 0.0, LON_MIN - RES_LON / 2,
    0.0, -abs(RES_LAT), LAT_MAX + abs(RES_LAT) / 2
)
COARSE_SHAPE = (len(LAT_COARSE), len(LON_COARSE))
NUM_POP_DATA_TO_PROCESS: int = 7

def load_raster(year_index: int) -> tuple[Any, Affine]:
        year: int = 1980 + year_index * 5
        print("Loading and processing raster for", year, "data. This will take several minutes")
        raster_path: str = POP_DATA_PATH_TEMPLATE.format(year=year)
        raster = rasterio.open(raster_path)
        transform: Affine = raster.transform
        return raster, transform

def process_raster(raster: Any, transform: Affine) -> NDArray:
        coarse_raster: NDArray = np.zeros(COARSE_SHAPE)
        print("    Reprojecting...")
        reproject(
                source=raster.read(1),
                destination=coarse_raster,
                src_transform=transform,
                src_crs=CRS.from_epsg(4326),
                dst_transform=COARSE_TRANSFORM,
                dst_crs=CRS.from_epsg(4326),
                resampling=Resampling.sum,
                num_threads=12
        )
        raster.close()
        return coarse_raster

def load_and_process_raster(year_index: int) -> NDArray:
        raster, transform = load_raster(year_index)
        processed_raster: NDArray = process_raster(raster, transform)
        return processed_raster

def compute_interpolated_raster(
        previous_year: NDArray, yearly_difference: NDArray, relative_year: int
) -> NDArray:
        return previous_year + yearly_difference * relative_year

def parallel_compute_interpolated_rasters(
        previous_year: NDArray, yearly_difference: NDArray, n_years: int
) -> list[NDArray]:
        if n_years > 0:
                years_range: range = range(1, n_years)
        else:
                years_range = range(n_years, 0)
        
        interpolated_rasters_generator = Parallel(n_jobs=-1, verbose=0)(
                delayed(compute_interpolated_raster)(
                        previous_year, yearly_difference, relative_year
                ) for relative_year in years_range
        )
        return [raster for raster in interpolated_rasters_generator] #type:ignore

def perform_linear_interpolation(
        base_year: NDArray, previous_year: NDArray, edge: bool = False
) -> tuple[list[NDArray], NDArray]:
        print("Performing linear interpolation")
        interpolated_rasters: list[NDArray] = []
        yearly_difference: NDArray = (base_year - previous_year) / 5
        
        if edge:
                interpolated_rasters += parallel_compute_interpolated_rasters(
                        previous_year, yearly_difference, -5
                )
        
        interpolated_rasters += [previous_year]
        interpolated_rasters += parallel_compute_interpolated_rasters(
                previous_year, yearly_difference, 5
        )
        return interpolated_rasters, base_year        

def interpolate_rasters(base_year_index: int, previous_year_raster: NDArray | None) -> tuple[list[NDArray], NDArray]:
        base_year_raster = load_and_process_raster(base_year_index)

        if previous_year_raster is None:
                previous_year_raster = load_and_process_raster(base_year_index-1)
                return perform_linear_interpolation(
                        base_year_raster, previous_year_raster, edge=True
                )
        
        return perform_linear_interpolation(
                base_year_raster, previous_year_raster
        )

def make_population_rasters() -> list[NDArray]:
        population_rasters: list[NDArray] = []
        previous_year_raster: NDArray | None = None

        for base_year_index in tqdm(range(NUM_POP_DATA_TO_PROCESS), total=NUM_POP_DATA_TO_PROCESS, desc="Creating all interpolated rasters..."):
                interpolated_rasters, previous_year_raster = interpolate_rasters(
                        base_year_index, previous_year_raster
                )
                population_rasters += interpolated_rasters

        if previous_year_raster is not None: population_rasters += [previous_year_raster]
        return population_rasters

def write_rasters_to_disk(population_rasters):
        raster_file = rasterio.open(
                POP_COARSE_PATH, "w", "GTiff", COARSE_SHAPE[1], COARSE_SHAPE[0],
                len(population_rasters), CRS.from_epsg(4326), COARSE_TRANSFORM,
                population_rasters[0].dtype
        )
        
        for band, raster in tqdm(enumerate(population_rasters), desc="Saving bands..."):
                raster_file.write(raster, band + 1)
        
        raster_file.close()
        
def main() -> None:
        population_rasters: list[NDArray] = make_population_rasters()
        write_rasters_to_disk(population_rasters)
        LON_LAT_NC_DATASET.close()
        return

if __name__ == "__main__":
        main()