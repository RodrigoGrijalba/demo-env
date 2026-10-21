import pandas as pd

RELEVANT_COLUMNS: list[str] = [
        "wbcode", "wbcode2", "year", "dem", "PopulationtotalSPPOPTOTL", "y"
]
RENAME_MAP: dict[str, str] = {
        "NRU": "NAU", "ROU": "ROM", "SCG": "SER", "SGP": "SIN",
        "TWN": "TAW", "ARE": "UAE", "COD": "ZAR", "SRB": "SER",
        "MNE": "SER"
}
DDCG_DATA_PATH: str = "data/input/DDCGdata_final.dta"
WB_REGIONS_PATH: str = "data/input/wb-regions.csv"
DEMOCRACY_DATA_PATH: str = "data/temp/pollution.dta"

def main() -> None:
        ddcg: pd.DataFrame = pd.read_stata(DDCG_DATA_PATH).query("year >= 1970")[RELEVANT_COLUMNS]
        regions: pd.DataFrame = pd.read_csv(WB_REGIONS_PATH)[["Code", "wb_region"]].replace(RENAME_MAP)
        merged: pd.DataFrame = ddcg.merge(regions, left_on="wbcode", right_on="Code").drop_duplicates()
        merged.to_stata(DEMOCRACY_DATA_PATH, write_index=False)
        return

if __name__ == "__main__":
        main()
