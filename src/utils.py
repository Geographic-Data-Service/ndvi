from pathlib import Path


class Paths:
    """Class to manage file paths used in the project."""

    GISDATA_DIR = Path("gisdata")
    RASTER_FOLDER = GISDATA_DIR / "13"
    NIDZ_SHAPEFILE = GISDATA_DIR / "geography-sdz2021-esri-shapefile.zip"
    SGDZ_SHAPEFILE = GISDATA_DIR / "SG_DataZoneBdry_2022.zip"
    LSOA_BOUNDARIES_SHAPEFILE = (
        GISDATA_DIR
        / "Lower_layer_Super_Output_Areas_2021_EW_BFC_V8_8154990398368723939.zip"
    )
    OUTPUT_PARQUET = GISDATA_DIR / "ndvi.parquet"
