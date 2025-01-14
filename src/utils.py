from pathlib import Path


class Paths:
    """Class to manage file paths used in the project."""

    DATA = Path("data")
    RASTER = DATA / "13"
    NIDZ = DATA / "geography-sdz2021-esri-shapefile.zip"
    SGDZ = DATA / "SG_DataZoneBdry_2022.zip"
    LSOA = (
        DATA / "Lower_layer_Super_Output_Areas_2021_EW_BFC_V8_8154990398368723939.zip"
    )
    OUTPUT = DATA / "ndvi.parquet"
