from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from rtree import index
from scipy.ndimage import median_filter
from scipy.spatial.distance import mahalanobis
from shapely.geometry import box
from tqdm import tqdm


def read_lsoa():
    """Reads and combines LSOA boundary data from multiple sources.

    Returns:
        GeoDataFrame: Combined LSOA boundary data.
    """
    nidz = (
        gpd.read_file(Path("gisdata/geography-dz2021-esri-shapefile.zip"))[
            ["DZ2021_cd", "geometry"]
        ]
        .to_crs(4326)
        .rename(columns={"DZ2021_cd": "LSOA21CD"})
    )
    sgdz = (
        gpd.read_file(
            Path("gisdata/SG_DataZoneBdry_2022.zip"), layer="SG_DataZoneBdry_2022_EoR"
        )[["DZCode", "geometry"]]
        .to_crs(4326)
        .rename(columns={"DZCode": "LSOA21CD"})
    )
    lsoa_boundaries = gpd.read_file(
        Path("gisdata/gov/LSOA2021/LSOA_2021_EW_BFC_V8.shp")
    ).to_crs(4326)[["LSOA21CD", "geometry"]]
    return pd.concat([lsoa_boundaries, sgdz, nidz])


def create_rtree_index(raster_files):
    """Creates an R-tree index for raster files.

    Args:
        raster_files (list): List of raster file paths.

    Returns:
        tuple: R-tree index and a dictionary of raster bounding boxes.
    """
    idx = index.Index()
    raster_bboxes = {}
    for i, raster_file in enumerate(raster_files):
        with rasterio.open(raster_file) as src:
            bounds = src.bounds
            raster_bbox = box(bounds[0], bounds[1], bounds[2], bounds[3])
            idx.insert(i, raster_bbox.bounds)
            raster_bboxes[i] = (raster_bbox, raster_file)
    return idx, raster_bboxes


def find_overlapping_rasters_lsoa(lsoa, rtree_idx, raster_bboxes):
    """Finds raster files overlapping with a given LSOA.

    Args:
        lsoa (GeoSeries): LSOA geometry.
        rtree_idx (Index): R-tree index of raster files.
        raster_bboxes (dict): Dictionary of raster bounding boxes.

    Returns:
        list: List of valid raster file paths.
    """
    lsoa_bbox = lsoa.geometry.bounds
    overlapping_rasters = list(rtree_idx.intersection(lsoa_bbox))

    valid_rasters = []
    for i in overlapping_rasters:
        raster_bbox, raster_file = raster_bboxes[i]
        if lsoa.geometry.intersects(raster_bbox):
            valid_rasters.append(raster_file)
    return valid_rasters


def detect_outliers(reshaped_nonzero):
    """Detects outliers in the data using Mahalanobis distance.

    Args:
        reshaped_nonzero (ndarray): Non-zero reshaped data.

    Returns:
        ndarray: Data with outliers removed.
    """
    mean_vector = np.mean(reshaped_nonzero, axis=0)
    cov_matrix = np.cov(reshaped_nonzero, rowvar=False)
    inv_cov_matrix = np.linalg.inv(cov_matrix)

    distances = np.array(
        [mahalanobis(pixel, mean_vector, inv_cov_matrix) for pixel in reshaped_nonzero]
    )

    threshold = np.percentile(distances, 95)
    outliers = distances > threshold
    return reshaped_nonzero[~outliers]


def calculate_ndvi(red_band, nir_band):
    """Calculates the Normalized Difference Vegetation Index (NDVI).

    Args:
        red_band (ndarray): Red band data.
        nir_band (ndarray): Near-infrared band data.

    Returns:
        ndarray: NDVI values.
    """
    ndvi = np.divide(
        nir_band - red_band,
        nir_band + red_band,
        where=(nir_band + red_band) != 0,
    )
    ndvi[(ndvi > 1) | (ndvi < 0.1)] = 0
    ndvi = ndvi[(~np.isnan(ndvi)) & (ndvi != 0)]
    ndvi = median_filter(ndvi, size=10)
    return ndvi


def calculate_evi(nir_band, red_band, green_band):
    """Calculates the Enhanced Vegetation Index (EVI).

    Args:
        nir_band (ndarray): Near-infrared band data.
        red_band (ndarray): Red band data.
        green_band (ndarray): Green band data.

    Returns:
        ndarray: EVI values.
    """
    G = 2.5
    C1 = 6
    C2 = 7.5
    L = 10000

    evi = G * ((nir_band - red_band) / (nir_band + C1 * red_band + C2 * green_band + L))
    return evi


def calculate_lsoa_stats(lsoa, rtree_idx, raster_bboxes):
    """Calculates statistics for a given LSOA.

    Args:
        lsoa (GeoSeries): LSOA geometry.
        rtree_idx (Index): R-tree index of raster files.
        raster_bboxes (dict): Dictionary of raster bounding boxes.

    Returns:
        dict: Statistics for the LSOA.
    """
    overlapping_rasters = find_overlapping_rasters_lsoa(lsoa, rtree_idx, raster_bboxes)
    all_ndvi_values, all_evi_values, all_raster_data = process_rasters(
        overlapping_rasters, lsoa
    )

    if all_ndvi_values:
        veg_fraction = np.sum(np.array(all_ndvi_values) > 0.2) / len(all_ndvi_values)
        cv_ndvi = np.nanstd(all_ndvi_values) / np.nanmean(all_ndvi_values)
        return compile_stats(
            lsoa,
            all_ndvi_values,
            all_evi_values,
            all_raster_data,
            veg_fraction,
            cv_ndvi,
        )
    else:
        return compile_empty_stats(lsoa)


def process_rasters(overlapping_rasters, lsoa):
    """Processes raster files to extract NDVI and EVI values.

    Args:
        overlapping_rasters (list): List of overlapping raster file paths.
        lsoa (GeoSeries): LSOA geometry.

    Returns:
        tuple: Lists of NDVI values, EVI values, and raster data.
    """
    all_ndvi_values = []
    all_evi_values = []
    all_raster_data = []

    for raster_file in overlapping_rasters:
        with rasterio.open(raster_file) as src:
            raster_data, _ = mask(src, [lsoa.geometry], crop=True)
            bands, height, width = raster_data.shape
            reshaped_data = raster_data.reshape(bands, height * width).T
            all_raster_data.extend(reshaped_data)

            zero_mask = np.all(reshaped_data != 0, axis=1)
            reshaped_nonzero = reshaped_data[zero_mask]
            reshaped_transposed = reshaped_nonzero.T

            green_band = reshaped_transposed[1]
            red_band = reshaped_transposed[2]
            nir_band = reshaped_transposed[3]

            ndvi = calculate_ndvi(red_band, nir_band)
            evi = calculate_evi(nir_band, red_band, green_band)
            if ndvi.size == 0:
                continue

            all_ndvi_values.extend(ndvi.flatten())
            all_evi_values.extend(evi.flatten())

    return all_ndvi_values, all_evi_values, all_raster_data


def compile_stats(
    lsoa, all_ndvi_values, all_evi_values, all_raster_data, veg_fraction, cv_ndvi
):
    """Compiles statistics for an LSOA with valid NDVI values.

    Args:
        lsoa (GeoSeries): LSOA geometry.
        all_ndvi_values (list): List of NDVI values.
        all_evi_values (list): List of EVI values.
        all_raster_data (list): List of raster data.
        veg_fraction (float): Vegetation fraction.
        cv_ndvi (float): Coefficient of variation for NDVI.

    Returns:
        dict: Compiled statistics.
    """
    return {
        "LSOA21CD": lsoa.LSOA21CD,
        "NDVI_MEAN": np.mean(all_ndvi_values),
        "NDVI_MEDIAN": np.median(all_ndvi_values),
        "NDVI_STD": np.nanstd(all_ndvi_values),
        "NDVI_MAX": np.max(all_ndvi_values),
        "NDVI_MIN": np.min(all_ndvi_values),
        "NDVI_CV": cv_ndvi,
        "VEG_FRAC": veg_fraction,
        "EVI_MEAN": np.mean(all_evi_values),
        "EVI_MEDIAN": np.median(all_evi_values),
        "EVI_STD": np.nanstd(all_evi_values),
        "EVI_MAX": np.max(all_evi_values),
        "EVI_MIN": np.min(all_evi_values),
        "PXL_COUNT": len(all_raster_data),
        "FINAL_PXL_COUNT": len(all_ndvi_values),
        "PCT_FILTERED": (len(all_ndvi_values) / len(all_raster_data)) * 100,
    }


def compile_empty_stats(lsoa):
    """Compiles empty statistics for an LSOA with no valid NDVI values.

    Args:
        lsoa (GeoSeries): LSOA geometry.

    Returns:
        dict: Compiled empty statistics.
    """
    return {
        "LSOA21CD": lsoa.LSOA21CD,
        "NDVI_MEAN": None,
        "NDVI_MEDIAN": None,
        "NDVI_MAX": None,
        "NDVI_MIN": None,
        "CV_NDVI": None,
        "VEG_FRAC": None,
        "EVI_MEAN": None,
        "EVI_MEDIAN": None,
        "EVI_STD": None,
        "EVI_MAX": None,
        "EVI_MIN": None,
        "PXL_COUNT": None,
        "FINAL_PXL_COUNT": None,
        "PCT_FILTERED": None,
    }


def main():
    raster_folder = Path("gisdata/13")
    raster_files = list(raster_folder.rglob("*.tif"))
    rtree_idx, raster_bboxes = create_rtree_index(raster_files)

    lsoa_boundaries = read_lsoa()

    results = [
        calculate_lsoa_stats(lsoa, rtree_idx, raster_bboxes)
        for lsoa in tqdm(lsoa_boundaries.itertuples())
    ]
    df = pd.DataFrame(results)

    df.to_parquet(Path("gisdata/ndvi.parquet"), index=False)


if __name__ == "__main__":
    main()
