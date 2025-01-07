import glob

import dask
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from dask.delayed import delayed
from dask.diagnostics.progress import ProgressBar
from rasterio.mask import mask
from rtree import index
from scipy.ndimage import median_filter
from scipy.spatial.distance import mahalanobis
from shapely.geometry import box
from tqdm import tqdm


def read_lsoa():
    nidz = (
        gpd.read_file("./gisdata/geography-dz2021-esri-shapefile.zip")[
            ["DZ2021_cd", "geometry"]
        ]
        .to_crs(4326)
        .rename(columns={"DZ2021_cd": "LSOA21CD"})
    )
    sgdz = (
        gpd.read_file(
            "./gisdata/SG_DataZoneBdry_2022.zip", layer="SG_DataZoneBdry_2022_EoR"
        )[["DZCode", "geometry"]]
        .to_crs(4326)
        .rename(columns={"DZCode": "LSOA21CD"})
    )
    lsoa_boundaries = gpd.read_file(
        "./gisdata/gov/LSOA2021/LSOA_2021_EW_BFC_V8.shp"
    ).to_crs(4326)[["LSOA21CD", "geometry"]]
    lsoa_boundaries = pd.concat([lsoa_boundaries, sgdz, nidz])
    return lsoa_boundaries


def create_rtree_index(raster_files):
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
    lsoa_bbox = lsoa.geometry.bounds
    overlapping_rasters = list(rtree_idx.intersection(lsoa_bbox))

    valid_rasters = []
    for i in overlapping_rasters:
        raster_bbox, raster_file = raster_bboxes[i]
        if lsoa.geometry.intersects(raster_bbox):
            valid_rasters.append(raster_file)
    return valid_rasters


def detect_outliers(reshaped_nonzero):
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
    ndvi = np.divide(
        nir_band - red_band,
        nir_band + red_band,
        where=(nir_band + red_band) != 0,
    )
    # values above 1 are definite outliers, below 0.1 are likely snow/water maybe some cloud
    ndvi[(ndvi > 1) | (ndvi < 0.1)] = 0
    # remove nans (dividing by zero) and values of exactly zero
    ndvi = ndvi[(~np.isnan(ndvi)) & (ndvi != 0)]
    # smooth NDVI values to reduce outliers
    ndvi = median_filter(ndvi, size=10)
    return ndvi


def calculate_evi(nir_band, red_band, green_band):
    # Constants for EVI calculation
    G = 2.5  # Gain factor
    C1 = 6  # Coefficient for Red band
    C2 = 7.5  # Coefficient for Green band
    L = 10000  # Canopy background factor

    # Calculate EVI using the three bands
    evi = G * ((nir_band - red_band) / (nir_band + C1 * red_band + C2 * green_band + L))

    return evi


def calculate_lsoa_stats(lsoa, rtree_idx, raster_bboxes):
    overlapping_rasters = find_overlapping_rasters_lsoa(lsoa, rtree_idx, raster_bboxes)
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

            # if (reshaped_nonzero.size == 0) | (reshaped_nonzero.shape[0] < 10):
            #     print("reshaped size", reshaped_nonzero.shape[0])
            #     continue

            # reshaped_filtered = detect_outliers(reshaped_nonzero)
            # reshaped_transposed = reshaped_filtered.T
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

    veg_fraction = np.sum(np.array(all_ndvi_values) > 0.2) / len(all_ndvi_values)
    cv_ndvi = np.nanstd(all_ndvi_values) / np.nanmean(all_ndvi_values)

    if all_ndvi_values:
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
    else:
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
    raster_folder = "gisdata/13/**/*.tif"
    raster_files = glob.glob(raster_folder, recursive=True)
    rtree_idx, raster_bboxes = create_rtree_index(raster_files)

    lsoa_boundaries = read_lsoa()

    tasks = [
        delayed(calculate_lsoa_stats)(lsoa, rtree_idx, raster_bboxes)
        for lsoa in lsoa_boundaries.copy().itertuples()
    ]

    with ProgressBar():
        results = dask.compute(*tasks)
    # results = [
    #     calculate_lsoa_stats(lsoa, rtree_idx, raster_bboxes)
    #     for lsoa in tqdm(lsoa_boundaries.itertuples())
    # ]
    df = pd.DataFrame(results)

    df.to_parquet("./gisdata/ndvi.parquet", index=False)


if __name__ == "__main__":
    main()
    df = pd.read_parquet("./gisdata/ndvi.parquet")
    ahah_gs = pd.read_csv("./gisdata/AHAH_V4.csv")[["LSOA21CD", "ah4gpas"]]
    ahah_gs.merge(df, on="LSOA21CD")[["NDVI_MEAN", "ah4gpas"]].corr()
    merged = lsoa_boundaries.merge(df)

    merged.plot("NDVI_MEAN", scheme="quantiles", legend=True)
    plt.show()
