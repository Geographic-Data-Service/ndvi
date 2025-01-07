# Sentinel Data Processing

This project processes Sentinel data to calculate vegetation indices and compile statistics for Local Super Output Areas (LSOAs). The output is saved as a Parquet file containing NDVI and EVI statistics for each LSOA.

## Prerequisites

- Python 3.8 or higher
- Required Python packages: geopandas, numpy, pandas, rasterio, rtree, scipy, shapely, tqdm, matplotlib

## Setup

1. Clone the repository and navigate to the project directory.

2. Install the required Python packages using `uv` (optional):

   ```bash
   uv sync
   ```

   Alternatively, you can use pip:

   ```bash
   pip install -r requirements.txt
   ```

3. Ensure the Sentinel raster data and shapefiles are placed in the appropriate directories as specified in `utils.py`.

## Data Files

- **Raster Data**: Place all Sentinel raster files in the `gisdata/13` directory.
- **Shapefiles**: Ensure the following shapefiles are available:
  - Northern Ireland Data Zones: `gisdata/geography-dz2021-esri-shapefile.zip`
  - Scottish Data Zones: `gisdata/SG_DataZoneBdry_2022.zip`
  - LSOA Boundaries: `gisdata/gov/LSOA2021/LSOA_2021_EW_BFC_V8.shp`

## Processing Steps

1. **Read LSOA Boundaries**: The script reads and combines LSOA boundary data from Northern Ireland, Scotland, and England/Wales.

2. **Create R-tree Index**: An R-tree index is created for efficient spatial querying of raster files.

3. **Find Overlapping Rasters**: For each LSOA, the script identifies raster files that overlap with the LSOA geometry.

4. **Process Rasters**: The script processes each overlapping raster to extract NDVI and EVI values:
   - NDVI is calculated using the red and NIR bands.
   - EVI is calculated using the NIR, red, and green bands.
   - Outliers are detected and removed using Mahalanobis distance.

5. **Compile Statistics**: For each LSOA, the script compiles statistics such as mean, median, standard deviation, and vegetation fraction for NDVI and EVI.

6. **Save Results**: The compiled statistics are saved to a Parquet file (`gisdata/ndvi.parquet`).

## Running the Script

To run the script and generate the output file, execute:

```bash
python main.py
```

## Visualizing Results

After running the script, you can visualize the results using the following command:

```bash
python -c "import matplotlib.pyplot as plt; plt.show()"
```

This will display a plot of the mean NDVI values for each LSOA.

## License

This project is licensed under the MIT License.
