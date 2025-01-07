import matplotlib.pyplot as plt
import pandas as pd

from main import read_lsoa

if __name__ == "__main__":
    lsoa_boundaries = read_lsoa()
    df = pd.read_parquet("./gisdata/ndvi.parquet")

    fig, ax = plt.subplots()
    ax.axis("off")
    lsoa_boundaries.merge(df, on="LSOA21CD").plot("NDVI_MEAN", ax=ax)
    plt.show()
