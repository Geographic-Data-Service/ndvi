import matplotlib.pyplot as plt
import pandas as pd

from src.main import read_lsoa


def plot_veg_frac(lsoa_boundaries, df):
    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))  # Adjust the figure size
    ax.axis("off")  # Remove axes for a clean look

    # Add title
    ax.set_title("Vegetation Fraction by LSOA", fontsize=16, fontweight="bold")

    # Plot the data
    cmap = "YlGn"  # Use a perceptually uniform colormap
    merged_data = lsoa_boundaries.merge(df, on="LSOA21CD")
    plot = merged_data.plot(
        "VEG_FRAC",
        ax=ax,
        legend=True,
        cmap=cmap,
        legend_kwds={
            "label": "Vegetation Fraction (%)",
            "orientation": "horizontal",
            "shrink": 0.6,
            "pad": 0.02,
        },
    )

    # Save the figure with a higher resolution
    plt.savefig("./img/veg_frac.png", dpi=300, bbox_inches="tight")

    # Display the plot
    plt.show()


if __name__ == "__main__":
    lsoa_boundaries = read_lsoa()
    df = pd.read_parquet("./gisdata/ndvi.parquet")

    plot_veg_frac(lsoa_boundaries, df)
