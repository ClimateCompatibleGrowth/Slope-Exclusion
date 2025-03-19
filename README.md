
# **SlopeEXCL - Slope-based Exclusion**

SlopeEXCL generates `.tif` files to exclude areas unsuitable for wind and solar energy placements based on slope thresholds.

- **Wind Turbines**: Applies a uniform slope threshold across all aspects.
- **Solar Panels**: Allows distinct thresholds for south-facing and north-east-west-facing slopes.

The generated .tif files can be integrated as exclusion layers in renewable energy planning tools like [GLAES](https://github.com/FZJ-IEK3-VSA/glaes). For example, in Glaes:

```python
ec.excludeRasterType(os.path.join(data_path, 'exclude_slope_pv.tif'), value=1, prewarp=True)
```

---

## **Getting Started**

This repository requires a 3-arc-second resolution Digital Elevation Model (DEM) from [HydroSHEDS](https://www.hydrosheds.org/hydrosheds-core-downloads). Place the downloaded DEM in the `data` folder.

> **Note:** DEMs for entire continents are often too large to process efficiently.  
> To avoid excessive file sizes, clip the DEM to your area of interest (e.g., a country).  
> You can download country boundary shapefiles from [Natural Earth](https://www.naturalearthdata.com/downloads/50m-cultural-vectors/50m-admin-0-countries-2/).  

### **How to Clip the DEM**

Before running the main exclusion script, clip the DEM to your country using the provided script: 

```bash
python clip_raster_to_boundary.py --raster data/dem.tif --boundary data/country_boundary.gpkg --output data/clipped_dem.tif
```


---

## **How to Run**

1. Set up Python and install the required dependencies (e.g., `rasterio`, `numpy`, `scipy`).
2. Execute the script from the terminal:
   ```bash
   python exclude_slope.py --type [solar|wind|both] --solar-nea [value] --solar-s [value] --wind-thresh [value] --output [filename]
   ```
   Example with custom values:
   ```bash
   python exclude_slope.py --type both --solar-nea 6.28 --solar-s 33 --wind-thresh 8.53 --output exclusion.tif
   ```

   To use the **default values**, omit the optional arguments:
   ```bash
   python exclude_slope.py --type both
   ```

   The script will then use:
   - `--solar-nea`: `6.28` (degrees)
   - `--solar-s`: `33` (degrees)
   - `--wind-thresh`: `8.53` (degrees)
   - `--sigma`: `1`
   - `--output`: `exclusion.tif`

---

## **Repository Structure**

- **`exclude_slope.py`**: Python script for generating exclusion files.
- **`exclude_slope.ipynb`**: Jupyter notebook for running the same code interactively.
- **`data`**: Directory for raw input file `dem.tif`.
- **`output`**: Directory where generated `.tif` files will be saved.
