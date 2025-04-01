
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
Download of files is indicated by 📋.

Code execution is indicated by 👉.

Set up Python and install the required dependencies (`rasterio`, `numpy`, `scipy`, `os`, `argparse`).

📋 This repository requires a 3-arc-second resolution Digital Elevation Model (DEM) from [HydroSHEDS](https://www.hydrosheds.org/hydrosheds-core-downloads). Place the downloaded DEM in the `data` folder and rename to "_fulL_dem_".

> **Note:** DEMs for entire continents are often too large to process efficiently. Clip them to your area of interest (see next step).

### **How to Clip the DEM**
📋 To avoid excessive file sizes, clip the DEM to your area of interest (e.g., a country). You can download country boundary GeoJson files from [opendatasoft](https://public.opendatasoft.com/explore/dataset/world-administrative-boundaries/export/) based on their ISO3 code. Place the .geojson file in the `data` folder and rename to "_country_boundary_".

👉 Clip the DEM to your country using the provided script: 

```bash
python clip_raster_to_boundary.py --raster data/full_dem.tif --boundary data/country_boundary.geojson --output data/dem.tif
```


---

## **How to Run**
👉 Execute the script from the terminal:
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
