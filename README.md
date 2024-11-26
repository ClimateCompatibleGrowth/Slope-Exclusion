
# **SlopeEXCL - Slope-based Exclusion**

SlopeEXCL generates `.tif` files to exclude areas unsuitable for wind and solar energy placements based on slope thresholds.

- **Wind Turbines**: Applies a uniform slope threshold across all aspects.
- **Solar Panels**: Allows distinct thresholds for south-facing and north-east-west-facing slopes.

---

## **About SlopeEXCL**

SlopeEXCL uses slope-based exclusions to identify suitable areas for wind turbines and PV panels. This methodology integrates slope gradients and aspect orientation to enhance spatial assessments for renewable energy planning, which is crucial in mountainous regions. It can be added as an additional exclusion in tools like [GLAES](https://github.com/FZJ-IEK3-VSA/glaes). For example:

```python
ec.excludeRasterType(os.path.join(data_path, 'Laos_slope_excluded_pv.tif'), value=1, prewarp=True)
current_time = time.time() - start_time
```

---

## **Getting Started**

This repository requires a 3-arc-second resolution Digital Elevation Model (DEM) from [HydroSHEDS](https://www.hydrosheds.org/hydrosheds-core-downloads). Place the downloaded DEM in the `data` folder.

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

   To use the **default values**, omit the optional arguments::
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
