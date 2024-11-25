import numpy as np
from osgeo import gdal
import rasterio
from rasterio.enums import Resampling
from scipy.ndimage import sobel, gaussian_filter


def dem_to_raster_pv_exclusion(dem_path, north_east_west_threshold, south_threshold, output_path, sigma=1):
    with rasterio.open(dem_path) as src:
        dem = src.read(1)
        profile = src.profile
        transform = src.transform
        nodata = src.nodata
        latitude = transform[5]

    # Calculate the approximate pixel size in meters for 3 arc-seconds
    pixel_size_m = 90 
    latitude_correction = np.cos(np.radians(latitude)) 
    pixel_size_x = pixel_size_m * latitude_correction
    pixel_size_y = pixel_size_m * latitude_correction

    # Create a mask for valid data (not NoData)
    if nodata is not None:
        valid_data_mask = dem != nodata
    else:
        valid_data_mask = np.ones_like(dem, dtype=bool)

    # Smooth the DEM using a Gaussian filter
    smoothed_dem = gaussian_filter(dem, sigma=sigma)

    # Calculate gradients using Sobel operators, incorporating pixel size
    gradient_x = sobel(smoothed_dem, axis=1) / (2 * pixel_size_x)  # X gradient (easting), adjusted for pixel size
    gradient_y = sobel(smoothed_dem, axis=0) / (2 * pixel_size_y)  # Y gradient (northing), adjusted for pixel size

    # Calculate the slope in degrees
    slope = np.degrees(np.arctan(np.sqrt(gradient_x**2 + gradient_y**2)))

    # Calculate the aspect in degrees
    aspect = np.degrees(np.arctan2(-gradient_y, gradient_x))
    aspect = np.where(aspect < 0, 360 + aspect, aspect)

    print(f"Slope min: {np.min(slope)}, max: {np.max(slope)}, mean: {np.mean(slope)}")
    print(f"Aspect min: {np.min(aspect)}, max: {np.max(aspect)}, mean: {np.mean(aspect)}")

    # Apply the slope thresholds based on aspect direction and valid data mask
    north_east_west_mask = (aspect >= 45) & (aspect <= 135) | (aspect >= 225) & (aspect <= 315)
    south_mask = (aspect > 135) & (aspect < 225)

    exclusion_mask = (
        ((north_east_west_mask) & (slope >= north_east_west_threshold)) |
        ((south_mask) & (slope >= south_threshold))
    ) & valid_data_mask

    # Calculate the percentage of points excluded
    total_points = np.sum(valid_data_mask)
    selected_points = np.sum(exclusion_mask)
    percentage_selected = (selected_points / total_points) * 100 if total_points > 0 else 0

    # Print the percentage of selected points
    print(f"Percentage of points excluded: {percentage_selected:.2f}%")

    # Update the profile to reflect the number of layers, dtype, and nodata value
    profile.update(
        dtype=rasterio.uint8,
        count=1,
        compress='deflate',
        nodata=0  # Set the NoData value to 0 for uint8
    )

    # Write the exclusion mask to a new raster file
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(exclusion_mask.astype(rasterio.uint8), 1)

    print(f"Output raster saved to {output_path}")


def dem_to_raster_wind_exlcusion(dem_path, threshold, output_path, sigma=1):
    # Read the DEM file
    with rasterio.open(dem_path) as src:
        dem = src.read(1)
        profile = src.profile
        transform = src.transform
        nodata = src.nodata

        # Extract latitude from transform (assuming it's a regular grid, so we use the top-left corner)
        latitude = transform[5]

    # Calculate the approximate pixel size in meters for 15 arc-seconds
    # It is 90m for 3 arc-seconds
    pixel_size_m = 90 # 450  # Approximate at the equator (in meters)
    latitude_correction = np.cos(np.radians(latitude))  # Correct for latitude
    pixel_size_x = pixel_size_m * latitude_correction
    pixel_size_y = pixel_size_m * latitude_correction

    # Create a mask for valid data (not NoData)
    if nodata is not None:
        valid_data_mask = dem != nodata
    else:
        valid_data_mask = np.ones_like(dem, dtype=bool)  # Assume all data is valid if no nodata is set

    # Smooth the DEM using a Gaussian filter
    smoothed_dem = gaussian_filter(dem, sigma=sigma)

    # Calculate gradients using Sobel operators, incorporating pixel size
    gradient_x = sobel(smoothed_dem, axis=1) / (2 * pixel_size_x)  # X gradient (easting), adjusted for pixel size
    gradient_y = sobel(smoothed_dem, axis=0) / (2 * pixel_size_y)  # Y gradient (northing), adjusted for pixel size

    # Calculate the slope in degrees
    slope = np.degrees(np.arctan(np.sqrt(gradient_x**2 + gradient_y**2)))

    # Debugging: Print slope statistics
    print(f"Slope min: {np.min(slope)}, max: {np.max(slope)}, mean: {np.mean(slope)}")

    # Apply the slope threshold and the valid data mask
    mask = (slope >= threshold) & valid_data_mask

    # Calculate the percentage of points with value 1
    total_points = np.sum(valid_data_mask)
    selected_points = np.sum(mask)
    percentage_selected = (selected_points / total_points) * 100 if total_points > 0 else 0

    # Print the percentage of selected points
    print(f"Percentage of points with slope >= {threshold} degrees: {percentage_selected:.2f}%")

    # Debugging: Print mask statistics
    print(f"Threshold: {threshold}, Mask sum: {selected_points}, Mask mean: {np.mean(mask)}")

    # Update the profile to reflect the number of layers, dtype, and nodata value
    profile.update(
        dtype=rasterio.uint8,
        count=1,
        compress='deflate',
        nodata=0  # Set the NoData value to 0 for uint8
    )

    # Write the mask to a new raster file
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(mask.astype(rasterio.uint8), 1)

    print(f"Output raster saved to {output_path}")