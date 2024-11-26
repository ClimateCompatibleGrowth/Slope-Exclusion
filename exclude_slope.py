import os
import argparse
import numpy as np
import rasterio
from scipy.ndimage import sobel, gaussian_filter


def process_dem(dem_path, sigma=1):
    """
    Process a DEM to calculate slope and aspect.

    This function reads a Digital Elevation Model (DEM) from a file, applies a Gaussian filter for smoothing, 
    and calculates the slope and aspect of the terrain. The slope is returned in degrees, and the aspect is 
    returned as compass directions in degrees.

    Args:
        dem_path (str): Path to the input DEM file (e.g., a GeoTIFF).
        sigma (float): Standard deviation for the Gaussian filter applied to smooth the DEM.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: The slope array in degrees.
            - numpy.ndarray: The aspect array in degrees (0-360, with 0 = north).
            - dict: The raster profile (metadata) of the input DEM.
            - int or float or None: The NoData value of the DEM, if defined.
    """
    with rasterio.open(dem_path) as src:
        dem = src.read(1)
        profile = src.profile
        transform = src.transform
        nodata = src.nodata

        latitude = transform[5]
        pixel_size_m = 90  # Approximate pixel size for 3 arc-seconds
        latitude_correction = np.cos(np.radians(latitude))
        pixel_size_x = pixel_size_m * latitude_correction
        pixel_size_y = pixel_size_m * latitude_correction

        smoothed_dem = gaussian_filter(dem, sigma=sigma)
        gradient_x = sobel(smoothed_dem, axis=1) / (2 * pixel_size_x)
        gradient_y = sobel(smoothed_dem, axis=0) / (2 * pixel_size_y)

        slope = np.degrees(np.arctan(np.sqrt(gradient_x**2 + gradient_y**2)))
        aspect = np.degrees(np.arctan2(-gradient_y, gradient_x))
        aspect = np.where(aspect < 0, 360 + aspect, aspect)

    return slope, aspect, profile, nodata


def calculate_percentage_excluded(mask):
    """
    Calculate the percentage of excluded points.

    This function computes the percentage of points in a given mask that are marked as excluded (True).

    Args:
        mask (numpy.ndarray): A boolean array where True indicates excluded points.

    Returns:
        float: The percentage of excluded points relative to the total number of points.
    """
    total_points = np.prod(mask.shape)  # Total number of points
    excluded_points = np.sum(mask)     # Points where exclusion is True
    percentage_excluded = (excluded_points / total_points) * 100
    return percentage_excluded


def create_combined_exclusion(dem_path, output_name, solar_nea, solar_s, wind_thresh, sigma=1):
    """
    Create and save a combined exclusion mask for solar and wind.

    This function generates exclusion masks for solar and wind based on terrain slope and aspect, 
    combines the masks, and saves the result to a raster file.

    Args:
        dem_path (str): Path to the input DEM file (e.g., a GeoTIFF).
        output_name (str): Name of the output raster file.
        solar_nea (float): Slope threshold for north, east, and west-facing slopes for solar PV.
        solar_s (float): Slope threshold for south-facing slopes for solar PV.
        wind_thresh (float): Slope threshold for wind exclusion.
        sigma (float): Standard deviation for the Gaussian filter applied to smooth the DEM.

    Returns:
        None
    """
    slope, aspect, profile, nodata = process_dem(dem_path, sigma)

    north_east_west_mask = (aspect >= 45) & (aspect <= 135) | (aspect >= 225) & (aspect <= 315)
    south_mask = (aspect > 135) & (aspect < 225)
    solar_exclusion = (
        ((north_east_west_mask) & (slope >= solar_nea)) |
        ((south_mask) & (slope >= solar_s))
    )
    wind_exclusion = slope >= wind_thresh
    combined_exclusion = solar_exclusion | wind_exclusion

    solar_percentage = calculate_percentage_excluded(solar_exclusion)
    wind_percentage = calculate_percentage_excluded(wind_exclusion)
    combined_percentage = calculate_percentage_excluded(combined_exclusion)

    print(f"Percentage of points excluded for solar: {solar_percentage:.2f}%")
    print(f"Percentage of points excluded for wind: {wind_percentage:.2f}%")
    print(f"Percentage of points excluded (combined): {combined_percentage:.2f}%")

    save_exclusion_raster(combined_exclusion, profile, output_name)


def create_solar_exclusion(dem_path, output_name, solar_nea, solar_s, sigma=1):
    """
    Create and save a solar exclusion mask.

    This function generates an exclusion mask for solar PV development based on terrain slope and aspect 
    and saves the result to a raster file.

    Args:
        dem_path (str): Path to the input DEM file (e.g., a GeoTIFF).
        output_name (str): Name of the output raster file.
        solar_nea (float): Slope threshold for north, east, and west-facing slopes for solar PV.
        solar_s (float): Slope threshold for south-facing slopes for solar PV.
        sigma (float): Standard deviation for the Gaussian filter applied to smooth the DEM.

    Returns:
        None
    """
    slope, aspect, profile, nodata = process_dem(dem_path, sigma)

    north_east_west_mask = (aspect >= 45) & (aspect <= 135) | (aspect >= 225) & (aspect <= 315)
    south_mask = (aspect > 135) & (aspect < 225)
    solar_exclusion = (
        ((north_east_west_mask) & (slope >= solar_nea)) |
        ((south_mask) & (slope >= solar_s))
    )

    solar_percentage = calculate_percentage_excluded(solar_exclusion)
    print(f"Percentage of points excluded for solar: {solar_percentage:.2f}%")

    save_exclusion_raster(solar_exclusion, profile, output_name)


def create_wind_exclusion(dem_path, output_name, wind_thresh, sigma=1):
    """
    Create and save a wind exclusion mask.

    This function generates an exclusion mask for wind development based on terrain slope 
    and saves the result to a raster file.

    Args:
        dem_path (str): Path to the input DEM file (e.g., a GeoTIFF).
        output_name (str): Name of the output raster file.
        wind_thresh (float): Slope threshold for wind exclusion.
        sigma (float): Standard deviation for the Gaussian filter applied to smooth the DEM.

    Returns:
        None
    """
    slope, _, profile, nodata = process_dem(dem_path, sigma)

    wind_exclusion = slope >= wind_thresh

    wind_percentage = calculate_percentage_excluded(wind_exclusion)
    print(f"Percentage of points excluded for wind: {wind_percentage:.2f}%")

    save_exclusion_raster(wind_exclusion, profile, output_name)


def save_exclusion_raster(mask, profile, output_name):
    """
    Save an exclusion mask to a raster file.

    This function writes a boolean mask to a GeoTIFF file, ensuring the mask is saved 
    as a single-layer raster with the appropriate metadata.

    Args:
        mask (numpy.ndarray): Boolean array representing the exclusion mask.
        profile (dict): Raster profile for the output file (e.g., spatial resolution, projection).
        output_name (str): Name of the output raster file.

    Returns:
        None
    """
    output_folder = "output"
    os.makedirs(output_folder, exist_ok=True)

    output_path = os.path.join(output_folder, output_name)

    profile.update(
        dtype=rasterio.uint8,
        count=1,
        compress='deflate',
        nodata=0
    )
    
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(mask.astype(rasterio.uint8), 1)

    print(f"Exclusion raster saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate exclusion rasters for solar and wind development.")
    parser.add_argument(
        "--type",
        required=True,
        choices=["solar", "wind", "both"],
        help="Type of exclusion to generate ('solar', 'wind', or 'both')."
    )
    parser.add_argument(
        "--solar-nea",
        type=float,
        default=6.28,  # Default value in degrees for north-east-west slopes
        help="Slope threshold for north, east, and west-facing slopes for solar PV (degrees)."
    )
    parser.add_argument(
        "--solar-s",
        type=float,
        default=33,  # Default value in degrees for south slopes
        help="Slope threshold for south-facing slopes for solar PV (degrees)."
    )
    parser.add_argument(
        "--wind-thresh",
        type=float,
        default=8.53,  # Default value in degrees for wind exclusion
        help="Slope threshold for wind exclusion (degrees)."
    )
    parser.add_argument(
        "--sigma",
        type=float,
        default=1,
        help="Smoothing factor for the Gaussian filter."
    )
    parser.add_argument(
        "--output",
        type=str,
        default="exclusion.tif",
        help="Filename for the exclusion raster (saved in the 'output' folder)."
    )

    args = parser.parse_args()

    dem_path = "data/dem.tif"

    if args.type == "both":
        create_combined_exclusion(
            dem_path=dem_path,
            output_name=args.output,
            solar_nea=args.solar_nea,
            solar_s=args.solar_s,
            wind_thresh=args.wind_thresh,
            sigma=args.sigma
        )
    elif args.type == "solar":
        create_solar_exclusion(
            dem_path=dem_path,
            output_name=args.output,
            solar_nea=args.solar_nea,
            solar_s=args.solar_s,
            sigma=args.sigma
        )
    elif args.type == "wind":
        create_wind