import argparse
import rasterio
import geopandas as gpd
from rasterio.mask import mask
from shapely.geometry import mapping

def clip_raster_to_boundary(raster_path, boundary_path, output_path):
    """
    Clips a raster to the boundary defined in a vector file (GeoPackage or Shapefile).

    Parameters:
    raster_path (str): Path to the input raster file (.tif).
    boundary_path (str): Path to the vector file (.gpkg, .shp) containing the boundary.
    output_path (str): Path to save the clipped raster file.
    """
    # Load the boundary file
    boundary = gpd.read_file(boundary_path)
    
    # Ensure the GeoDataFrame is in the same coordinate reference system (CRS) as the raster
    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
        boundary = boundary.to_crs(raster_crs)

        # Get the geometry as a GeoJSON-like dict
        geometry = [mapping(boundary.geometry.union_all())]

        # Mask the raster with the boundary
        out_image, out_transform = mask(src, geometry, crop=True)
        out_meta = src.meta.copy()

    # Update metadata
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform
    })

    # Save the clipped raster
    with rasterio.open(output_path, "w", **out_meta) as dest:
        dest.write(out_image)
    
    print(f"Clipping completed. The new file is saved as '{output_path}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clip a raster to a boundary.")
    parser.add_argument("--raster", required=True, help="Path to the input raster file (.tif)")
    parser.add_argument("--boundary", required=True, help="Path to the vector file (.gpkg or .shp)")
    parser.add_argument("--output", required=True, help="Path to save the clipped raster file")
    args = parser.parse_args()
    
    clip_raster_to_boundary(args.raster, args.boundary, args.output)