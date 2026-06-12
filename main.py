#import libraries
import ee
import geemap
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import rasterio

#Configure the GEE API

TRAIN = False #Takes data from the geemap and saves it.

OUTPUT_DIR = 'outputs'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Initialize the Earth Engine API
ee.Initialize(project='inbound-lattice-403009')
ROI = ee.Geometry.Rectangle([116.0, 0.5, 118.5, 2.5])
YEAR = 2023

print(f"ROI: {ROI.getInfo()['coordinates']}")
print(f"Year: {YEAR}")

#Functions
def get_sentinal2(year, roi):
    '''
    We need all bands to compute multiple indices.
    This function loads Sentinal-2 data that is cloud-masked for the specific ROI and year.
    Returns image for all spectral bands retained for multi-index compuatation.
    '''
    start = f"{year}-01-01"
    end = f"{year}-12-31"

    #Make it cloud-free 
    def mask_s2_clouds(image):
        qa =image.select('QA60') #QA60 is for sentinal-2 image
        cloud_bit_mask = 1 << 10 #bit 10 is clouds 
        cirrus_bit_mask = 1 << 11 # bit 11 is cirrus clouds
        
        #Both flags should be set to zero, indicating clear conditions.
        mask = (
            qa.bitwiseAnd(cloud_bit_mask).eq(0)
            .And(qa.bitwiseAnd(cirrus_bit_mask).eq(0)) 
        )
        return image.updateMask(mask).divide(10000)
    
    #collect the images
    collection = (
        ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(roi)
        .filterDate(start, end)
        #Pre-filter to get less cloudy granules.
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30))
        .map(mask_s2_clouds)
    )

    count = collection.size().getInfo()
    print(f"Found {count} images for {year}")
    if count == 0:
        raise ValueError(f"No images found for the {year}. Try relaxing cloud filter.")
    '''
    In remote sensing specifically:

    Median → compositing (what I'm doing here) — removes cloud contamination
    Mean → calculating anomalies, e.g. "how much warmer is this pixel than the 10-year average temperature?" Here you want all values to contribute
    Mean → when you've already masked all bad pixels perfectly and want maximum use of clean data
    '''
    return collection.median().clip(roi) #single composite image representing the whole year

def compute_indices(image):
    '''
    compte NDVI, EVI, NDWI
    NDVI = (NIR - R)/(NIR + R)
    Bands used = B8 & B4
    High values indicates dense vegetation/forest
    low values indicates bare soil/water/urban areas
    '''   
    ndvi = image.normalizedDifference(["B8", "B4"]).rename('NDVI')
    evi = image.expression(
        "2.5 * ((n - r) / (n + 6 * r - 7.5 * b + 1))",
        {
            "n": image.select("B8"),
            "r": image.select("B4"),
            "b": image.select("B2")
        }
    ).rename("EVI")

    ndwi = image.normalizedDifference(["B3", "B8"]).rename("NDWI")

    return ee.Image.cat([ndvi, evi, ndwi])

def classify_vegetation(ndvi, evi, ndwi):
    '''
    Multi-index threshold classification into 5 land cover classes.
    Using multiple indices gives more accurate classification:

    Class 0 - Water:
      NDWI > 0.1 reliably identifies water bodies
      (NDWI is designed specifically for water detection)

    Class 1 - Bare soil / Urban:
      Low NDVI (< 0.15) and not water
      Little to no vegetation signal

    Class 2 - Sparse / Degraded vegetation:
      NDVI 0.15-0.35, low EVI
      Grass, scrub, recently cleared land

    Class 3 - Agriculture / Plantation:
      NDVI 0.35-0.6, moderate EVI (0.1-0.35)
      Crops, oil palm, young plantations
      EVI helps distinguish from natural forest

    Class 4 - Dense Forest:
      High NDVI (> 0.6) and high EVI (> 0.35)
      Mature tropical forest canopy
    '''
    water = ndvi.gt(0.1).multiply(0)
    bare = ndvi.lte(0.1).And(ndvi.lt(0.15)).multiply(1)
    sparse = ndvi.gte(0.15).And(ndvi.lt(0.35)).multiply(2)
    agri = ndvi.gte(0.36).And(ndvi.lt(0.6)).And(evi.lt(0.35)).multiply(3)
    forest = ndvi.gte(0.6).And(evi.gte(0.35)).multiply(4)

    #Combine so that each pixel falls into one class
    classified = water.add(bare).add(sparse).add(agri).add(forest).rename("class")
    return classified

# Main pipeline
if TRAIN:
    print(f"Get sentinal-2 images for the year {YEAR}....")
    image = get_sentinal2(YEAR, ROI)

    print(f"Compute the indices (NDVI, EVI, NDWI)...")
    indices = compute_indices(image)

    ndvi = indices.select("NDVI")
    evi = indices.select("EVI")
    ndwi = indices.select("NDWI")

    print(f"Classifying vegetation/landcover....")
    classified = classify_vegetation(ndvi, evi, ndwi)

    print("EXPORTING GeoTIFFS...")
    geemap.ee_export_image(
        ndvi, filename=f"{OUTPUT_DIR}/ndvi_{YEAR}.tif",
        scale=100, region=ROI, file_per_band=False
    )
    geemap.ee_export_image(
        evi, filename=f"{OUTPUT_DIR}/evi_{YEAR}.tif",
        scale=200, region=ROI, file_per_band=False
    )
    geemap.ee_export_image(
        ndwi, filename=f"{OUTPUT_DIR}/ndwi_{YEAR}.tif",
        scale=100, region=ROI, file_per_band=False
    )
    geemap.ee_export_image(
        classified, filename=f"{OUTPUT_DIR}/classified_{YEAR}.tif",
        scale=100, region=ROI, file_per_band=False
    )

    print("Export complete. Set TRAIN=False for visualizing.")

else:
    print("Loading the saved outputs for visualization...")
    files = {
        "ndvi": f"{OUTPUT_DIR}/ndvi_{YEAR}.tif",
        "evi": f"{OUTPUT_DIR}/evi_{YEAR}.tif",
        "ndwi": f"{OUTPUT_DIR}/ndwi_{YEAR}.tif",
        "classified": f"{OUTPUT_DIR}/classified_{YEAR}.tif"
    }

    #Check if any files are missing
    missing = [k for k, v in files.items() if not os.path.exists(v)]
    if missing:
        print("Missing files: {missing}")
        print("First run with TRAIN=True")
        exit()

    def read_tif(path):
        with rasterio.open(path) as src:
            return src.read(1) #read the first band of the GeoTIFF file
        
    ndvi = read_tif(files["ndvi"])
    evi = read_tif(files["evi"])
    ndwi = read_tif(files["ndwi"])
    classified = read_tif(files["classified"])

    #Plot 
    class_cmap = mcolors.ListedColormap(["#2166ac", "#d9b365", "#f4e04d", "#74c476", "#005a32"])
    class_norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], class_cmap.N)
    class_labels = ["Water", "Bare/Urban", "Sparse Vegetation", "Agriculture", "Forest"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle(f"Crop & Vegetation Mapping — East Kalimantan, Borneo\n"
                 f"{YEAR} | Sentinel-2 | NDVI + EVI + NDWI",
                 fontsize=13, fontweight="bold")
    
    #NDVI
    im0 = axes[0, 0].imshow(ndvi, cmap="RdYlGn", vmin=-0.2, vmax=0.8)
    axes[0, 0].set_title(f"NDVI {YEAR}")
    plt.colorbar(im0, ax=axes[0, 0], fraction=0.046)

    # EVI
    im1 = axes[0, 1].imshow(evi, cmap="RdYlGn", vmin=-0.2, vmax=0.6)
    axes[0, 1].set_title(f"EVI {YEAR}")
    plt.colorbar(im1, ax=axes[0, 1], fraction=0.046)

    # NDWI
    im2 = axes[1, 0].imshow(ndwi, cmap="RdBu", vmin=-0.5, vmax=0.5)
    axes[1, 0].set_title(f"NDWI {YEAR}")
    plt.colorbar(im2, ax=axes[1, 0], fraction=0.046)

    # Classification
    im3 = axes[1, 1].imshow(classified, cmap=class_cmap, norm=class_norm)
    axes[1, 1].set_title(f"Land Cover Classification {YEAR}")
    cbar = plt.colorbar(im3, ax=axes[1, 1], fraction=0.046,
                        ticks=[0, 1, 2, 3, 4])
    cbar.set_ticklabels(class_labels)
    for ax in axes.flat:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/vegetation_map_{YEAR}.png", dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {OUTPUT_DIR}/vegetation_map_{YEAR}.png")

    # STATS 
    print(f"\nLand Cover Statistics {YEAR} ")
    total = classified.size
    for i, label in enumerate(class_labels):
        pct = np.sum(classified == i) / total * 100
        print(f"{label:20s}  {pct:5.1f}%")


        
    
        
    





