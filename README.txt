# IWD Plugin for QGIS: Illegal Waste Detection via NDVI

## Overview

This QGIS plugin is designed to **detect illegal waste dumps** using **NDVI (Normalized Difference Vegetation Index)** maps generated from **multispectral drone imagery**. It works by identifying vegetation anomalies in NDVI maps and highlighting regions that likely represent illegal dumping sites.

The plugin is open-source and compatible with **QGIS 3.4+**.

---

## 🔧 Installation Instructions

### 1. Unpack the Plugin

- Download and unzip the file `Illegal waste dumping.zip`.

### 2. Copy to QGIS Plugins Directory

Copy the entire unzipped folder (e.g., `Illegal waste dumping` or similar) into your **QGIS plugin directory**.

**Typical locations:**
- **Windows:** `C:\Users\<YourUsername>\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\`
- **Linux:** `/home/<yourusername>/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`
- **MacOS:** `/Users/<yourusername>/Library/Application Support/QGIS/QGIS3/profiles/default/python/plugins/`

> Restart QGIS after copying.

### 3. Activate the Plugin

- Open **QGIS**.
- Go to **Plugins > Manage and Install Plugins...**
- Find **IWD** or **Illegal Waste Detection** in the list.
- Click **Install** or **Enable**.

---

## ▶️ How to Use the Plugin

1. **Load NDVI Raster**
   - Load your pre-processed NDVI orthomosaic (.tif) into QGIS.
   - This raster must already represent vegetation index (NDVI), where low values indicate possible disturbed areas (e.g., bare soil, waste).

2. **Open Plugin**
   - Navigate to **Plugins > IWD** and click **Detect Dumps**.

3. **Set Parameters**
   - **Layer Selection:** Select your NDVI raster from the list or browse to add a new one.
   - **Use Current Extent:** Optionally restrict analysis to your current map view.
   - **NDVI Threshold:** Choose the NDVI threshold below which an area is considered suspicious.
   - **Minimum Area (m²):** Exclude small polygons (e.g., < 5 m²).
   - **Add Results to Map:** Toggle to automatically add output as a layer.

4. **Run Detection**
   - Click **Run Detection** to process.
   - After completion, a message shows how many potential dump sites were found and what NDVI threshold was used.

5. **Review Dumps**
   - A simple GUI allows reviewing detected polygons one-by-one.
   - Confirm if the object is a dump (`Yes/No`). 


---

## 📁 Input Data Requirements

- **NDVI raster (.tif)** – GeoTIFF format preferred, georeferenced. Generated from a multispectral drone image (e.g., using RED + NIR bands).
- **Coordinate Reference System (CRS):** Preferably projected in meters (e.g., EPSG:3912 or UTM) for accurate area calculation.

---

## 📄 Output

- **GeoJSON** file with detected dump polygons.

---

## 🛠 Notes

- The plugin does not compute NDVI from raw imagery. It expects the user to **pre-calculate NDVI** before use.
- Threshold tuning may be necessary depending on vegetation season, terrain, and dump characteristics.

---

## 👤 Authors & Contact

Developed by: OneDrone d.o.o. & ICAERUS Project contributors  
Contact: info@onedrone.com 

