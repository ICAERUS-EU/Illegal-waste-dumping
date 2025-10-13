"""
Dump Detector – v1.5 with built-in Otsu
– ROIcropping with crop=True to minimize memory use
– Automatic NDVI threshold if user didn’t specify
– Progress dialog support
– QGIS 3.22+ V3 vector writer
– Proper rasterio → QGIS CRS conversion
"""

from typing import Optional
import numpy as np
import rasterio
from rasterio.features import shapes
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
import rasterio.mask as _mask
from shapely import wkt as _wkt
from shapely.geometry import shape

from qgis.core import (
    QgsProject,
    QgsVectorLayer,
    QgsFeature,
    QgsGeometry,
    QgsFields,
    QgsField,
    QgsSingleSymbolRenderer,
    QgsSymbol,
    QgsVectorFileWriter,
    QgsCoordinateReferenceSystem
)
from PyQt5.QtCore import QVariant
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QMessageBox, QProgressDialog, QFileDialog


def threshold_otsu(a: np.ndarray, nbins: int = 256) -> float:
    """Compute Otsu’s threshold on a 1D array (ignores NaNs)."""
    arr = a[~np.isnan(a)]
    hist, bin_edges = np.histogram(arr, bins=nbins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    weight1 = np.cumsum(hist)
    weight2 = weight1[-1] - weight1
    cumsum_mean = np.cumsum(hist * bin_centers)
    mean1 = cumsum_mean / weight1
    mean2 = (cumsum_mean[-1] - cumsum_mean) / weight2
    var_between = weight1 * weight2 * (mean1 - mean2) ** 2
    return bin_centers[np.nanargmax(var_between)]


def _compute_vari(rgb: np.ndarray) -> np.ndarray:
    """Compute VARI = (G - R) / (G + R - B)."""
    r, g, b = (rgb[i].astype(np.float32) for i in range(3))
    denom = g + r - b
    denom[denom == 0] = 1e-6
    return (g - r) / denom


def detect_illegal_dumps(
    ndvi_path: str,
    rgb_path: str,
    ndvi_threshold: Optional[float],
    shp_out: str,
    vari_threshold: float = 0.0,
    min_area: float = 0.25,
    roi_wkt: Optional[str] = None,
    *,
    method: str = "clip",
    progress: Optional[QProgressDialog] = None
) -> tuple[Optional[QgsVectorLayer], float]:
    """
    Detect potential illegal-dump polygons based on NDVI and RGB rasters.

    Parameters
    ----------
    ndvi_path : str
        Path to a single-band NDVI raster.  Values outside ``[-1,1]`` are ignored.
    rgb_path : str
        Path to a three-band RGB raster to derive the VARI index.
    ndvi_threshold : float or None
        Threshold for NDVI; pixels with NDVI < threshold are considered
        candidate waste.  If ``None`` then an Otsu threshold is computed
        automatically on the valid NDVI pixels.
    shp_out : str
        Destination path for the output shapefile.  Will be overwritten.
    vari_threshold : float, optional
        Threshold for VARI; pixels with VARI < threshold are retained.  Defaults
        to 0.0.
    min_area : float, optional
        Minimum polygon area (in map units) to keep.  Defaults to 0.25.
    roi_wkt : str, optional
        An optional WKT polygon string defining a region of interest.  Only
        polygons inside this geometry are considered.  If ``None`` then
        the entire NDVI raster is processed.
    method : {'clip','extract'}, optional
        How to handle the ROI geometry when provided:

        * ``'clip'`` (default): crop the NDVI raster to the bounding box of the ROI.
          This reduces processing to the minimal extent.  The output polygons are
          derived from the cropped image and coordinates reflect the cropped
          extent.
        * ``'extract'``: mask pixels outside the ROI but keep the raster at its
          original extent.  This preserves the full image dimensions so that
          polygon coordinates align with the full raster.  Pixels outside the
          ROI are set to NaN and ignored during detection.
    progress : QProgressDialog, optional
        A progress dialog for feedback and cancellation.

    Returns
    -------
    layer : QgsVectorLayer or None
        A memory vector layer containing candidate dump polygons with fields
        ``id``, ``area_m2`` and ``status``, or ``None`` if no candidates
        were found.  The layer is also added to the current QGIS project.
    used_threshold : float
        The NDVI threshold that was applied.  When automatic thresholding
        is enabled this will be the computed Otsu value, otherwise it
        matches the provided ``ndvi_threshold`` argument.
    """
    # 0% load NDVI
    if progress:
        progress.setValue(0)
        progress.setLabelText("Loading NDVI…")

    with rasterio.open(ndvi_path) as ndvi_ds:
        # Read the full NDVI raster into memory as float32.  Values outside
        # the valid range will remain as-is and can later be thresholded
        # appropriately.  Retain the raster's CRS and transform for later use.
        ndvi_full: np.ndarray = ndvi_ds.read(1).astype(np.float32)
        raster_crs = ndvi_ds.crs
        dst_transform = ndvi_ds.transform
        h, w = ndvi_ds.height, ndvi_ds.width

        # If an ROI is provided then optionally clip or mask the NDVI.  When
        # cropping, the transform and dimensions are updated to match the
        # cropped area; when masking, the transform and dimensions remain
        # unchanged but pixels outside the ROI are filled with NaN so they
        # fall below any sensible threshold and are ignored.
        if roi_wkt:
            roi_geom = _wkt.loads(roi_wkt)
            if method not in {"clip", "extract"}:
                raise ValueError(f"Invalid method '{method}'; expected 'clip' or 'extract'.")
            if method == "clip":
                clipped, out_transform = _mask.mask(
                    ndvi_ds,
                    [roi_geom.__geo_interface__],
                    crop=True,
                    filled=True,
                    nodata=np.nan
                )
                ndvi = clipped[0].astype(np.float32)
                dst_transform = out_transform
                h, w = ndvi.shape
            else:  # extract
                masked, _ = _mask.mask(
                    ndvi_ds,
                    [roi_geom.__geo_interface__],
                    crop=False,
                    filled=True,
                    nodata=np.nan
                )
                ndvi = masked[0].astype(np.float32)
                # keep original transform and dimensions
        else:
            ndvi = ndvi_full

    # auto-Otsu if needed
    if ndvi_threshold is None:
        valid = ndvi[~np.isnan(ndvi)]
        if valid.size == 0:
            raise ValueError("No valid NDVI pixels for auto threshold!")
        ndvi_threshold = float(threshold_otsu(valid))

    # rasterio CRS → QGIS CRS
    epsg = raster_crs.to_epsg()
    if epsg:
        qgs_crs = QgsCoordinateReferenceSystem(f"EPSG:{epsg}")
    else:
        qgs_crs = QgsCoordinateReferenceSystem()
        qgs_crs.createFromWkt(raster_crs.to_wkt())

    # 10% warp RGB
    if progress:
        progress.setValue(10)
        progress.setLabelText("Warping RGB…")
        if progress.wasCanceled():
            raise RuntimeError("User aborted.")

    with rasterio.open(rgb_path) as rgb_src:
        if rgb_src.crs != raster_crs:
            raise ValueError("CRS mismatch between NDVI and RGB.")
        with WarpedVRT(
            rgb_src,
            crs=raster_crs,
            transform=dst_transform,
            width=w, height=h,
            resampling=Resampling.bilinear
        ) as vrt:
            if vrt.count < 3:
                raise ValueError(f"RGB needs ≥3 bands; got {vrt.count}")
            rgb = vrt.read([1, 2, 3]).astype(np.float32)

    # 30% compute indices
    if progress:
        progress.setValue(30)
        progress.setLabelText("Computing indices…")
        if progress.wasCanceled():
            raise RuntimeError("User aborted.")

    vari = _compute_vari(rgb)
    # For 'extract' method, mask VARI outside the ROI where NDVI is NaN to
    # ensure those pixels are not considered during thresholding
    if roi_wkt and method == "extract":
        vari[np.isnan(ndvi)] = np.nan
    mask = (ndvi < ndvi_threshold) & (vari < vari_threshold)

    # 50% extract polygons
    if progress:
        progress.setValue(50)
        progress.setLabelText("Extracting polygons…")

    polys = []
    for geom_dict, _ in shapes(mask.astype(np.uint8), mask=mask, transform=dst_transform):
        poly = shape(geom_dict)
        if poly.area >= min_area:
            polys.append(poly)
        if progress and progress.wasCanceled():
            raise RuntimeError("User aborted.")

    if not polys:
        if progress:
            progress.close()
        # Even if no polygons are found, return the used threshold for
        # consistency with the return type.
        return None, float(ndvi_threshold)

    # 80% build memory layer
    if progress:
        progress.setValue(80)
        progress.setLabelText("Building layer…")

    fields = QgsFields()
    fields.append(QgsField("id", QVariant.Int))
    fields.append(QgsField("area_m2", QVariant.Double))
    fields.append(QgsField("status", QVariant.String))

    layer = QgsVectorLayer(f"Polygon?crs={qgs_crs.authid()}", "IWD_polygons", "memory")
    pr = layer.dataProvider()
    pr.addAttributes(fields)
    layer.updateFields()

    for idx, poly in enumerate(polys, start=1):
        # Compute the minimum rotated rectangle around the polygon
        rect = poly.minimum_rotated_rectangle
        # Filter out highly elongated rectangles (aspect ratio > 5)
        minx, miny, maxx, maxy = rect.bounds
        width = maxx - minx
        height = maxy - miny
        # Avoid division by zero
        if width == 0 or height == 0:
            continue
        ratio = max(width, height) / min(width, height)
        if ratio > 5.0:
            continue
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromWkt(rect.wkt))
        # Store area based on the original polygon, not the bounding box
        feat.setAttributes([idx, poly.area, None])
        pr.addFeature(feat)

    layer.updateExtents()
    QgsProject.instance().addMapLayer(layer)

    # Store the RGB source on the layer so that review can pick it up for PDF rendering
    try:
        layer.setCustomProperty('iwd_rgb_src', rgb_path)
    except Exception:
        pass

    # save with V3 writer
    tx = QgsProject.instance().transformContext()
    opts = QgsVectorFileWriter.SaveVectorOptions()
    opts.driverName = "ESRI Shapefile"
    opts.fileEncoding = "UTF-8"
    res = QgsVectorFileWriter.writeAsVectorFormatV3(layer, shp_out, tx, opts)
    err, msg = res[0], (res[1] if len(res)>1 else "")
    if err != QgsVectorFileWriter.NoError:
        raise RuntimeError(f"Failed to write {shp_out}: {msg}")

    if progress:
        progress.setValue(100)
        progress.close()

    return layer, float(ndvi_threshold)


def review_dumps(layer: QgsVectorLayer, iface):
    """Interactively review candidate dumps, save verified shapefile+.

    The user is prompted for each candidate polygon whether it represents waste.
    The prompt dialog is positioned in the top-right corner of the QGIS main
    window so that it does not obscure the map view.  After reviewing all
    candidates, a verified shapefile is written to disk and individual PDF
    reports are generated for each confirmed dump.
    """
    idx = layer.fields().indexOf("status")
    if idx == -1:
        QMessageBox.warning(iface.mainWindow(), "IWD", "Layer lacks 'status' field.")
        return

    # Preserve the original renderer and draw polygons with red outlines
    orig = layer.renderer().clone()
    sym = QgsSymbol.defaultSymbol(layer.geometryType())
    sym.setColor(QColor(0, 0, 0, 0))
    sl = sym.symbolLayer(0)
    sl.setStrokeColor(QColor(255, 0, 0))
    sl.setStrokeWidth(0.8)
    layer.setRenderer(QgsSingleSymbolRenderer(sym))
    layer.triggerRepaint()

    canvas = iface.mapCanvas()
    layer.startEditing()
    aborted = False

    # Iterate features for review
    for feat in layer.getFeatures():
        # Skip those already reviewed
        if feat[idx]:
            continue
        # Zoom to the feature's bounding box
        canvas.setExtent(feat.geometry().boundingBox())
        canvas.refresh()
        # Build a custom message box to reposition
        msg = QMessageBox(iface.mainWindow())
        msg.setWindowTitle("Dump review")
        msg.setIcon(QMessageBox.Question)
        msg.setText(
            f"Feature {feat['id']} (area {feat['area_m2']:.2f} m²)\nIs this waste?"
        )
        msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
        # Move the message box to the top-right corner of the main window
        main_win = iface.mainWindow()
        # Ensure size is calculated before moving
        msg.adjustSize()
        # Compute coordinates with 50 px margin
        x = main_win.x() + main_win.width() - msg.width() - 50
        y = main_win.y() + 50
        msg.move(x, y)
        res = msg.exec()
        if res == QMessageBox.Cancel:
            aborted = True
            break
        layer.changeAttributeValue(
            feat.id(), idx, "waste" if res == QMessageBox.Yes else "clean"
        )

    # Commit edits and restore renderer
    layer.commitChanges()
    layer.setRenderer(orig)
    layer.triggerRepaint()

    if aborted:
        QMessageBox.information(
            iface.mainWindow(), "IWD", "Review aborted; partial results saved."
        )
        return

    # Build a confirmed layer containing only features marked as waste
    confirmed = QgsVectorLayer(
        f"Polygon?crs={layer.crs().authid()}", "IWD_confirmed", "memory"
    )
    pr2 = confirmed.dataProvider()
    pr2.addAttributes(layer.fields())
    confirmed.updateFields()
    for feat in layer.getFeatures():
        if feat[idx] == "waste":
            pr2.addFeature(feat)
    confirmed.updateExtents()
    QgsProject.instance().addMapLayer(confirmed)

    # Ask the user where to save the verified shapefile
    save_path, _ = QFileDialog.getSaveFileName(
        iface.mainWindow(), "Save verified shapefile", "", "ESRI Shapefile (*.shp)"
    )
    if not save_path:
        return
    if not save_path.lower().endswith(".shp"):
        save_path += ".shp"

    tx = QgsProject.instance().transformContext()
    opts = QgsVectorFileWriter.SaveVectorOptions()
    opts.driverName = "ESRI Shapefile"
    opts.fileEncoding = "UTF-8"
    res = QgsVectorFileWriter.writeAsVectorFormatV3(confirmed, save_path, tx, opts)
    err, msg = res[0], (res[1] if len(res) > 1 else "")
    if err != QgsVectorFileWriter.NoError:
        raise RuntimeError(f"Failed to write {save_path}: {msg}")

   

    # Restore original canvas extent and refresh the view.  This prevents
    # altering the user's map view during report generation.
    try:
        canvas.setExtent(original_extent)
        canvas.refresh()
    except Exception:
        pass
