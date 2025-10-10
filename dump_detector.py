#"""
#Dump Detector – v0.8  (QGIS-3.22 safe)
#=======================================
#• Uses transparent fill + red stroke instead of setFillStyle (not available in 3.22).
#• Progress dialog unchanged (Optional[QProgressDialog]).
#"""
#
#from typing import Optional
#import numpy as np
#from qgis.core import (
#    QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry,
#    QgsFields, QgsField, QgsVectorFileWriter, QgsSingleSymbolRenderer, QgsSymbol
#)
#from PyQt5.QtCore import QVariant
#from PyQt5.QtGui import QColor
#from PyQt5.QtWidgets import QMessageBox, QProgressDialog
#
## -----------------------------------------------------------------------------
## Helper
## -----------------------------------------------------------------------------
#
#def _compute_vari(rgb: np.ndarray) -> np.ndarray:
#    r, g, b = (rgb[i].astype(np.float32) for i in range(3))
#    denom = r + g - b
#    denom[denom == 0] = 1e-6
#    return (g - r) / denom
#
## -----------------------------------------------------------------------------
## Detection (progress-aware)
## -----------------------------------------------------------------------------
#
#def detect_illegal_dumps(ndvi_path: str, rgb_path: str, ndvi_thr: float,
#                         shp_out: str, vari_threshold: float = 0.0,
#                         min_area: float = 0.25,
#                         roi_wkt: str | None = None,
#                         progress: Optional[QProgressDialog] = None):
#    if progress:
#        progress.setValue(0); progress.setLabelText("Loading NDVI…")
#
#    import rasterio
#    from rasterio.features import shapes
#    from rasterio.vrt import WarpedVRT
#    from rasterio.enums import Resampling
#    from shapely.geometry import shape
#
#    with rasterio.open(ndvi_path) as ndvi_ds:
#        ndvi = ndvi_ds.read(1).astype(np.float32)
#        dst_crs = ndvi_ds.crs
#        dst_transform = ndvi_ds.transform
#        h, w = ndvi_ds.height, ndvi_ds.width
#
#    if progress:
#        progress.setValue(10); progress.setLabelText("Warping RGB…")
#        if progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#
#    with rasterio.open(rgb_path) as rgb_src:
#        if rgb_src.crs != dst_crs:
#            raise ValueError("CRS mismatch between NDVI and RGB rasters.")
#        with WarpedVRT(rgb_src, crs=dst_crs, transform=dst_transform,
#                       width=w, height=h, resampling=Resampling.bilinear) as vrt:
#            rgb = vrt.read([1, 2, 3])
#
#    if progress:
#        progress.setValue(30); progress.setLabelText("Computing indices…")
#        if progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#
#    vari = _compute_vari(rgb)
#    mask = (ndvi < ndvi_thr) & (vari < vari_threshold)
#
#    if progress:
#        progress.setValue(50); progress.setLabelText("Extracting polygons…")
#
#    geoms = (shape(g) for g, _ in shapes(mask.astype(np.uint8), mask=mask, transform=dst_transform))
#    polys = []
#    for g in geoms:
#        if progress and progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#        if g.area >= min_area:
#            polys.append(g)
#
#    if not polys:
#        if progress: progress.close()
#        return None
#
#    if progress:
#        progress.setValue(80); progress.setLabelText("Writing outputs…")
#
#    fields = QgsFields(); fields.append(QgsField("id", QVariant.Int)); fields.append(QgsField("area_m2", QVariant.Double)); fields.append(QgsField("status", QVariant.String))
#    layer = QgsVectorLayer(f"Polygon?crs={dst_crs.to_string()}", "IWD_polygons", "memory")
#    pr = layer.dataProvider(); pr.addAttributes(fields); layer.updateFields()
#    for idx, geom in enumerate(polys, 1):
#        feat = QgsFeature(); feat.setGeometry(QgsGeometry.fromWkt(geom.wkt)); feat.setAttributes([idx, geom.area, None]); pr.addFeature(feat)
#    layer.updateExtents(); QgsProject.instance().addMapLayer(layer)
#    QgsVectorFileWriter.writeAsVectorFormat(layer, shp_out, "UTF-8", layer.crs(), "ESRI Shapefile")
#
#    if progress:
#        progress.setValue(100); progress.close()
#    return layer
#
## -----------------------------------------------------------------------------
## Review – 3.22-compatible outline
## -----------------------------------------------------------------------------
#
#def review_dumps(layer: QgsVectorLayer, iface):
#    idx = layer.fields().indexOf("status")
#    if idx == -1:
#        QMessageBox.warning(iface.mainWindow(), "IWD", "Layer lacks 'status' field.")
#        return
#
#    orig = layer.renderer().clone()
#    sym = QgsSymbol.defaultSymbol(layer.geometryType())
#    sym.setColor(QColor(0, 0, 0, 0))           # transparent fill
#    sl = sym.symbolLayer(0)
#    sl.setStrokeColor(QColor(255, 0, 0))
#    sl.setStrokeWidth(0.8)
#    layer.setRenderer(QgsSingleSymbolRenderer(sym)); layer.triggerRepaint()
#
#    canvas = iface.mapCanvas(); layer.startEditing(); aborted = False
#    for f in layer.getFeatures():
#        if f[idx]:
#            continue
#        canvas.setExtent(f.geometry().boundingBox()); canvas.refresh()
#        res = QMessageBox.question(
#            iface.mainWindow(), "Dump review",
#            f"Feature {f['id']} (area {f['area_m2']:.2f} m²)\nIs this waste?",
#            QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
#        if res == QMessageBox.Cancel:
#            aborted = True; break
#        layer.changeAttributeValue(f.id(), idx, "waste" if res == QMessageBox.Yes else "clean")
#
#    layer.commitChanges(); layer.setRenderer(orig); layer.triggerRepaint()
#    if aborted:
#        QMessageBox.information(iface.mainWindow(), "IWD", "Review aborted; partial results saved.")
#        return
#
#    confirmed = QgsVectorLayer("Polygon?crs=" + layer.crs().to_string(), "IWD_confirmed", "memory")
#    pr2 = confirmed.dataProvider(); pr2.addAttributes(layer.fields()); confirmed.updateFields()
#    for f in layer.getFeatures():
#        if f[idx] == "waste":
#            pr2.addFeature(f)
#    confirmed.updateExtents(); QgsProject.instance().addMapLayer(confirmed)
#    out = layer.source().replace('.shp', '_verified.shp')
#    QgsVectorFileWriter.writeAsVectorFormat(confirmed, out, "UTF-8", confirmed.crs(), "ESRI Shapefile")
#"""
#Dump Detector – final v1.0
#==========================
#✔ ROI NDVI clipping via _mask
#✔ No local assignment to rasterio
#✔ Progress and review unchanged
#"""
#
#from typing import Optional
#import numpy as np
#from qgis.core import (
#    QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry,
#    QgsFields, QgsField, QgsVectorFileWriter, QgsSingleSymbolRenderer, QgsSymbol
#)
#from PyQt5.QtCore import QVariant
#from PyQt5.QtGui import QColor
#from PyQt5.QtWidgets import QMessageBox, QProgressDialog
#
#from shapely import wkt as _wkt
#from shapely.geometry import shape
#import rasterio
#from rasterio.features import shapes
#from rasterio.vrt import WarpedVRT
#from rasterio.enums import Resampling
#import rasterio.mask as _mask
#
#
#def _compute_vari(rgb: np.ndarray) -> np.ndarray:
#    r, g, b = (rgb[i].astype(np.float32) for i in range(3))
#    denom = g + r - b
#    denom[denom == 0] = 1e-6
#    return (g - r) / denom
#
#
#def detect_illegal_dumps(
#    ndvi_path: str,
#    rgb_path: str,
#    ndvi_threshold: float,
#    shp_out: str,
#    vari_threshold: float = 0.0,
#    min_area: float = 0.25,
#    roi_wkt: Optional[str] = None,
#    progress: Optional[QProgressDialog] = None
#) -> Optional[QgsVectorLayer]:
#    """Detect dumps with optional ROI clipping and update progress."""
#    # 0%: load NDVI
#    if progress:
#        progress.setValue(0);
#        progress.setLabelText("Loading NDVI…")
#
#    with rasterio.open(ndvi_path) as ndvi_ds:
#        ndvi = ndvi_ds.read(1).astype(np.float32)
#        dst_crs = ndvi_ds.crs
#        dst_transform = ndvi_ds.transform
#        h, w = ndvi_ds.height, ndvi_ds.width
#
#        # Clip NDVI to ROI if given
#        if roi_wkt:
#            roi_geom = _wkt.loads(roi_wkt)
#            clipped, _ = _mask.mask(
#                ndvi_ds,
#                [roi_geom.__geo_interface__],
#                crop=False,
#                filled=True,
#                nodata=np.nan
#            )
#            ndvi = clipped[0].astype(np.float32)
#
#    # 10%: warp RGB
#    if progress:
#        progress.setValue(10);
#        progress.setLabelText("Warping RGB…");
#        if progress.wasCanceled(): raise RuntimeError("User aborted.")
#
#    with rasterio.open(rgb_path) as rgb_src:
#        if rgb_src.crs != dst_crs:
#            raise ValueError("CRS mismatch between NDVI and RGB.")
#        with WarpedVRT(
#            rgb_src,
#            crs=dst_crs,
#            transform=dst_transform,
#            width=w,
#            height=h,
#            resampling=Resampling.bilinear
#        ) as vrt:
#            rgb = vrt.read([1, 2, 3]).astype(np.float32)
#
#    # 30%: compute indices
#    if progress:
#        progress.setValue(30);
#        progress.setLabelText("Computing indices…");
#        if progress.wasCanceled(): raise RuntimeError("User aborted.")
#
#    vari = _compute_vari(rgb)
#    mask = (ndvi < ndvi_threshold) & (vari < vari_threshold)
#
#    # 50%: extract polygons
#    if progress:
#        progress.setValue(50);
#        progress.setLabelText("Extracting polygons…")
#
#    polys = []
#    for geom_dict, val in shapes(mask.astype(np.uint8), mask=mask, transform=dst_transform):
#        poly = shape(geom_dict)
#        if poly.area >= min_area:
#            polys.append(poly)
#        if progress and progress.wasCanceled(): raise RuntimeError("User aborted.")
#
#    if not polys:
#        if progress: progress.close()
#        return None
#
#    # 80%: build layer
#    if progress:
#        progress.setValue(80);
#        progress.setLabelText("Writing outputs…")
#
#    fields = QgsFields()
#    fields.append(QgsField("id", QVariant.Int))
#    fields.append(QgsField("area_m2", QVariant.Double))
#    fields.append(QgsField("status", QVariant.String))
#
#    layer = QgsVectorLayer(f"Polygon?crs={dst_crs.to_string()}", "IWD_polygons", "memory")
#    pr = layer.dataProvider(); pr.addAttributes(fields); layer.updateFields()
#
#    for idx, poly in enumerate(polys, 1):
#        f = QgsFeature()
#        f.setGeometry(QgsGeometry.fromWkt(poly.wkt))
#        f.setAttributes([idx, poly.area, None])
#        pr.addFeature(f)
#
#    layer.updateExtents(); QgsProject.instance().addMapLayer(layer)
#    QgsVectorFileWriter.writeAsVectorFile(layer,shp_out,layer.crs(),"ESRI Shapefile",onlySelected=False,encoding="UTF-8")
#
#    if progress:
#        progress.setValue(100);
#        progress.close()
#
#    return layer
#
#
#def review_dumps(layer: QgsVectorLayer, iface):
#    idx = layer.fields().indexOf("status")
#    if idx == -1:
#        QMessageBox.warning(iface.mainWindow(), "IWD", "Layer lacks 'status' field.")
#        return
#
#    orig = layer.renderer().clone()
#    sym = QgsSymbol.defaultSymbol(layer.geometryType())
#    sym.setColor(QColor(0, 0, 0, 0))
#    sl = sym.symbolLayer(0); sl.setStrokeColor(QColor(255, 0, 0)); sl.setStrokeWidth(0.8)
#    layer.setRenderer(QgsSingleSymbolRenderer(sym)); layer.triggerRepaint()
#
#    canvas = iface.mapCanvas(); layer.startEditing(); aborted = False
#    for feat in layer.getFeatures():
#        if feat[idx]: continue
#        canvas.setExtent(feat.geometry().boundingBox()); canvas.refresh()
#        res = QMessageBox.question(iface.mainWindow(), "Dump review",
#                                   f"Feature {feat['id']} (area {feat['area_m2']:.2f} m²)\nIs this waste?",
#                                   QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel)
#        if res == QMessageBox.Cancel:
#            aborted = True; break
#        layer.changeAttributeValue(feat.id(), idx, "waste" if res == QMessageBox.Yes else "clean")
#
#    layer.commitChanges(); layer.setRenderer(orig); layer.triggerRepaint()
#    if aborted:
#        QMessageBox.information(iface.mainWindow(), "IWD", "Review aborted; partial results saved.")
#        return
#
#    confirmed = QgsVectorLayer(f"Polygon?crs={layer.crs().authid()}", "IWD_confirmed", "memory")
#    pr2 = confirmed.dataProvider(); pr2.addAttributes(layer.fields()); confirmed.updateFields()
#    for feat in layer.getFeatures():
#        if feat[idx] == "waste": pr2.addFeature(feat)
#    confirmed.updateExtents(); QgsProject.instance().addMapLayer(confirmed)
#    out = layer.source().replace('.shp', '_verified.shp')
#    QgsVectorFileWriter.writeAsVectorFormat(confirmed, out, "UTF-8", confirmed.crs(), "ESRI Shapefile")
#
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
#"""
#Dump Detector – v1.4
#– ROI clipping on NDVI
#– Progress dialog support
#– QGIS 3.22+ V3 vector writer
#– Proper rasterio → QGIS CRS conversion
#"""
#
#from typing import Optional
#import numpy as np
#from qgis.core import (
#    QgsProject,
#    QgsVectorLayer,
#    QgsFeature,
#    QgsGeometry,
#    QgsFields,
#    QgsField,
#    QgsSingleSymbolRenderer,
#    QgsSymbol,
#    QgsVectorFileWriter,
#    QgsCoordinateReferenceSystem
#)
#from PyQt5.QtCore import QVariant
#from PyQt5.QtGui import QColor
#from PyQt5.QtWidgets import QMessageBox, QProgressDialog, QFileDialog
#
#from shapely import wkt as _wkt
#from shapely.geometry import shape
#import rasterio
#from rasterio.features import shapes
#from rasterio.vrt import WarpedVRT
#from rasterio.enums import Resampling
#import rasterio.mask as _mask
#
#
#def _compute_vari(rgb: np.ndarray) -> np.ndarray:
#    """Compute VARI = (G - R) / (G + R - B)."""
#    r, g, b = (rgb[i].astype(np.float32) for i in range(3))
#    denom = g + r - b
#    denom[denom == 0] = 1e-6
#    return (g - r) / denom
#
#
#def detect_illegal_dumps(
#    ndvi_path: str,
#    rgb_path: str,
#    ndvi_threshold: float,
#    shp_out: str,
#    vari_threshold: float = 0.0,
#    min_area: float = 0.25,
#    roi_wkt: Optional[str] = None,
#    progress: Optional[QProgressDialog] = None
#) -> Optional[QgsVectorLayer]:
#    """Detect potential illegal-dump polygons and save to shapefile."""
#    # 0%: Load NDVI
#    if progress:
#        progress.setValue(0)
#        progress.setLabelText("Loading NDVI…")
#
#    with rasterio.open(ndvi_path) as ndvi_ds:
#        ndvi = ndvi_ds.read(1).astype(np.float32)
#        raster_crs = ndvi_ds.crs
#        dst_transform = ndvi_ds.transform
#        h, w = ndvi_ds.height, ndvi_ds.width
#
#    # Convert rasterio CRS to QGIS CRS
#    epsg = raster_crs.to_epsg()
#    if epsg:
#        qgs_crs = QgsCoordinateReferenceSystem(f"EPSG:{epsg}")
#    else:
#        qgs_crs = QgsCoordinateReferenceSystem()
#        qgs_crs.createFromWkt(raster_crs.to_wkt())
#
#    # Clip NDVI to ROI if provided
#    if roi_wkt:
#        roi_geom = _wkt.loads(roi_wkt)
#        with rasterio.open(ndvi_path) as ds:
#            clipped, out_transform = _mask.mask(
#                ds,
#                [roi_geom.__geo_interface__],
#                crop=True,
#                filled=True,
#                nodata=np.nan
#            )
#            ndvi = clipped[0].astype(np.float32)
#            # zdaj osveži transform in velikost
#            dst_transform = out_transform
#h, w = ndvi.shape
#
#
#    # 10%: Warp RGB onto NDVI grid
#    if progress:
#        progress.setValue(10)
#        progress.setLabelText("Warping RGB…")
#        if progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#
#    with rasterio.open(rgb_path) as rgb_src:
#        if rgb_src.crs != raster_crs:
#            raise ValueError("CRS mismatch between NDVI and RGB rasters.")
#        with WarpedVRT(
#            rgb_src,
#            crs=raster_crs,
#            transform=dst_transform,
#            width=w, height=h,
#            resampling=Resampling.bilinear
#        ) as vrt:
#            if vrt.count < 3:
#                raise ValueError(f"RGB input needs ≥3 bands; got {vrt.count}")
#            rgb = vrt.read([1, 2, 3]).astype(np.float32)
#
#    # 30%: Compute indices
#    if progress:
#        progress.setValue(30)
#        progress.setLabelText("Computing indices…")
#        if progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#
#    vari = _compute_vari(rgb)
#    mask = (ndvi < ndvi_threshold) & (vari < vari_threshold)
#
#    # 50%: Extract polygons
#    if progress:
#        progress.setValue(50)
#        progress.setLabelText("Extracting polygons…")
#
#    polys = []
#    for geom_dict, _ in shapes(mask.astype(np.uint8), mask=mask, transform=dst_transform):
#        poly = shape(geom_dict)
#        if poly.area >= min_area:
#            polys.append(poly)
#        if progress and progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#
#    if not polys:
#        if progress:
#            progress.close()
#        return None
#
#    # 80%: Build memory layer
#    if progress:
#        progress.setValue(80)
#        progress.setLabelText("Building layer…")
#
#    fields = QgsFields()
#    fields.append(QgsField("id", QVariant.Int))
#    fields.append(QgsField("area_m2", QVariant.Double))
#    fields.append(QgsField("status", QVariant.String))
#
#    layer = QgsVectorLayer(f"Polygon?crs={qgs_crs.authid()}", "IWD_polygons", "memory")
#    pr = layer.dataProvider()
#    pr.addAttributes(fields)
#    layer.updateFields()
#
#    for idx, poly in enumerate(polys, start=1):
#        feat = QgsFeature()
#        feat.setGeometry(QgsGeometry.fromWkt(poly.wkt))
#        feat.setAttributes([idx, poly.area, None])
#        pr.addFeature(feat)
#
#    layer.updateExtents()
#    QgsProject.instance().addMapLayer(layer)
#
#    # Save to shapefile using V3 writer
#    tx = QgsProject.instance().transformContext()
#    opts = QgsVectorFileWriter.SaveVectorOptions()
#    opts.driverName = "ESRI Shapefile"
#    opts.fileEncoding = "UTF-8"
#
#     # writeAsVectorFormatV3 can return a tuple of >2 items
#    res = QgsVectorFileWriter.writeAsVectorFormatV3(
#        layer, shp_out, tx, opts
#    )
#    err = res[0]
#    msg = res[1] if len(res) > 1 else ""
#    if err != QgsVectorFileWriter.NoError:
#        raise RuntimeError(f"Failed to write {shp_out}: {msg}")
#
#    if progress:
#        progress.setValue(100)
#        progress.close()
#
#    return layer
#
#
#def review_dumps(layer: QgsVectorLayer, iface):
#    """Prompt user to mark each candidate as waste or clean and save verified."""
#    idx = layer.fields().indexOf("status")
#    if idx == -1:
#        QMessageBox.warning(iface.mainWindow(), "IWD", "Layer lacks 'status' field.")
#        return
#
#    # Outline-only style
#    orig = layer.renderer().clone()
#    sym = QgsSymbol.defaultSymbol(layer.geometryType())
#    sym.setColor(QColor(0, 0, 0, 0))
#    sl = sym.symbolLayer(0)
#    sl.setStrokeColor(QColor(255, 0, 0))
#    sl.setStrokeWidth(0.8)
#    layer.setRenderer(QgsSingleSymbolRenderer(sym))
#    layer.triggerRepaint()
#
#    canvas = iface.mapCanvas()
#    layer.startEditing()
#    aborted = False
#
#    for feat in layer.getFeatures():
#        if feat[idx]:
#            continue
#        canvas.setExtent(feat.geometry().boundingBox())
#        canvas.refresh()
#        res = QMessageBox.question(
#            iface.mainWindow(), "Dump review",
#            f"Feature {feat['id']} (area {feat['area_m2']:.2f} m²)\nIs this waste?",
#            QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
#        )
#        if res == QMessageBox.Cancel:
#            aborted = True
#            break
#        val = "waste" if res == QMessageBox.Yes else "clean"
#        layer.changeAttributeValue(feat.id(), idx, val)
#
#    layer.commitChanges()
#    layer.setRenderer(orig)
#    layer.triggerRepaint()
#
#    if aborted:
#        QMessageBox.information(iface.mainWindow(), "IWD", "Review aborted; partial results saved.")
#        return
#
#    # Build confirmed layer
#    confirmed = QgsVectorLayer(
#        f"Polygon?crs={layer.crs().authid()}", "IWD_confirmed", "memory"
#    )
#    pr2 = confirmed.dataProvider()
#    pr2.addAttributes(layer.fields())
#    confirmed.updateFields()
#
#    for feat in layer.getFeatures():
#        if feat[idx] == "waste":
#            pr2.addFeature(feat)
#
#    confirmed.updateExtents()
#    QgsProject.instance().addMapLayer(confirmed)
#
#    # Ask user where to save
#    save_path, _ = QFileDialog.getSaveFileName(
#        iface.mainWindow(),
#        "Save verified shapefile",
#        "",
#        "ESRI Shapefile (*.shp)"
#    )
#    if not save_path:
#        return
#    if not save_path.lower().endswith(".shp"):
#        save_path += ".shp"
#
#    # Save via V3 writer
#    tx = QgsProject.instance().transformContext()
#    opts = QgsVectorFileWriter.SaveVectorOptions()
#    opts.driverName = "ESRI Shapefile"
#    opts.fileEncoding = "UTF-8"
#
#    res = QgsVectorFileWriter.writeAsVectorFormatV3(
#        confirmed, save_path, tx, opts
#    )
#    err = res[0]
#    msg = res[1] if len(res) > 1 else ""
#    if err != QgsVectorFileWriter.NoError:
#        raise RuntimeError(f"Failed to write {save_path}: {msg}")
#
# -*- coding: utf-8 -*-
#"""
#Dump Detector – v1.5
#– ROI‐cropping with crop=True to minimize memory use
#– Progress dialog support
#– QGIS 3.22+ V3 vector writer
#– Proper rasterio → QGIS CRS conversion
#"""
#
#from typing import Optional
#import numpy as np
#import rasterio
#from rasterio.features import shapes
#from rasterio.vrt import WarpedVRT
#from rasterio.enums import Resampling
#import rasterio.mask as _mask
#from shapely import wkt as _wkt
#from shapely.geometry import shape
#from skimage.filters import threshold_otsu
#
#
#from qgis.core import (
#    QgsProject,
#    QgsVectorLayer,
#    QgsFeature,
#    QgsGeometry,
#    QgsFields,
#    QgsField,
#    QgsSingleSymbolRenderer,
#    QgsSymbol,
#    QgsVectorFileWriter,
#    QgsCoordinateReferenceSystem
#)
#from PyQt5.QtCore import QVariant
#from PyQt5.QtGui import QColor
#from PyQt5.QtWidgets import QMessageBox, QProgressDialog, QFileDialog
#
#
#def _compute_vari(rgb: np.ndarray) -> np.ndarray:
#    """Compute VARI = (G - R) / (G + R - B)."""
#    r, g, b = (rgb[i].astype(np.float32) for i in range(3))
#    denom = g + r - b
#    denom[denom == 0] = 1e-6
#    return (g - r) / denom
#
#
#def detect_illegal_dumps(
#    ndvi_path: str,
#    rgb_path: str,
#    ndvi_threshold: float,
#    shp_out: str,
#    vari_threshold: float = 0.0,
#    min_area: float = 0.25,
#    roi_wkt: Optional[str] = None,
#    progress: Optional[QProgressDialog] = None
#) -> Optional[QgsVectorLayer]:
#    """Detect potential illegal-dump polygons and save to shapefile."""
#
#    # 0%: load NDVI
#    if progress:
#        progress.setValue(0)
#        progress.setLabelText("Loading NDVI…")
#
#    with rasterio.open(ndvi_path) as ndvi_ds:
#        ndvi_full = ndvi_ds.read(1).astype(np.float32)
#        raster_crs = ndvi_ds.crs
#        dst_transform = ndvi_ds.transform
#        h, w = ndvi_ds.height, ndvi_ds.width
#
#        # crop to ROI if requested
#        if roi_wkt:
#            roi_geom = _wkt.loads(roi_wkt)
#            clipped, out_transform = _mask.mask(
#                ndvi_ds,
#                [roi_geom.__geo_interface__],
#                crop=True,
#                filled=True,
#                nodata=np.nan
#            )
#            ndvi = clipped[0].astype(np.float32)
#            dst_transform = out_transform
#            h, w = ndvi.shape
#        else:
#            ndvi = ndvi_full
#
#    # convert rasterio CRS → QGIS CRS
#    epsg = raster_crs.to_epsg()
#    if epsg:
#        qgs_crs = QgsCoordinateReferenceSystem(f"EPSG:{epsg}")
#    else:
#        qgs_crs = QgsCoordinateReferenceSystem()
#        qgs_crs.createFromWkt(raster_crs.to_wkt())
#
#    # 10%: warp RGB onto NDVI grid (ROI-cropped or full)
#    if progress:
#        progress.setValue(10)
#        progress.setLabelText("Warping RGB…")
#        if progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#
#    with rasterio.open(rgb_path) as rgb_src:
#        if rgb_src.crs != raster_crs:
#            raise ValueError("CRS mismatch between NDVI and RGB.")
#        with WarpedVRT(
#            rgb_src,
#            crs=raster_crs,
#            transform=dst_transform,
#            width=w, height=h,
#            resampling=Resampling.bilinear
#        ) as vrt:
#            if vrt.count < 3:
#                raise ValueError(f"RGB input needs ≥3 bands; got {vrt.count}")
#            rgb = vrt.read([1, 2, 3]).astype(np.float32)
#
#    # 30%: compute indices
#    if progress:
#        progress.setValue(30)
#        progress.setLabelText("Computing indices…")
#        if progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#
#    vari = _compute_vari(rgb)
#    mask = (ndvi < ndvi_threshold) & (vari < vari_threshold)
#
#    # 50%: extract polygons
#    if progress:
#        progress.setValue(50)
#        progress.setLabelText("Extracting polygons…")
#
#    polys = []
#    for geom_dict, _ in shapes(mask.astype(np.uint8), mask=mask, transform=dst_transform):
#        poly = shape(geom_dict)
#        if poly.area >= min_area:
#            polys.append(poly)
#        if progress and progress.wasCanceled():
#            raise RuntimeError("User aborted.")
#
#    if not polys:
#        if progress:
#            progress.close()
#        return None
#
#    # 80%: build memory layer
#    if progress:
#        progress.setValue(80)
#        progress.setLabelText("Building layer…")
#
#    fields = QgsFields()
#    fields.append(QgsField("id", QVariant.Int))
#    fields.append(QgsField("area_m2", QVariant.Double))
#    fields.append(QgsField("status", QVariant.String))
#
#    layer = QgsVectorLayer(f"Polygon?crs={qgs_crs.authid()}", "IWD_polygons", "memory")
#    pr = layer.dataProvider()
#    pr.addAttributes(fields)
#    layer.updateFields()
#
#    for idx, poly in enumerate(polys, start=1):
#        rect = poly.minimum_rotated_rectangle
#        feat = QgsFeature()
#        feat.setGeometry(QgsGeometry.fromWkt(rect.wkt))
#        feat.setAttributes([idx, poly.area, None])
#        pr.addFeature(feat)
#
#
#    layer.updateExtents()
#    QgsProject.instance().addMapLayer(layer)
#
#    # save to shapefile via V3 writer
#    tx = QgsProject.instance().transformContext()
#    opts = QgsVectorFileWriter.SaveVectorOptions()
#    opts.driverName = "ESRI Shapefile"
#    opts.fileEncoding = "UTF-8"
#
#    res = QgsVectorFileWriter.writeAsVectorFormatV3(layer, shp_out, tx, opts)
#    err = res[0]
#    msg = res[1] if len(res) > 1 else ""
#    if err != QgsVectorFileWriter.NoError:
#        raise RuntimeError(f"Failed to write {shp_out}: {msg}")
#
#    if progress:
#        progress.setValue(100)
#        progress.close()
#
#    return layer
#
#
#def review_dumps(layer: QgsVectorLayer, iface):
#    """Prompt user to mark dumps and then save verified shapefile."""
#    idx = layer.fields().indexOf("status")
#    if idx == -1:
#        QMessageBox.warning(iface.mainWindow(), "IWD", "Layer lacks 'status' field.")
#        return
#
#    # outline-only style
#    orig = layer.renderer().clone()
#    sym = QgsSymbol.defaultSymbol(layer.geometryType())
#    sym.setColor(QColor(0, 0, 0, 0))
#    sl = sym.symbolLayer(0)
#    sl.setStrokeColor(QColor(255, 0, 0))
#    sl.setStrokeWidth(0.8)
#    layer.setRenderer(QgsSingleSymbolRenderer(sym))
#    layer.triggerRepaint()
#
#    canvas = iface.mapCanvas()
#    layer.startEditing()
#    aborted = False
#
#    for feat in layer.getFeatures():
#        if feat[idx]:
#            continue
#        canvas.setExtent(feat.geometry().boundingBox())
#        canvas.refresh()
#        res = QMessageBox.question(
#            iface.mainWindow(), "Dump review",
#            f"Feature {feat['id']} (area {feat['area_m2']:.2f} m²)\nIs this waste?",
#            QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel
#        )
#        if res == QMessageBox.Cancel:
#            aborted = True
#            break
#        val = "waste" if res == QMessageBox.Yes else "clean"
#        layer.changeAttributeValue(feat.id(), idx, val)
#
#    layer.commitChanges()
#    layer.setRenderer(orig)
#    layer.triggerRepaint()
#
#    if aborted:
#        QMessageBox.information(iface.mainWindow(), "IWD", "Review aborted; partial results saved.")
#        return
#
#    # build confirmed layer
#    confirmed = QgsVectorLayer(f"Polygon?crs={layer.crs().authid()}", "IWD_confirmed", "memory")
#    pr2 = confirmed.dataProvider()
#    pr2.addAttributes(layer.fields())
#    confirmed.updateFields()
#
#    for feat in layer.getFeatures():
#        if feat[idx] == "waste":
#            pr2.addFeature(feat)
#
#    confirmed.updateExtents()
#    QgsProject.instance().addMapLayer(confirmed)
#
#    # prompt save
#    save_path, _ = QFileDialog.getSaveFileName(
#        iface.mainWindow(),
#        "Save verified shapefile",
#        "",
#        "ESRI Shapefile (*.shp)"
#    )
#    if not save_path:
#        return
#    if not save_path.lower().endswith(".shp"):
#        save_path += ".shp"
#
#    # save via V3 writer
#    tx = QgsProject.instance().transformContext()
#    opts = QgsVectorFileWriter.SaveVectorOptions()
#    opts.driverName = "ESRI Shapefile"
#    opts.fileEncoding = "UTF-8"
#
#    res = QgsVectorFileWriter.writeAsVectorFormatV3(confirmed, save_path, tx, opts)
#    err = res[0]
#    msg = res[1] if len(res) > 1 else ""
#    if err != QgsVectorFileWriter.NoError:
#        raise RuntimeError(f"Failed to write {save_path}: {msg}")
#
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Dump Detector – v1.5 with built-in Otsu
– ROI‐cropping with crop=True to minimize memory use
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
    """Interactively review candidate dumps, save verified shapefile and PDFs.

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

    # Generate PDF reports for each confirmed dump
    import os
    output_dir = os.path.dirname(save_path)
    try:
        # Pass along the RGB source stored on the original layer for high‑quality maps
        background_src = layer.customProperty('iwd_rgb_src', None)
        _generate_pdf_reports(confirmed, iface, output_dir, background_src)
    except Exception as e:
        # Log but don't interrupt the user
        iface.messageBar().pushWarning("IWD", f"PDF generation failed: {e}")


def _generate_pdf_reports(
    layer: QgsVectorLayer,
    iface,
    output_dir: str,
    background_src: str | None = None,
) -> None:
    """
    Generate a multi‑page PDF report for all confirmed dumps in a vector layer.

    Each page of the report summarises a single dump with its ID, area, status
    and centroid, followed by a map excerpt showing the dump geometry overlaid
    on the orthophoto or current map canvas.  The map extent is expanded by
    25 % around the feature to provide context.  The resulting PDF is saved
    as ``iwd_report.pdf`` in the provided output directory.  The report is
    always a single file containing one page per dump.

    Parameters
    ----------
    layer : QgsVectorLayer
        The layer containing confirmed waste polygons.  All features should
        have the attribute ``status`` set to ``'waste'``.
    iface : QgsInterface
        Reference to the QGIS interface, used to access the map canvas.
    output_dir : str
        Directory where the generated report will be saved.
    background_src : str, optional
        Path to an RGB orthophoto used as background.  If provided and
        successfully loaded, a temporary raster layer is rendered beneath
        the vector layer for the map excerpt.  Otherwise, the current
        map canvas contents are captured as a fallback.
    """
    import os
    from PyQt5.QtGui import QPdfWriter, QPainter, QFont, QFontMetrics
    from PyQt5.QtCore import QPointF, QRectF, QSize
    from PyQt5.QtWidgets import QApplication
    from qgis.core import (
        QgsRectangle,
        QgsRasterLayer,
        QgsMapSettings,
        QgsMapRendererParallelJob,
        QgsMapLayer,
    )

    # Access the QGIS map canvas and store its original extent so it can be
    # restored after report generation.
    canvas = iface.mapCanvas()
    original_extent = canvas.extent()

    # Attempt to load a background raster layer from ``background_src``.  If
    # loading fails or no source is provided, ``bg_layer`` will remain ``None``
    # and the current map canvas will be used instead.
    bg_layer: QgsRasterLayer | None = None
    if background_src:
        try:
            candidate = QgsRasterLayer(background_src, "IWD_Background")
            if candidate.isValid():
                bg_layer = candidate
        except Exception:
            bg_layer = None

    # Prepare the PDF writer.  Use a high DPI (300) so that inserted map
    # images appear crisp.  The default page size of QPdfWriter is A4
    # portrait when no page size is explicitly set.
    report_path = os.path.join(output_dir, "iwd_report.pdf")
    writer = QPdfWriter(report_path)
    writer.setResolution(300)
    painter = QPainter(writer)
    # Choose a reasonable font size; metrics will determine line spacing.
    font = QFont("Arial", 12)
    painter.setFont(font)
    metrics = QFontMetrics(font)

    # Cache resolution and compute margin and inter‑element spacing in
    # device pixels.  Margins and spacing are defined in inches to produce
    # consistent sizing regardless of PDF page dimensions.  A 0.5 inch margin
    # around the page and a 0.1 inch gap between text and map are used.
    dpi = writer.resolution()
    margin_px = int(0.5 * dpi)
    gap_px = int(0.1 * dpi)

    from qgis.core import QgsProject  # imported here to avoid cyclic imports
    first_page = True
    # Iterate through features in the confirmed layer.  Each feature produces
    # one page in the report.
    for feat in layer.getFeatures():
        dump_id = feat.get("id")
        area = feat.get("area_m2")
        # Compute centroid using feature geometry
        try:
            centroid_pt = feat.geometry().centroid().asPoint()
        except Exception:
            centroid_pt = None
        # Build a context rectangle by expanding the feature's bounding box
        bbox = feat.geometry().boundingBox()
        dx = bbox.width() * 0.25
        dy = bbox.height() * 0.25
        exp_rect = QgsRectangle(
            bbox.xMinimum() - dx,
            bbox.yMinimum() - dy,
            bbox.xMaximum() + dx,
            bbox.yMaximum() + dy,
        )

        # Determine the dimensions of the PDF page in device pixels.  These
        # values change automatically if the writer's page size is modified
        # (e.g., if the user prints to a different format).  The page width
        # and height returned by QPdfWriter include margins defined via
        # ``writer.setPageMargins()``; however we compute our own margins.
        page_width_px = writer.width()
        page_height_px = writer.height()

        # Render the map excerpt.  Prefer a high‑quality rendering using
        # QgsMapRendererParallelJob when a background raster is available.
        image = None
        if bg_layer:
            # Temporarily add the background layer to the project if it is
            # not already present.  QgsMapSettings expects layer IDs which
            # must belong to layers registered in the project.  Remove
            # afterwards to avoid polluting the user's layer list.
            bg_added = False
            try:
                if not QgsProject.instance().mapLayer(bg_layer.id()):
                    QgsProject.instance().addMapLayer(bg_layer, addToLegend=False)
                    bg_added = True
                settings = QgsMapSettings()
                # Set layers by ID – both raster and vector must be known to
                # the project.  Use the vector layer's ID since it is
                # already registered when calling review.
                settings.setLayers([bg_layer.id(), layer.id()])
                settings.setExtent(exp_rect)
                # Determine available map area after margins and text.
                text_block_height = 4 * metrics.height() + gap_px
                avail_width = page_width_px - 2 * margin_px
                avail_height = page_height_px - (margin_px + text_block_height) - margin_px
                if avail_width <= 0 or avail_height <= 0:
                    avail_width = max(avail_width, 1)
                    avail_height = max(avail_height, 1)
                settings.setOutputSize(QSize(avail_width, avail_height))
                settings.setOutputDpi(dpi)
                # Set white background
                from PyQt5.QtGui import QColor
                settings.setBackgroundColor(QColor(255, 255, 255))
                job = QgsMapRendererParallelJob(settings)
                job.start()
                job.waitForFinished()
                image = job.renderedImage()
            except Exception:
                image = None
            finally:
                # Remove the background layer if it was added
                if bg_added:
                    try:
                        QgsProject.instance().removeMapLayer(bg_layer.id())
                    except Exception:
                        pass
        # Fallback: capture the current map canvas.  This is used when
        # ``background_src`` is not provided or failed to load, or when
        # rendering via QgsMapRendererParallelJob fails for some reason.  The
        # canvas is zoomed to the expanded bounding box before capture.
        if not image:
            try:
                canvas.setExtent(exp_rect)
                canvas.refresh()
                QApplication.processEvents()
                # Grab returns a pixmap of the canvas; convert to QImage
                image = canvas.grab().toImage()
            except Exception:
                image = None

        # If multiple pages, explicitly call newPage() before drawing on
        # subsequent pages.  QPdfWriter draws into the same QPainter; it
        # automatically clears the page when newPage() is invoked.
        if not first_page:
            writer.newPage()
        first_page = False

        # Start drawing the page contents.  Begin at the top margin and
        # increment Y coordinate after each line.  Use the font metrics to
        # compute line spacing.  ``y`` holds the top of the next line.
        y = margin_px
        x = margin_px
        # Compose the lines of text to draw
        lines = [
            f"Dump ID: {dump_id}",
            f"Area (m²): {area:.2f}" if area is not None else "Area (m²): N/A",
            "Status: Confirmed",
            (
                f"Centroid: {centroid_pt.x():.2f}, {centroid_pt.y():.2f}"
                if centroid_pt
                else "Centroid: N/A"
            ),
        ]
        for idx, line in enumerate(lines):
            # Draw the baseline of the text at y + ascent to align properly
            painter.drawText(QPointF(x, y + metrics.ascent()), line)
            y += metrics.height()
        # Leave a small gap below the last line before drawing the map
        y += gap_px

        # Compute the available rectangle for the map image
        avail_width = page_width_px - 2 * margin_px
        avail_height = page_height_px - y - margin_px
        if avail_width <= 0 or avail_height <= 0:
            # Avoid zero or negative sizes
            avail_width = max(avail_width, 1)
            avail_height = max(avail_height, 1)
        dest_rect = QRectF(x, y, avail_width, avail_height)

        # If no image was generated, skip drawing the map on this page
        if image is not None and not image.isNull():
            # Calculate aspect ratios to preserve the map's proportions
            img_ratio = image.width() / image.height() if image.height() != 0 else 1.0
            dest_ratio = dest_rect.width() / dest_rect.height() if dest_rect.height() != 0 else 1.0
            # Determine the width and height of the drawn image while
            # maintaining aspect ratio.  We fit the image within dest_rect
            # without cropping.  Excess space remains blank.
            draw_rect = QRectF(dest_rect)
            if img_ratio > dest_ratio:
                # Image is relatively wider than destination.  Fit width.
                new_height = dest_rect.width() / img_ratio
                draw_rect.setHeight(new_height)
            else:
                # Image is relatively taller than destination.  Fit height.
                new_width = dest_rect.height() * img_ratio
                draw_rect.setWidth(new_width)
            painter.drawImage(draw_rect, image)

    # Finish writing the PDF
    painter.end()

    # Restore original canvas extent and refresh the view.  This prevents
    # altering the user's map view during report generation.
    try:
        canvas.setExtent(original_extent)
        canvas.refresh()
    except Exception:
        pass







































































































































































































































































