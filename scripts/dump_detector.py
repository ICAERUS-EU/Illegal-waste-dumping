import numpy as np
import os
from qgis.core import (
    QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry, QgsPointXY,
    QgsFields, QgsField, QgsRasterLayer, QgsCoordinateReferenceSystem,
    QgsCoordinateTransformContext, QgsVectorFileWriter, QgsWkbTypes
)
from PyQt5.QtCore import QVariant

def detect_illegal_dumps(ndvi_layer_path, rgb_layer_path, threshold, output_path, min_size_sqm=0.25):
    """
    ndvi_layer_path: path to NDVI GeoTIFF
    rgb_layer_path: path to RGB GeoTIFF
    threshold: NDVI threshold (e.g. 0.2)
    output_path: path to save the output shapefile (dots)
    min_size_sqm: minimum area in square meters (e.g. 0.25 m² for 50x50 cm)
    """
    import rasterio
    from rasterio.features import shapes
    from shapely.geometry import shape, Point

    # Open NDVI raster
    with rasterio.open(ndvi_layer_path) as src:
        ndvi = src.read(1)
        transform = src.transform
        crs = src.crs

        # Mask areas below threshold
        mask = ndvi < threshold

        # Extract polygons from mask
        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for s, v in shapes(mask.astype(np.uint8), mask=mask, transform=transform)
        )

        # Filter small areas and convert to points
        spots = []
        for result in results:
            geom = shape(result["geometry"])
            area = geom.area  # in CRS units (likely meters²)
            if area >= min_size_sqm:
                centroid = geom.centroid
                spots.append((centroid.x, centroid.y, area))

    # Create vector layer with dots
    fields = QgsFields()
    fields.append(QgsField("id", QVariant.Int))
    fields.append(QgsField("area_m2", QVariant.Double))

    mem_layer = QgsVectorLayer("Point?crs=" + crs.to_string(), "Illegal_Dumps", "memory")
    mem_layer.dataProvider().addAttributes(fields)
    mem_layer.updateFields()

    for i, (x, y, area) in enumerate(spots):
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
        feat.setAttributes([i + 1, area])
        mem_layer.dataProvider().addFeature(feat)

    mem_layer.updateExtents()
    QgsProject.instance().addMapLayer(mem_layer)

    # Optionally export as Shapefile
    QgsVectorFileWriter.writeAsVectorFormat(mem_layer, output_path, "UTF-8", mem_layer.crs(), "ESRI Shapefile")

    return mem_layer
