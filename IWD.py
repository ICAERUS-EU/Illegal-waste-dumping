#from qgis.PyQt.QtWidgets import (
#    QAction, QDialog, QFileDialog, QVBoxLayout, QHBoxLayout, QLabel,
#    QLineEdit, QPushButton, QDoubleSpinBox, QProgressDialog
#)
#from qgis.core import QgsProject
#from .dump_detector import detect_illegal_dumps, review_dumps
#
## ───────────── param-dialog ─────────────
#class DumpDlg(QDialog):
#    def __init__(self, parent=None):
#        super().__init__(parent)
#        self.setWindowTitle("Detect illegal dumps")
#        lay = QVBoxLayout(self)
#
#        def row(txt, *w):
#            h = QHBoxLayout(); h.addWidget(QLabel(txt))
#            for x in w: h.addWidget(x)
#            return h
#
#        self.ndvi_le = QLineEdit(); b1 = QPushButton("…")
#        b1.clicked.connect(lambda: self._pick(self.ndvi_le))
#        lay.addLayout(row("NDVI raster:", self.ndvi_le, b1))
#
#        self.rgb_le = QLineEdit(); b2 = QPushButton("…")
#        b2.clicked.connect(lambda: self._pick(self.rgb_le))
#        lay.addLayout(row("RGB raster:", self.rgb_le, b2))
#
#        self.thr_sb = QDoubleSpinBox(); self.thr_sb.setDecimals(2)
#        self.thr_sb.setSingleStep(0.01); self.thr_sb.setValue(0.30)
#        lay.addLayout(row("NDVI threshold:", self.thr_sb))
#
#        self.out_le = QLineEdit(); b3 = QPushButton("…")
#        b3.clicked.connect(lambda: self._save(self.out_le))
#        lay.addLayout(row("Output shapefile:", self.out_le, b3))
#
#        run = QPushButton("Run"); run.clicked.connect(self.accept)
#        cancel = QPushButton("Cancel"); cancel.clicked.connect(self.reject)
#        lay.addLayout(row("", run, cancel))
#
#    def _pick(self, line):
#        f, _ = QFileDialog.getOpenFileName(self, "Select raster", "", "TIFF (*.tif)")
#        if f: line.setText(f)
#
#    def _save(self, line):
#        f, _ = QFileDialog.getSaveFileName(self, "Save shapefile", "", "Shapefile (*.shp)")
#        if f:
#            if not f.lower().endswith(".shp"): f += ".shp"
#            line.setText(f)
#
## ───────────── plugin ─────────────
#class IWD:
#    def __init__(self, iface):
#        self.iface = iface
#
#    def initGui(self):
#        a1 = QAction("Detect dumps…", self.iface.mainWindow())
#        a1.triggered.connect(self.run_detect)
#        a2 = QAction("Review dumps…", self.iface.mainWindow())
#        a2.triggered.connect(self.run_review)
#        self.iface.addPluginToMenu("&IWD", a1)
#        self.iface.addPluginToMenu("&IWD", a2)
#        self.iface.addToolBarIcon(a1)
#        self._a1, self._a2 = a1, a2               # keep refs
#
#    def unload(self):
#        self.iface.removePluginMenu("&IWD", self._a1)
#        self.iface.removePluginMenu("&IWD", self._a2)
#        self.iface.removeToolBarIcon(self._a1)
#
#    # --- detect ---
#    def run_detect(self):
#        dlg = DumpDlg(self.iface.mainWindow())
#        if dlg.exec() != QDialog.Accepted:
#            return
#        if not dlg.ndvi_le.text() or not dlg.rgb_le.text() or not dlg.out_le.text():
#            self.iface.messageBar().pushWarning("IWD", "All paths required."); return
#
#        prog = QProgressDialog("Detecting…", "Abort", 0, 100, self.iface.mainWindow())
#        prog.setMinimumDuration(0); prog.setAutoClose(False)
#
#        try:
#            layer = detect_illegal_dumps(
#                dlg.ndvi_le.text(), dlg.rgb_le.text(), dlg.thr_sb.value(),
#                dlg.out_le.text(), progress=prog)
#        except RuntimeError:
#            self.iface.messageBar().pushWarning("IWD", "Detection aborted."); return
#        except Exception as e:
#            self.iface.messageBar().pushWarning("IWD", str(e)); return
#
#        if layer:
#            self.iface.messageBar().pushSuccess(
#                "IWD", f"Found {layer.featureCount()} candidate dumps.")
#
#    # --- review ---
#    def run_review(self):
#        active = self.iface.activeLayer()
#        layer = active if (active and active.name() == "IWD_polygons") else None
#        if not layer:
#            lst = QgsProject.instance().mapLayersByName("IWD_polygons")
#            layer = lst[0] if lst else None
#        if not layer:
#            self.iface.messageBar().pushWarning("IWD", "No polygon layer."); return
#        review_dumps(layer, self.iface)
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
IWD QGIS Plugin – Main entrypoint
"""

from qgis.PyQt.QtWidgets import QAction, QDialog, QProgressDialog
from qgis.core import QgsProject
from .IWD_dialog import IWDDialog
from .dump_detector import detect_illegal_dumps, review_dumps


class IWD:
    def __init__(self, iface):
        self.iface = iface
        self.act_detect = None
        self.act_review = None

    def initGui(self):
        # Detect button
        self.act_detect = QAction("Detect dumps…", self.iface.mainWindow())
        self.act_detect.triggered.connect(self.run_detect)
        self.iface.addPluginToMenu("&IWD", self.act_detect)
        self.iface.addToolBarIcon(self.act_detect)

        # Review button
        self.act_review = QAction("Review dumps…", self.iface.mainWindow())
        self.act_review.triggered.connect(self.run_review)
        self.iface.addPluginToMenu("&IWD", self.act_review)

    def unload(self):
        self.iface.removePluginMenu("&IWD", self.act_detect)
        self.iface.removeToolBarIcon(self.act_detect)
        self.iface.removePluginMenu("&IWD", self.act_review)

    def run_detect(self):
        dlg = IWDDialog(self.iface, self.iface.mainWindow())
        if dlg.exec() != QDialog.Accepted:
            return

        # Collect parameters
        params = dlg.values()
        ndvi_src      = params['ndvi_src']
        rgb_src       = params['rgb_src']
        roi_src       = params['roi_src']
        use_extent    = params.get('use_extent', False)
        method        = params['method']
        threshold     = params['ndvi_threshold']
        min_area      = params['min_area']
        output_folder = params['output_folder']
        add_to_map    = params['add_to_map']

        # Convert NDVI and RGB inputs into file paths if the user selected
        # loaded raster layers by name.  When a layer name or ID is provided
        # instead of a filename, resolve it to its data source.  The data
        # source may include provider options separated by '|' – use only
        # the part before the delimiter as the filesystem path.
        import os
        from qgis.core import QgsProject

        def _resolve_raster(src: str) -> str:
            """Attempt to convert a raster layer name or ID to its file path.

            If ``src`` corresponds to a raster layer currently loaded in the
            project then return its data source path.  Otherwise return
            ``src`` unchanged.
            """
            # If the string looks like an existing file, return as‑is
            if os.path.exists(src):
                return src
            # Try to find a layer by ID
            layer = QgsProject.instance().mapLayer(src)
            if not layer:
                # Try by name
                candidates = QgsProject.instance().mapLayersByName(src)
                layer = candidates[0] if candidates else None
            if layer and hasattr(layer, 'source'):
                path = layer.source()
                # Split at provider delimiter '|' to obtain file path
                return path.split('|')[0]
            return src

        ndvi_src = _resolve_raster(ndvi_src)
        rgb_src = _resolve_raster(rgb_src)

        # Basic validation
        if not ndvi_src or not rgb_src or not output_folder:
            self.iface.messageBar().pushWarning(
                "IWD", "NDVI, RGB and Output folder are required."
            )
            return

        # Build ROI WKT based on the selected option
        roi_wkt = None
        if use_extent:
            # Use the current map canvas extent as the ROI
            extent = self.iface.mapCanvas().extent()
            from shapely.geometry import box
            roi_wkt = box(extent.xMinimum(), extent.yMinimum(),
                          extent.xMaximum(), extent.yMaximum()).wkt
        elif roi_src:
            # Use the selected polygon layer as ROI (unify all polygons)
            candidates = QgsProject.instance().mapLayersByName(roi_src)
            if candidates:
                layer = candidates[0]
                feats = list(layer.getFeatures())
                if feats:
                    from shapely import wkt as _wkt
                    from shapely.ops import unary_union
                    geoms = [_wkt.loads(f.geometry().asWkt()) for f in feats]
                    roi_wkt = unary_union(geoms).wkt

        # Show progress dialog
        prog = QProgressDialog("Detecting dumps…", "Abort", 0, 100,
                               self.iface.mainWindow())
        prog.setMinimumDuration(0)
        prog.setAutoClose(False)
        prog.show()

        try:
            # detect_illegal_dumps returns a tuple (layer, used_threshold)
            layer, used_thr = detect_illegal_dumps(
                ndvi_path      = ndvi_src,
                rgb_path       = rgb_src,
                ndvi_threshold = threshold,
                shp_out        = f"{output_folder}/dumps.shp",
                vari_threshold = 0.0,
                min_area       = min_area,
                roi_wkt        = roi_wkt,
                method         = method,
                progress       = prog
            )

        except RuntimeError:
            self.iface.messageBar().pushWarning("IWD", "Detection aborted.")
            return
        except Exception as e:
            self.iface.messageBar().pushWarning("IWD", str(e))
            return

        # Optionally add the resulting layer to the map
        if layer and add_to_map:
            QgsProject.instance().addMapLayer(layer)

        # Report results: either the count of dumps or a message that none were found
        thr_msg = f"{used_thr:.2f}" if used_thr is not None else "?"
        if layer:
            cnt = layer.featureCount()
            self.iface.messageBar().pushSuccess(
                "IWD", f"Found {cnt} dumps. NDVI threshold was {thr_msg}."
            )
        else:
            # No polygons detected
            self.iface.messageBar().pushInfo(
                "IWD", f"No dumps detected. NDVI threshold was {thr_msg}."
            )

    def run_review(self):
        layer = self.iface.activeLayer()
        if not layer or layer.name() != "IWD_polygons":
            layers = QgsProject.instance().mapLayersByName("IWD_polygons")
            layer = layers[0] if layers else None
        if not layer:
            self.iface.messageBar().pushWarning(
                "IWD", "No polygon layer to review."
            )
            return

        review_dumps(layer, self.iface)
