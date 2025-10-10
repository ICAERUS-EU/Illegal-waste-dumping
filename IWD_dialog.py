## -*- coding: utf-8 -*-
#"""
#IWD QGIS Plugin – Simplified Detect Dialog
#Streamlined UI: wider inputs, small browse buttons, removed split and NoData.
#"""
#
#from qgis.PyQt.QtWidgets import (
#    QDialog, QTabWidget, QWidget, QVBoxLayout, QFormLayout, QHBoxLayout,
#    QLabel, QLineEdit, QComboBox, QLineEdit, QPushButton, QDoubleSpinBox,
#    QDialogButtonBox, QFileDialog, QTextBrowser, QSizePolicy, QCheckBox
#)
#
#from qgis.core import QgsProject, QgsRasterLayer, QgsVectorLayer, QgsWkbTypes
#
#
#class IWDDialog(QDialog):
#    """Parameter dialog: raster inputs, ROI, method, threshold, output."""
#    def __init__(self, iface, parent=None):
#        super().__init__(parent)
#        self.iface = iface
#        self.setWindowTitle("Detect illegal dumps")
#        self.resize(500, 360)
#
#        main_layout = QVBoxLayout(self)
#        tabs = QTabWidget(); main_layout.addWidget(tabs)
#
#        # Parameters tab
#        param_tab = QWidget()
#        form = QFormLayout(param_tab)
#        tabs.addTab(param_tab, "Parameters")
#
#        # NDVI index raster
#        self.ndvi_cb = QComboBox(); self.ndvi_cb.setEditable(True)
#        self.ndvi_cb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
#        self._populate_layers(self.ndvi_cb, raster=True)
#        ndvi_btn = QPushButton("…"); ndvi_btn.setFixedWidth(20)
#        ndvi_btn.clicked.connect(lambda: self._browse(self.ndvi_cb, raster=True))
#        h_ndvi = QHBoxLayout(); h_ndvi.addWidget(self.ndvi_cb); h_ndvi.addWidget(ndvi_btn)
#        form.addRow("NDVI index raster:", h_ndvi)
#
#        # RGB orthomosaic input
#        self.rgb_cb = QComboBox(); self.rgb_cb.setEditable(True)
#        self.rgb_cb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
#        self._populate_layers(self.rgb_cb, raster=True)
#        rgb_btn = QPushButton("…"); rgb_btn.setFixedWidth(20)
#        rgb_btn.clicked.connect(lambda: self._browse(self.rgb_cb, raster=True))
#        h_rgb = QHBoxLayout(); h_rgb.addWidget(self.rgb_cb); h_rgb.addWidget(rgb_btn)
#        form.addRow("RGB orthomosaic:", h_rgb)
#
#        # Coverage ROI polygon
#        self.roi_cb = QComboBox(); self.roi_cb.setEditable(True)
#        self.roi_cb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
#        self.roi_cb.addItem("<none>")
#        self._populate_layers(self.roi_cb, raster=False, geometry=QgsWkbTypes.PolygonGeometry)
#        roi_btn = QPushButton("…"); roi_btn.setFixedWidth(20)
#        roi_btn.clicked.connect(lambda: self._browse(self.roi_cb, raster=False))
#        h_roi = QHBoxLayout(); h_roi.addWidget(self.roi_cb); h_roi.addWidget(roi_btn)
#        form.addRow("Coverage area (polygon):", h_roi)
#
#        # Method: Clip or Extract
#        self.method_clip = QPushButton("Clip to ROI")
#        self.method_clip.setCheckable(True); self.method_clip.setChecked(True)
#        self.method_extract = QPushButton("Extract within ROI")
#        self.method_extract.setCheckable(True)
#        self.method_clip.clicked.connect(lambda: self.method_extract.setChecked(False))
#        self.method_extract.clicked.connect(lambda: self.method_clip.setChecked(False))
#        h_method = QHBoxLayout(); h_method.addWidget(self.method_clip); h_method.addWidget(self.method_extract)
#        form.addRow("Method:", h_method)
#
#        # NDVI threshold
#        self.thr_sb = QDoubleSpinBox(); self.thr_sb.setRange(-1.0, 1.0)
#        self.thr_sb.setDecimals(2); self.thr_sb.setSingleStep(0.01); self.thr_sb.setValue(0.30)
#        form.addRow("NDVI threshold:", self.thr_sb)
#        
#        # Minimum dump area (m²)
#        self.area_sb = QDoubleSpinBox()
#        self.area_sb.setRange(0.0, 100000.0)      # from 0 up to 100 000 m²
#        self.area_sb.setDecimals(2)
#        self.area_sb.setSingleStep(0.25)
#        self.area_sb.setValue(0.25)               # default 0.25 m²
#        form.addRow("Min. dump area (m²):", self.area_sb)
#
#
#        # Output folder
#        self.folder_le = QLineEdit(); self.folder_le.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
#        folder_btn = QPushButton("…"); folder_btn.setFixedWidth(20)
#        folder_btn.clicked.connect(lambda: self._browse_folder(self.folder_le))
#        h_folder = QHBoxLayout(); h_folder.addWidget(self.folder_le); h_folder.addWidget(folder_btn)
#        form.addRow("Output folder:", h_folder)
#
#        # Add results checkbox
#        self.add_to_map = QPushButton("Add results to map")
#        self.add_to_map.setCheckable(True)
#        form.addRow("", self.add_to_map)
#
#        # Help tab
#        help_tab = QWidget(); help_layout = QVBoxLayout(help_tab)
#        tabs.addTab(help_tab, "Help")
#        help_browser = QTextBrowser()
#        help_browser.setHtml(
#            "<h2>Detect Illegal Dumps</h2>"
#            "<p>Select rasters, ROI, method, set threshold, then Run.</p>"
#        )
#        help_layout.addWidget(help_browser)
#
#        # Buttons: Run, Close, Help
#        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Close | QDialogButtonBox.Help)
#        btns.button(QDialogButtonBox.Ok).setText("Run")
#        btns.button(QDialogButtonBox.Close).setText("Close")
#        btns.button(QDialogButtonBox.Help).setText("Help")
#        btns.accepted.connect(self.accept); btns.rejected.connect(self.reject)
#        btns.helpRequested.connect(lambda: tabs.setCurrentWidget(help_tab))
#        main_layout.addWidget(btns)
#
#    def _populate_layers(self, combo, raster=True, geometry=None):
#        for lyr in QgsProject.instance().mapLayers().values():
#            if raster and isinstance(lyr, QgsRasterLayer):
#                combo.addItem(lyr.name(), lyr.id())
#            if not raster and isinstance(lyr, QgsVectorLayer):
#                if geometry is None or lyr.geometryType() == geometry:
#                    combo.addItem(lyr.name(), lyr.id())
#
#    def _browse(self, combo, raster=True):
#        filt = "GeoTIFF (*.tif *.tiff)" if raster else "Vector (*.shp *.geojson)"
#        path, _ = QFileDialog.getOpenFileName(self, "Select file", "", filt)
#        if path: combo.setEditText(path)
#
#    def _browse_folder(self, lineedit):
#        path = QFileDialog.getExistingDirectory(self, "Select folder", "")
#        if path: lineedit.setText(path)
#
#    def values(self):
#        return {
#            'ndvi_src': self.ndvi_cb.currentText(),
#            'rgb_src': self.rgb_cb.currentText(),
#            'roi_src': None if self.roi_cb.currentText() == "<none>" else self.roi_cb.currentText(),
#            'method': 'clip' if self.method_clip.isChecked() else 'extract',
#            'threshold': self.thr_sb.value(),
#            'min_area': self.area_sb.value(),
#            'output_folder': self.folder_le.text(),
#            'add_to_map': self.add_to_map.isChecked()
#        }
#
# -*- coding: utf-8 -*-
"""
IWD QGIS Plugin – Simplified Detect Dialog
Streamlined UI with automatic or manual NDVI threshold
"""

from qgis.PyQt.QtWidgets import (
    QDialog, QTabWidget, QWidget, QVBoxLayout, QFormLayout, QHBoxLayout,
    QLabel, QLineEdit, QComboBox, QPushButton, QDoubleSpinBox,
    QDialogButtonBox, QFileDialog, QTextBrowser, QSizePolicy, QCheckBox
)
from qgis.core import QgsProject, QgsRasterLayer, QgsVectorLayer, QgsWkbTypes

class IWDDialog(QDialog):
    """Parameter dialog: raster inputs, ROI, method, threshold, output."""
    def __init__(self, iface, parent=None):
        super().__init__(parent)
        self.iface = iface
        self.setWindowTitle("Detect illegal dumps")
        self.resize(500, 360)

        main_layout = QVBoxLayout(self)
        tabs = QTabWidget(); main_layout.addWidget(tabs)

        # ---------------------------------------------------------------------
        # Parameters tab
        # ---------------------------------------------------------------------
        param_tab = QWidget()
        form = QFormLayout(param_tab)
        tabs.addTab(param_tab, "Parameters")

        # NDVI index raster
        self.ndvi_cb = QComboBox(); self.ndvi_cb.setEditable(True)
        self.ndvi_cb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._populate_layers(self.ndvi_cb, raster=True)
        ndvi_btn = QPushButton("…"); ndvi_btn.setFixedWidth(20)
        ndvi_btn.clicked.connect(lambda: self._browse(self.ndvi_cb, raster=True))
        h_ndvi = QHBoxLayout(); h_ndvi.addWidget(self.ndvi_cb); h_ndvi.addWidget(ndvi_btn)
        form.addRow("NDVI index raster:", h_ndvi)

        # RGB orthomosaic input
        self.rgb_cb = QComboBox(); self.rgb_cb.setEditable(True)
        self.rgb_cb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._populate_layers(self.rgb_cb, raster=True)
        rgb_btn = QPushButton("…"); rgb_btn.setFixedWidth(20)
        rgb_btn.clicked.connect(lambda: self._browse(self.rgb_cb, raster=True))
        h_rgb = QHBoxLayout(); h_rgb.addWidget(self.rgb_cb); h_rgb.addWidget(rgb_btn)
        form.addRow("RGB orthomosaic:", h_rgb)

        # Coverage ROI polygon or current extent
        # Allow the user to either select a polygon layer/file or use the current map extent
        self.roi_cb = QComboBox()
        self.roi_cb.setEditable(True)
        self.roi_cb.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.roi_cb.addItem("<none>")
        # Populate with existing polygon layers
        self._populate_layers(self.roi_cb, raster=False, geometry=QgsWkbTypes.PolygonGeometry)
        roi_btn = QPushButton("…"); roi_btn.setFixedWidth(20)
        roi_btn.clicked.connect(lambda: self._browse(self.roi_cb, raster=False))
        h_roi = QHBoxLayout(); h_roi.addWidget(self.roi_cb); h_roi.addWidget(roi_btn)
        form.addRow("Coverage area (polygon):", h_roi)

        # Checkbox to use the current map extent as ROI
        self.use_extent_cb = QCheckBox("Use current map extent")
        self.use_extent_cb.setChecked(False)
        # Disable ROI combobox when using current extent
        self.use_extent_cb.toggled.connect(self.roi_cb.setDisabled)
        form.addRow("", self.use_extent_cb)

        # Method: Clip or Extract, use radio buttons for clarity
        from PyQt5.QtWidgets import QRadioButton, QButtonGroup
        self.method_clip = QRadioButton("Clip to ROI")
        self.method_extract = QRadioButton("Extract within ROI")
        # Group the radio buttons to allow only one selection
        self.method_group = QButtonGroup(self)
        self.method_group.addButton(self.method_clip)
        self.method_group.addButton(self.method_extract)
        # Default method is clip
        self.method_clip.setChecked(True)
        h_method = QHBoxLayout(); h_method.addWidget(self.method_clip); h_method.addWidget(self.method_extract)
        form.addRow("Method:", h_method)

        # Automatic NDVI threshold?
        self.auto_cb = QCheckBox("Automatic NDVI threshold (Otsu)")
        self.auto_cb.setChecked(True)
        form.addRow("", self.auto_cb)

        # NDVI threshold (manual)
        self.thr_sb = QDoubleSpinBox()
        self.thr_sb.setRange(-1.0, 1.0)
        self.thr_sb.setDecimals(2)
        self.thr_sb.setSingleStep(0.01)
        self.thr_sb.setValue(0.30)
        self.thr_sb.setEnabled(False)
        # toggle enable when auto_cb toggles
        self.auto_cb.toggled.connect(self.thr_sb.setDisabled)
        form.addRow("NDVI threshold:", self.thr_sb)

        # Minimum dump area (m²)
        self.area_sb = QDoubleSpinBox()
        self.area_sb.setRange(0.0, 100000.0)
        self.area_sb.setDecimals(2)
        self.area_sb.setSingleStep(0.25)
        self.area_sb.setValue(0.25)
        form.addRow("Min. dump area (m²):", self.area_sb)

        # Output folder
        self.folder_le = QLineEdit(); self.folder_le.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        folder_btn = QPushButton("…"); folder_btn.setFixedWidth(20)
        folder_btn.clicked.connect(lambda: self._browse_folder(self.folder_le))
        h_folder = QHBoxLayout(); h_folder.addWidget(self.folder_le); h_folder.addWidget(folder_btn)
        form.addRow("Output folder:", h_folder)

        # Add results checkbox: use QCheckBox instead of a push button
        self.add_to_map = QCheckBox("Add results to map")
        self.add_to_map.setChecked(True)
        form.addRow("", self.add_to_map)


        # ---------------------------------------------------------------------
        # Help tab
        # ---------------------------------------------------------------------
        help_tab = QWidget(); help_layout = QVBoxLayout(help_tab)
        tabs.addTab(help_tab, "Help")
        help_browser = QTextBrowser()
        help_browser.setHtml(
            "<h2>Detect Illegal Dumps</h2>"
            "<p>Select rasters, ROI, choose auto or manual NDVI threshold, "
            "set minimum area, then Run.</p>"
        )
        help_layout.addWidget(help_browser)


        # ---------------------------------------------------------------------
        # Buttons: Run, Close, Help
        # ---------------------------------------------------------------------
        #
        # Button box – expose on the instance for external access
        #
        # Assign the button box to an instance attribute so that tests or
        # other code can retrieve individual buttons (e.g. OK/Cancel).  In
        # earlier versions of this dialog the button box was stored on
        # ``self.button_box``; to preserve API compatibility we reassign it
        # here.
        btns = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Close | QDialogButtonBox.Help
        )
        btns.button(QDialogButtonBox.Ok).setText("Run")
        btns.button(QDialogButtonBox.Close).setText("Close")
        btns.button(QDialogButtonBox.Help).setText("Help")
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        btns.helpRequested.connect(lambda: tabs.setCurrentWidget(help_tab))
        self.button_box = btns  # expose for tests
        main_layout.addWidget(btns)

    def _populate_layers(self, combo, raster=True, geometry=None):
        for lyr in QgsProject.instance().mapLayers().values():
            if raster and isinstance(lyr, QgsRasterLayer):
                combo.addItem(lyr.name(), lyr.id())
            if not raster and isinstance(lyr, QgsVectorLayer):
                if geometry is None or lyr.geometryType() == geometry:
                    combo.addItem(lyr.name(), lyr.id())

    def _browse(self, combo, raster=True):
        filt = "GeoTIFF (*.tif *.tiff)" if raster else "Vector (*.shp *.geojson)"
        path, _ = QFileDialog.getOpenFileName(self, "Select file", "", filt)
        if path:
            combo.setEditText(path)

    def _browse_folder(self, lineedit):
        path = QFileDialog.getExistingDirectory(self, "Select folder", "")
        if path:
            lineedit.setText(path)

    def values(self):
        """Return a dictionary with all user-specified parameters.

        - ndvi_src/rgb_src: either a loaded layer name or a file path
        - roi_src: name/path of a polygon layer if provided
        - use_extent: whether to use the current map extent instead of ROI polygon
        - method: 'clip' or 'extract' for ROI clipping/selection
        - ndvi_threshold: None when automatic threshold is requested
        - min_area: minimum area for candidate polygons
        - output_folder: folder where results (shapefile, PDFs) will be saved
        - add_to_map: whether to add results to the map canvas automatically
        """
        roi_text = self.roi_cb.currentText()
        return {
            'ndvi_src': self.ndvi_cb.currentText().strip(),
            'rgb_src': self.rgb_cb.currentText().strip(),
            'roi_src': None if roi_text == "<none>" or self.use_extent_cb.isChecked() else roi_text,
            'use_extent': self.use_extent_cb.isChecked(),
            'method': 'clip' if self.method_clip.isChecked() else 'extract',
            'ndvi_threshold': None if self.auto_cb.isChecked() else self.thr_sb.value(),
            'min_area': self.area_sb.value(),
            'output_folder': self.folder_le.text().strip(),
            'add_to_map': self.add_to_map.isChecked()
        }
