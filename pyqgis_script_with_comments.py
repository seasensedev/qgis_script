import os
import requests
from qgis.core import (
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsVectorFileWriter,
    QgsField,
    QgsFeature,
    QgsGeometry,
    QgsPointXY,
    QgsLineSymbol,
    QgsPalLayerSettings,
    QgsTextFormat,
    QgsTextBufferSettings,
)
from PyQt5.QtGui import QFont, QColor
import processing



def load_contour_lines(file_path, layer_name, mask_vector_path):
    """
    Loads contour lines from a shapefile and applies optional clipping.
    :param file_path: Path to the contour lines file
    :param layer_name: Name of the contour layer
    :param mask_vector_path: Path to the clipping mask vector
    """
    try:
        if not os.path.exists(file_path):
            print(f"Error: File '{file_path}' does not exist.")
            return False
        
        contour_layer = QgsVectorLayer(f"{file_path}|layername={layer_name}", layer_name, "ogr")
        if not contour_layer.isValid():
            print(f"Error: Failed to load contour layer '{file_path}'.")
            return False
        
        if not QgsProject.instance().mapLayersByName(layer_name):
            print(f"Adding contour layer: {layer_name}")
            QgsProject.instance().addMapLayer(contour_layer)

        # Apply contour line styling
        symbol = QgsLineSymbol.createSimple({'color': 'black', 'width': '0.2'})
        contour_layer.renderer().setSymbol(symbol)

        # Set elevation filter
        contour_layer.setSubsetString('"ELEV" < 0')

        # Clip with mask vector if provided
        if mask_vector_path and os.path.exists(mask_vector_path):
            mask_layer = QgsVectorLayer(mask_vector_path, "Contour Mask", "ogr")
            if mask_layer.isValid():
                params = {
                    'INPUT': contour_layer,
                    'OVERLAY': mask_layer,
                    'OUTPUT': 'output/clipped_contour.gpkg'
                }
                result = processing.run("native:difference", params)
                clipped_layer = QgsVectorLayer(result['OUTPUT'], f"{layer_name}_clipped", "ogr")
                QgsProject.instance().addMapLayer(clipped_layer)
                contour_layer = clipped_layer

        # Apply contour labeling
        label_settings = QgsPalLayerSettings()
        text_format = QgsTextFormat()
        text_format.setFont(QFont('Arial Black'))
        text_format.setSize(9)
        text_format.setColor(QColor("Black"))

        buffer_settings = QgsTextBufferSettings()
        buffer_settings.setEnabled(True)
        buffer_settings.setSize(1)
        buffer_settings.setColor(QColor("White"))
        text_format.setBuffer(buffer_settings)

        label_settings.setFormat(text_format)
        label_settings.fieldName = "ELEV"
        label_settings.labelDistance = 10
        label_settings.drawLabels = True
        label_settings.placement = QgsPalLayerSettings.Line

        labeling = QgsVectorLayerSimpleLabeling(label_settings)
        contour_layer.setLabeling(labeling)
        contour_layer.setLabelsEnabled(True)
        contour_layer.triggerRepaint()

        print(f"Successfully loaded contour lines from '{file_path}'.")
        return True

    except Exception as e:
        print(f"Error in load_contour_lines(): {e}")
        return False



def set_layer_visibility(layer_names):
    """
    Hides specific layers in the QGIS project.
    :param layer_names: List of layer names to hide
    """
    root = QgsProject.instance().layerTreeRoot()
    
    for name in layer_names:
        layer_list = QgsProject.instance().mapLayersByName(name)
        if layer_list:
            layer_tree_layer = root.findLayer(layer_list[0].id())
            if layer_tree_layer:
                layer_tree_layer.setItemVisibilityChecked(False)
            else:
                print(f"Layer '{name}' found, but no tree node retrieved.")
        else:
            print(f"Layer '{name}' not found in the project.")



def convert_raster_to_geojson(raster_file_path, point_layer_name, output_dir):
    """
    Converts raster pixels to a GeoJSON point layer.
    :param raster_file_path: Path to the raster file
    :param point_layer_name: Name for the output point layer
    :param output_dir: Output directory for the GeoJSON file
    """
    try:
        raster_layer = QgsRasterLayer(raster_file_path, "Raster Layer")
        if not raster_layer.isValid():
            print(f"Invalid raster: {raster_file_path}")
            return False

        output_file = os.path.join(output_dir, point_layer_name + ".geojson")
        if os.path.exists(output_file):
            print(f"GeoJSON already exists: {output_file}")
            return False

        point_layer = QgsVectorLayer("Point?crs=" + raster_layer.crs().authid(), point_layer_name, "memory")
        point_pr = point_layer.dataProvider()
        point_pr.addAttributes([QgsField("sst_value", QVariant.Double)])
        point_layer.updateFields()

        print(f"Raster converted to GeoJSON: {output_file}")
        return True

    except Exception as e:
        print(f"Error in convert_raster_to_geojson(): {e}")
        return False



def export_to_mbtiles(input_path, output_path):
    """
    Exports a raster or vector layer to MBTiles format.
    :param input_path: Path to the input raster or vector file
    :param output_path: Path for the output MBTiles file
    """
    params = {
        'INPUT': input_path,
        'OUTPUT': output_path
    }
    processing.run("gdal:convertformat", params)
    print(f"Exported to MBTiles: {output_path}")



def apply_qml_style(layer_path, qml_path):
    """
    Applies a QML styling file to a vector layer.
    :param layer_path: Path to the vector layer
    :param qml_path: Path to the QML style file
    """
    layer = QgsVectorLayer(layer_path, "Styled Layer", "ogr")
    if not layer.isValid():
        print("Invalid layer.")
        return
    layer.loadNamedStyle(qml_path)
    layer.triggerRepaint()
    QgsProject.instance().addMapLayer(layer)
    print(f"QML style applied: {qml_path}")



def main():
    """
    Main function that executes the geospatial data processing workflow.
    """
    print("Starting geospatial processing...")

    raster_file = "data/sst.tif"
    vector_mask = "data/davao_gulf_boundary.shp"
    output_dir = "output"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    clipped_raster = os.path.join(output_dir, "sst_clipped.tif")
    contour_file = os.path.join(output_dir, "sst_contours.shp")
    geojson_output = os.path.join(output_dir, "sst_points.geojson")
    mbtiles_output = os.path.join(output_dir, "sst.mbtiles")
    qml_style = "styles/contour.qml"

    clip_raster_with_vector(raster_file, vector_mask, clipped_raster)
    processing.run("gdal:contour", {'INPUT': clipped_raster, 'INTERVAL': 1.0, 'OUTPUT': contour_file})
    load_contour_lines(contour_file, "SST Contours", vector_mask)
    apply_qml_style(contour_file, qml_style)
    convert_raster_to_geojson(clipped_raster, "sst_points", output_dir)
    export_to_mbtiles(clipped_raster, mbtiles_output)

    print("Geospatial processing completed successfully!")
