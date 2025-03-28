from qgis.core import (QgsProject, QgsRasterLayer, QgsRaster, QgsStyle, QgsColorRampShader, QgsSingleBandPseudoColorRenderer,
    QgsRasterShader, QgsVectorLayer, QgsPalLayerSettings, QgsTextFormat, QgsLineSymbol, 
    QgsVectorLayerSimpleLabeling, QgsApplication, QgsPointXY, QgsFeature, QgsGeometry,
    QgsField, QgsFields, QgsVectorFileWriter, QgsWkbTypes, QgsProcessingParameters,
    QgsCoordinateReferenceSystem, QgsTextBufferSettings,
)
from qgis.analysis import QgsNativeAlgorithms, QgsRasterCalculator, QgsRasterCalculatorEntry
from PyQt5.QtGui import QColor, QFont
from PyQt5.QtCore import QVariant

import sys
import os
import subprocess
import requests

QgsApplication.setPrefixPath("D:\PROGRAM_FILES\QGIS\apps\qgis-ltr", True)
qgs = QgsApplication([], False)
qgs.initQgis()

sys.path.append(r'D:\PROGRAM_FILES\QGIS\apps\qgis-ltr\python\plugins') # Folder where Processing is located
from processing.core.Processing import Processing
Processing.initialize()
import processing

#variables
clippedraster_qml_style_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\styles\clipped_raster\v2_clipped_raster_styling.qml"
mbtiles_output_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\mbtile\sst_davaogulf_maptile.mbtiles"
sst_geojson = "sst_davaogulf_points"
# clipped_rasterName = "D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\clipped\clipped_sst_davao_gulf.tif"
output_geojson_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\geojson\sst_davaogulf_points.geojson"
output_geojson_dir = "D:/PROGRAM_FILES/QGIS_WORKS/SST_AUG_5_2024/MAIN/geojson/"
layers_to_hide = ['sst_davao_gulf_celsius', 'sst_davao_gulf_rasterized', 'davao_gulf_vector', 'sst_davaogulf_points']
# contour_line_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\contour\ng_elevation_contours.gpkg" #Keep
v2_contourline_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\contour\v2_davaogulf_contour.gpkg"
v2_contour_qml_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\contour\v2_davaogulf_contour_style.qml"
contour_layerName = "ng_elevation_contours" #Keep
contour_lines_mbtile = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\contour\sst_davaogulf_contour_lines.gpkg"
ex_davao_gulf_rst = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\sst_davao_gulf.tif"
rasterL_name = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\sst_davao_gulf_celsius.tif"
vectorL_name = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\vector\davao_gulf_except_coverage_for_samal.gpkg"
output_clippedraster_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\clipped\clipped_sst_davao_gulf.tif"
layers_to_mbtiles = [output_clippedraster_path, contour_lines_mbtile]
files_to_delete = [ex_davao_gulf_rst, rasterL_name, output_clippedraster_path, output_geojson_path, mbtiles_output_path]
#variables



def inspect_nc_dataset():
    netcdf_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\dataset\20241020090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc"
    
    try:
        result = subprocess.run(['gdalinfo', netcdf_path], capture_output=True, text=True, check=True)
        
        print(result.stdout)
        
    except subprocess.CalledProcessError as e:
        print(f"Error inspecting dataset: {e}")


def load_nc_dataset():
    netcdf_directory = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\dataset"
    outp_raster_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\sst_davao_gulf.tif"
    
    if not os.path.exists(netcdf_directory):
        print(f"Dataset directory does not exist. Aborting load_nc_dataset...\n")
        return False
        
    netcdf_files = [os.path.join(netcdf_directory, f) for f in os.listdir(netcdf_directory) if f.endswith(".nc")]
    
    if not netcdf_files:
        print(f"No NetCDF files detected in this directory ({netcdf_directory}). Aborting load_nc_dataset...\n")
        return False
    
    latest_file  = max(netcdf_files, key=os.path.getctime)
    print(f"Latest dataset detected: {latest_file}")
    
    if os.path.exists(latest_file):
        variable = 'NETCDF:"{}":analysed_sst'.format(latest_file)
        params = {
            'INPUT': variable,
            'BAND': 1,
            'OUTPUT': outp_raster_path,
            'OUTPUT_TYPE': 5
            }
                
        try:
            processing.run("gdal:translate", params)
            print('GeoTIFF added successfully!')
        except Exception as e:
            print(f"Error creating GeoTIFF: {e}\n")
            return False
    else:
        print('File does not exist!')


def convert_raster_to_celsius(ex_rasterLayer, retries=1):
    # Check if the raster file exists
    if not os.path.exists(ex_rasterLayer):
        print("Raster file does not exist. Loading layer... Executing load_nc_dataset()!")
        
        if retries > 0:
            load_nc_dataset()  # Ensure this method loads the necessary raster file
            print(f"Restarting convert_raster_to_celsius after loading the dataset...\n")
            return convert_raster_to_celsius(ex_rasterLayer, retries - 1)
        else:
            print("SST raster layer still does not exist. Aborting convert_raster_to_celsius method execution!")
            return False
    
    print(f"Raster file '{ex_rasterLayer}' found. Proceeding with conversion...\n")
    
    # Load the raster layer from the file path
    raster_layer = QgsRasterLayer(ex_rasterLayer, "sst_davao_gulf_rasterized")
    
    if not raster_layer.isValid():
        print(f"Failed to load raster layer from '{ex_rasterLayer}'. Aborting conversion!")
        return False

    output_raster_path = r"D:\PROGRAM_FILES\QGIS_WORKS\SST_AUG_5_2024\MAIN\raster\sst_davao_gulf_celsius.tif"
    
    # Optionally remove the old output file
    if os.path.exists(output_raster_path):
        os.remove(output_raster_path)

    entry = QgsRasterCalculatorEntry()
    entry.ref = 'raster@1'
    entry.raster = raster_layer
    entry.bandNumber = 1
    
    expression = 'raster@1 - 273.15'
    calc = QgsRasterCalculator(expression, output_raster_path, 'GTiff', raster_layer.extent(), raster_layer.width(), raster_layer.height(), [entry])
    
    if calc.processCalculation() == 0:
        print(f"Converted raster saved successfully at '{output_raster_path}'!\n")
        return output_raster_path  # Return the path of the output raster
    else:
        print('Error during raster calculation!')
        return False



def clip_raster_with_vector(raster_file_path, vector_file_path, output_path):
    # Check if clipped raster file already exists
    if os.path.exists(output_path):
        print(f"Clipped raster file '{output_path}' already exists. Skipping clip_raster_with_vector() method execution.\n")
        return False
    
    # Load the raster layer from the file path
    raster_layer = QgsRasterLayer(raster_file_path, "sst_davao_gulf_rasterized")
    if not raster_layer.isValid():
        print(f"1. Raster layer does not exist! Check the file path: '{raster_file_path}'\n")
        return False
    
    # Set the CRS for the raster layer (if needed)
    raster_layer.setCrs(QgsCoordinateReferenceSystem("EPSG:4326"))
    
    # Load the vector layer from the file path
    vector_layer = QgsVectorLayer(vector_file_path, "davao_gulf_vector", "ogr")
    if not vector_layer.isValid():
        print(f"1. Vector layer does not exist! Check the file path: '{vector_file_path}'\n")
        return False
    
    # Remove the output file if it exists
    if os.path.exists(output_path):
        os.remove(output_path)
        
    # Mask layer clipping parameters
    params = {
        'INPUT': raster_layer.source(),
        'MASK': vector_layer.source(),
        'SOURCE_CRS': vector_layer.crs(),
        'TARGET_CRS': vector_layer.crs(),
        'OUTPUT': output_path
    }
    
    try:
        processing.run("gdal:cliprasterbymasklayer", params)
        # Load the clipped layer but do not add it to the project yet
        clipped_layer = QgsRasterLayer(output_path, "clipped_sst_davaogulf_celcius")
        if clipped_layer.isValid():
            print(f"Successfully created the clipped raster layer: {clipped_layer}\n")
            # Return the output path for later use
            return output_path
        else:
            print(f"2. Failed to load the clipped raster layer: !!{clipped_layer}\n")
            return False
    except Exception as e:
        print(f"E. Error during clipping process: {e}\n")
        return False
        


def apply_singleband_pseudocolor(raster_file_path, color_ramp='Turbo', min_value=29, max_value=32, classes=8):
    try:
        # Load raster layer from file path
        raster_layer = QgsRasterLayer(raster_file_path, "SST Raster Layer")
        if not raster_layer.isValid():
            print(f"Layer from file path '{raster_file_path}' is not valid. Aborting apply_singleband_pseudocolor() execution.\n")
            return False
        
        # Set the CRS if needed (modify according to your requirements)
        raster_layer.setCrs(QgsCoordinateReferenceSystem("EPSG:4326"))

        # Check if the color ramp exists
        ramp = QgsStyle().defaultStyle().colorRamp(color_ramp)
        if ramp is None:
            print(f"Error: Color ramp '{color_ramp}' does not exist or is invalid. Aborting apply_singleband_pseudocolor() method execution.\n")
            return False

        color_ramp_shader = QgsColorRampShader(min_value, max_value, ramp, QgsColorRampShader.Interpolated, QgsColorRampShader.Quantile)

        # Specified values for color ramp items
        specified_values = [29.1840153, 29.7679239, 29.781001, 29.8110074, 29.9179559, 30.2429578, 30.4368844, float('inf')]
        
        color_ramp_items = []
        for i, value in enumerate(specified_values):
            color = ramp.color(float(i) / (classes - 1))
            label = f"<= {value:.2f}" if i == 0 else f"> {value:.2f}" if i == len(specified_values) - 1 else f"{value:.2f}"
            color_ramp_items.append(QgsColorRampShader.ColorRampItem(value, color, label))

        color_ramp_shader.setColorRampItemList(color_ramp_items)

        print("Color Ramp Items Generated:")
        for item in color_ramp_items:
            print(f"Value: {item.value}, Color: {item.color.name()}, Label: {item.label}")

        shader = QgsRasterShader()
        shader.setRasterShaderFunction(color_ramp_shader)

        # Create a pseudo color renderer
        renderer = QgsSingleBandPseudoColorRenderer(raster_layer.dataProvider(), raster_layer.type(), shader)

        # Set the renderer and trigger repaint
        raster_layer.setRenderer(renderer)
        raster_layer.triggerRepaint()

        # Optionally add the raster layer to the project if needed
        QgsProject.instance().addMapLayer(raster_layer)

        print(f"Color ramp '{color_ramp}' applied to layer from file path '{raster_file_path}' with min value of {min_value} and max value of {max_value}\n")
        return True

    except Exception as e:
        print(f"An error occurred: {e}\n")
        return False

def clippedraster_apply_qml_style(raster_filepath, qml_file_path):
    try:
        raster_layer = QgsRasterLayer(raster_filepath, "clipped_sst_davaogulf")
        if not raster_layer:
            print(f"Raster layer from file path; {raster_filepath} does not exist! Aborting clippedraster_apply_qml_style...\n")
            return False
            
        raster_layer.setCrs(QgsCoordinateReferenceSystem("EPSG:4326"))
        
        if not raster_layer.loadNamedStyle(qml_file_path):
            print(f"failed to apply qml styling on clipped raster layer. Aborting clippedraster_apply_qml_style()...\n")
            return False
        
        raster_layer.triggerRepaint()
        
        QgsProject.instance().addMapLayer(raster_layer)
        
        print(f"QML style successfully applied to clipped raster layer.\n")
        return True
        
    except Exception as e:
        print(f"An error occured at")

def load_contour_lines(file_path, layer_name, mask_vector_path):
    try:
        # Check if the contour lines file exists
        if not os.path.exists(file_path):
            print(f"Error: The file '{file_path}' does not exist.\n")
            return False
        
        # Load the contour layer from the file path
        contour_layer = QgsVectorLayer(f"{file_path}|layername={layer_name}", layer_name, "ogr")
        
        # Check if the contour layer is valid
        if not contour_layer.isValid():
            print(f"Error: Failed to load contour layer '{file_path}'.\n")
            return False
        
        # Add the contour layer to the project if it does not already exist
        existing_layer = QgsProject.instance().mapLayersByName(layer_name)
        if not existing_layer:
            print(f"Adding contour line {layer_name}\n")
            QgsProject.instance().addMapLayer(contour_layer)
        
        # Set the symbol for the contour lines
        symbol = QgsLineSymbol.createSimple({'color': 'black', 'width': '0.2'})
        contour_layer.renderer().setSymbol(symbol)
        
        # Set the subset string to filter elevation
        contour_layer.setSubsetString('"ELEV" < 0')
        
        # Check and load the mask vector layer if the path is provided
        if mask_vector_path and os.path.exists(mask_vector_path):
            mask_layer = QgsVectorLayer(mask_vector_path, "davaogulf_contour_mask", "ogr")
            if mask_layer.isValid():
                #QgsProject.instance().addMapLayer(mask_layer)
                
                # Clip the contour lines with the mask
                params = {
                    'INPUT': contour_layer,
                    'OVERLAY': mask_layer,
                    'OUTPUT': 'D:\\PROGRAM_FILES\\QGIS_WORKS\\SST_AUG_5_2024\\MAIN\\raster\\contour\\sst_davaogulf_contour_lines.gpkg'
                }
                
                result = processing.run("native:difference", params)
                clipped_layer = QgsVectorLayer(result['OUTPUT'], f"{layer_name}_clipped", "ogr")
                QgsProject.instance().addMapLayer(clipped_layer)
                contour_layer = clipped_layer
        
        # Contour labeling
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
        
        # Trigger repaint to update the layer display
        contour_layer.triggerRepaint()
        
        print(f"Successfully loaded the contour lines from '{file_path}'.\n")
        return True
        
    except Exception as e:
        print(f"Error executing load_contour_lines(): {e}\n")
        return False



def set_layer_visibility(layer_names):
    root = QgsProject.instance().layerTreeRoot()
    
    for name in layer_names:
        # Get the layer by name
        layer_list = QgsProject.instance().mapLayersByName(name)
        
        if layer_list:
            layer_tree_layer = root.findLayer(layer_list[0].id())
            if layer_tree_layer:
                # Hide the layer by setting visibility to False
                layer_tree_layer.setItemVisibilityChecked(False)
            else:
                print(f"Layer '{name}' found, but could not retrieve layer tree node.\n")
        else:
            print(f"Layer '{name}' not found in the project.\n")



def convert_raster_to_geojson(raster_file_path, point_layer_name, output_dir):
    try:
        # Load raster layer from file path
        tgt_layer = QgsRasterLayer(raster_file_path, "Raster Layer")
        if not tgt_layer.isValid():
            print(f"Raster layer '{raster_file_path}' is invalid. Aborting conversion!\n")
            return False

        # Check if point layer already exists in the specified output location
        output_file = os.path.join(output_dir, point_layer_name + ".geojson")
        if os.path.exists(output_file):
            print(f"GeoJSON file '{output_file}' already exists. Skipping conversion.\n")
            return False

        # Create an in-memory point layer
        point_layer = QgsVectorLayer("Point?crs=" + tgt_layer.crs().authid(), point_layer_name, "memory")
        point_pr = point_layer.dataProvider()
        
        # Add SST value attribute to the point layer
        point_pr.addAttributes([QgsField("sst_value", QVariant.Double)])
        point_layer.updateFields()

        # Retrieve raster properties
        raster = tgt_layer.dataProvider()
        extent = tgt_layer.extent()
        width = tgt_layer.width()
        height = tgt_layer.height()
        
        # Convert each pixel to a point if it has a valid value
        for row in range(height):
            for col in range(width):
                pixel_value = raster.identify(
                    QgsPointXY(extent.xMinimum() + col * extent.width() / width,
                               extent.yMaximum() - row * extent.height() / height),
                    QgsRaster.IdentifyFormatValue
                ).results()

                if pixel_value:
                    value = pixel_value[1]
                    if value is not None:
                        x = extent.xMinimum() + col * extent.width() / width
                        y = extent.yMaximum() - row * extent.height() / height
                        point_feature = QgsFeature()
                        point_feature.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
                        point_feature.setAttributes([value])
                        point_pr.addFeatures([point_feature])

        # Save the point layer to GeoJSON format without adding it to the map
        QgsVectorFileWriter.writeAsVectorFormat(
            point_layer, output_file, "utf-8", tgt_layer.crs(), "GeoJSON"
        )
        
        print(f"Raster '{raster_file_path}' successfully converted to GeoJSON as '{output_file}'.\n")
        return True

    except Exception as e:
        print(f"An error occurred during the conversion process: {e}")
        return False
    


def convert_to_mbtiles(layer_paths, output_mbtiles_path):
    try:
        layers = []
        
        # Load each layer from file path
        for path in layer_paths:
            if path.lower().endswith('.tif'):
                layer = QgsRasterLayer(path, os.path.basename(path))
            elif path.lower().endswith('.gpkg'):
                layer = QgsVectorLayer(path, os.path.basename(path), "ogr")
            else:
                print(f"Unsupported file type for '{path}'. Skipping.\n")
                continue
            
            if layer.isValid():
                layers.append(layer)
            else:
                print(f"Layer at '{path}' is invalid or could not be loaded.\n")
        
        if not layers:
            print("No valid layers found. Aborting MBTiles generation.\n")
            return False
        
        # Set extent from the first layer's extent
        extent = layers[0].extent()
        
        # Define parameters for MBTiles export
        params = {
            'LAYERS': [layer.id() for layer in layers],
            'ZOOM_MIN': 9,
            'ZOOM_MAX': 12,
            'DPI': 96,
            'TILE_FORMAT': 0,  # 0 for PNG, 1 for JPEG
            'TMS_CONVENTION': False,
            'EXTENT': extent,
            'METATILESIZE': 4,
            'OUTPUT_FILE': output_mbtiles_path,
            'COMPRESS': True
        }
        
        # Execute MBTiles creation
        processing.run("native:tilesxyzmbtiles", params)
        
        print(f"Successfully converted map layers to MBTiles at '{output_mbtiles_path}'.\n")
        return True

    except Exception as e:
        print(f"Error converting map layers to MBTiles: {e}\n")
        return False
        

def apply_qml_style_to_contour(contour_file_path, qml_style_file_path):
    try:
        contour_layer = QgsVectorLayer(contour_file_path, "Contour Layer", "ogr")
        if not contour_layer.isValid():
            print("Failed to load the contour layer. Aborting apply_qml_style_to_contour()...\n")
            return False
            
        
        if os.path.exists(qml_style_file_path):
            contour_layer.loadNamedStyle(qml_style_file_path)
            contour_layer.triggerRepaint()
            print(f"QML style '{qml_style_file_path}' applied to contour layer.\n")
            QgsProject.instance().addMapLayer(contour_layer)
            
        else:
            print(f"QML style '{qml_style_file_path}' does not exist.\n")
            return False
            
            
        return True
    except Exception as e:
        print(f"Error applying QML style at apply_qml_style_to_contour(): {e}\n")
        return False


    

def delete_files(files):
    for file_path in files:
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path}\n")
            except Exception as e:
                print(f"Error deleting: {file_path}\n {e}")
        else:
            print(f"File does not exist: {file_path}")

def load_basemap():
    # Google Satellite TMS layer
    google_maps_layer_name = "Google Satellite"
    
    # Check if Google Maps layer already exists in the project
    if QgsProject.instance().mapLayersByName(google_maps_layer_name):
        print(f"Layer '{google_maps_layer_name}' already exists. Skipping Google Maps satellite layer loading...\n")
    else:
        # Define the URL for Google Satellite
        url = "mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}"
        google_maps = f"type=xyz&zmin=0&zmax=21&url=https://{requests.utils.quote(url)}"
        
        # Create the TMS layer
        tms_layer = QgsRasterLayer(google_maps, google_maps_layer_name, "wms")
        
        if tms_layer.isValid():
            QgsProject.instance().addMapLayer(tms_layer)
            print(f"Google satellite layer '{google_maps_layer_name}' loaded successfully!\n")
        else:
            print(f"Failed to load Google satellite layer!\n")

    # Davao Gulf Vector Layer
    davao_gulf_vector_layer_name = "davao_gulf_vector"
    
    # Check if Davao Gulf vector layer already exists
    if QgsProject.instance().mapLayersByName(davao_gulf_vector_layer_name):
        print(f"Layer '{davao_gulf_vector_layer_name}' already exists. Skipping Davao Gulf vector layer loading...\n")
    else:
        # Load vector layer from file path
        vector_layer_path = r"D:/PROGRAM_FILES/QGIS_WORKS/SST_AUG_5_2024/MAIN/vector/davao_gulf_except_coverage_for_samal.gpkg"
        davao_gulf_vector = QgsVectorLayer(vector_layer_path, davao_gulf_vector_layer_name, "ogr")
        
        if davao_gulf_vector.isValid():
            # QgsProject.instance().addMapLayer(davao_gulf_vector)
            print(f"Vector layer '{davao_gulf_vector_layer_name}' added successfully!\n")
        else:
            print(f"Failed to load vector layer for Davao Gulf!\n")


# inspect_nc_dataset()
def main():
    delete_files(files_to_delete)
    load_nc_dataset()
    load_basemap()
    if not convert_raster_to_celsius(ex_davao_gulf_rst):
        print("convert_raster_to_celsius() method terminated!\n")
    else:
        print("Layer already exist \n")
    
    apply_qml_style_to_contour(v2_contourline_path, v2_contour_qml_path)
    clip_raster_with_vector(rasterL_name, vectorL_name, output_clippedraster_path)
    # apply_singleband_pseudocolor(output_clippedraster_path, color_ramp='Turbo', min_value=29.1840153, max_value=30.7619877, classes=8)
    # load_contour_lines(contour_line_path, contour_layerName, vectorL_name)
    clippedraster_apply_qml_style(output_clippedraster_path, clippedraster_qml_style_path)
    convert_raster_to_geojson(output_clippedraster_path, sst_geojson, output_geojson_dir)
    set_layer_visibility(layers_to_hide)
    convert_to_mbtiles(layers_to_mbtiles, mbtiles_output_path)
    
    # Refresh map
    if 'iface' in globals():
        iface.mapCanvas().refresh()
        

if __name__ == "__main__":
    main()

    qgs.exitQgis()

# Processing tasks:
# 1. load_basemap()
# 2.1. convert_raster_to_celsius()
# 2.2. convert_raster_to_celsius() -> load_nc_dataset() -> Repeat {This only happens when the dataset raster layer is not loaded yet.}
# 3. clipped_raster_with_vector(rasterLayerName, vectorLayerName, outputPath)