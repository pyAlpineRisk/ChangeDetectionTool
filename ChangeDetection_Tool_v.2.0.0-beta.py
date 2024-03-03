# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   QGIS ChangeDetection Tool v.2.0.0-beta                                *
#   QQGIS-Version: 3.28.15-Firenze                                        * 
*   pyAlpineRisk                                                          *
*   Nicole Kamp                                                           *
*   niki.kamp@gmail.com                                                   *
*   www.nicolekamp.comb                                                   *
*   October 2021                                                          *
*   revised 03/2024                                                       *
*                                                                         *
*   Bugs:                                                                 *
*        1. Visualization does not match with legend                      *
*        2. Error with group visualization                                *
***************************************************************************
"""


from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProject,
                       QgsVectorLayer,
                       QgsTextFormat,
                       QgsExpression,
                       QgsFeatureRequest,
                       QgsFeature,
                       QgsGeometry,
                       QgsPoint,
                       QgsPointXY,
                       QgsVectorFileWriter,
                       QgsRasterBandStats,
                       QgsColorRampShader,
                       QgsRasterTransparency,
                       QgsFillSymbol,
                       QgsRasterShader,
                       QgsSingleBandPseudoColorRenderer,
                       QgsWkbTypes,
                       QgsVectorLayerSimpleLabeling,
                       QgsPalLayerSettings,
                       QgsFeatureSink,
                       QgsProcessingException,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterFile,
                       QgsRasterLayer,
                       QgsCoordinateReferenceSystem,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterExtent,
                       QgsProcessingMultiStepFeedback,
                       QgsProcessingParameterEnum,
                       QgsMessageLog)
from qgis import processing
from qgis.analysis import QgsRasterCalculatorEntry, QgsRasterCalculator
from qgis.PyQt.QtWidgets import QApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtGui import QColor
from qgis.utils import iface

import osgeo.gdal as gdal
import numpy as np
import csv

from osgeo import ogr, osr
from shapely.geometry import Polygon
from matplotlib import rcParams
import matplotlib.pyplot as plt 
import pandas as pd
import subprocess
from subprocess import call

import string, os, sys, copy, shutil, math, numpy, time, datetime
from time import *
from sys import *

## ----------------------------------------------------------------------------------------------------------##
## ----------------------------------------------------------------------------------------------------------##
# FUNCTIONS
## ----------------------------------------------------------------------------------------------------------##
# Raster to array
def raster_to_array(raster):
    in_rds = gdal.Open(raster)
    rds_band = in_rds.GetRasterBand(1)
    rds_array = rds_band.ReadAsArray()
    return rds_array

## ----------------------------------------------------------------------------------------------------------##
#Analyse thresholds
def process_threshold(threshold, diff_array, uc_array, cs, width, height, countAOI, analyse, analyse_csv, final_path, gt, epsg):
    max_diff = 10000
    min_uc = 0
    max_uc = 1
    volErrorWithinUC = 0
        
    if cs == 0:
        raise ValueError("Zellengröße (cs) darf nicht 0 sein.")
        
    countPosGes = 0
    countNegGes = 0
    countPos = 0
    countNeg = 0
    countPosArea = 0
    countNegArea = 0
    countAcc = 0
    countCell = 0
    countCellVol = 0
    
    # Create difference raster
    # Save thresholded difference raster
    if threshold != "UC":  
        threshold_value = float(threshold)  
        # Create a corrected difference array based on the threshold value
        corrected_diff = np.where(np.abs(diff_array) <= threshold_value, 0, diff_array)
    if threshold == "UC":  
        corrected_diff = np.where((np.abs(diff_array) <= uc_array), 0, diff_array)
    # Masks for calculating volumes and areas
    pos_mask = corrected_diff > 0
    neg_mask = corrected_diff < 0
    volPos = np.sum(corrected_diff[pos_mask]) * cs * cs
    volNeg = np.sum(np.abs(corrected_diff[neg_mask])) * cs * cs
    # Write difference raster
    threshold_label = 'uc' if threshold == "UC" else ('raw' if threshold == 0.0 else f'{threshold_value}m')
    output_path = os.path.join(final_path, f'corrected_diff_{threshold_label}.tif')
    save_grid_as_raster(output_path, corrected_diff, gt, epsg)
    
    # Detailed analysis
    for row in range(height):
        for col in range(width):
            diff = diff_array[row, col]
            diff = min(max(diff, -max_diff), max_diff) 
            uc = uc_array[row, col] if uc_array is not None and threshold == "UC" else 0
            
            if np.isinf(diff):
                continue
                            
            if threshold == "UC":
                ES = np.clip(uc, min_uc, max_uc) 
                if abs(diff) <= ES:
                    countCell += 1
                    cell_vol = cs * cs * abs(diff)
                    countCellVol += cell_vol
                    countAcc += cell_vol
                elif diff > ES:
                    volPos = cs * cs * (diff - ES) #neu
                    pos_value = cs * cs * (diff - ES)
                    countPos += pos_value
                    countPosArea += cs * cs
                elif diff < -ES:
                    volNeg = cs * cs * (abs(diff) - ES)
                    neg_value = cs * cs * (abs(diff) - ES)
                    countNeg += neg_value
                    countNegArea += cs * cs
            else:
                threshold_value = float(threshold)
                if -threshold_value < diff < threshold_value:
                    countCell += 1
                    cell_vol = cs * cs * abs(diff)
                    countCellVol += cell_vol
                    countAcc += cell_vol
                if abs(diff) > threshold_value:
                    if diff < -threshold_value:
                        volNeg = cs * cs * (abs(diff) - threshold_value)
                        countNeg += volNeg
                        countNegArea += cs * cs
                        countNegGes += cs * cs * abs(diff)
                    elif diff > threshold_value:
                        volPos = cs * cs * (diff - threshold_value)
                        countPos += volPos
                        countPosArea += cs * cs
                        countPosGes += cs * cs * diff

    # Write reults for threshold to TXT/CSV
    write_analysis(analyse, analyse_csv, threshold, countAOI, countPosArea, countNegArea, countPos, countNeg, countAcc, countCell, countCellVol, countPosGes, countNegGes)
        
    return {'volPos': round(countPos, 2), 'volNeg': round(countNeg, 2)}

## ----------------------------------------------------------------------------------------------------------##
# Write results to TXT/CSV
def write_analysis(analyse, analyse_csv, threshold, countAOI, countPosArea, countNegArea, countPos, countNeg, countAcc, countCell, countCellVol, countPosGes, countNegGes):
    analyse.write(f"Analysis with Threshold of +/- {threshold}\n")
    analyse.write(f"Thresholded Area [m2] of Detectable Change: {round(countPosArea+countNegArea, 2)}\n")
    
    # Check whether countAOI is greater than 0 to avoid division by zero
    if countAOI > 0:
        percentDetectableChange = round(((countPosArea+countNegArea)/countAOI)*100, 2)
    else:
        percentDetectableChange = 0 
    
    analyse.write(f"Thresholded Area of Interest [%] with Detectable Change: {percentDetectableChange}\n")
    analyse.write(f"Thresholded Volume [m3] of Surface Lowering: {round(countNeg, 2)}\n")
    analyse.write(f"Thresholded Volume [m3] of Surface Raising: {round(countPos, 2)}\n")
    analyse.write(f"Thresholded Volume [m3] of Difference: {round(countPos+countNeg, 2)}\n")
    analyse.write(f"Thresholded Net Volume [m3] of Difference: {round(countPos-countNeg, 2)}\n")
    analyse.write(f"Volume [m3] of Error within Threshold of {threshold}: {round(countAcc, 2)}\n")
    analyse.write(f"Count of Cells within {threshold}: {countCell}\n")
    analyse.write(f"Volume [m3] of Cells between {threshold}: {round(countCellVol, 2)}\n")
    analyse.write("\n\n")

    # Customize the CSV
    if threshold == "UC":
        analyse_csv.write(f"{threshold}; {countAOI}; {countPosArea+countNegArea}; {percentDetectableChange}; ; ; {countNeg}; {countPos}; {countPos+countNeg}; {countPos-countNeg}\n")
    else:
        analyse_csv.write(f"{threshold}; {countAOI}; {countPosArea+countNegArea}; {percentDetectableChange}; {countNegArea}; {countPosArea}; {countNeg}; {countPos}; {countPos+countNeg}; {countPos-countNeg}\n")

## ----------------------------------------------------------------------------------------------------------##
#Plot Data
plt.style.use('seaborn-darkgrid')  # Ein dunklerer Stil mit Gitterlinien

def plot_data(results, thresholds, final_path):
    # Preparing the data for plotting
    plotdata = pd.DataFrame({
        'Surface Raising': [results[str(th)]['volPos'] for th in thresholds],
        'Surface Lowering': [-results[str(th)]['volNeg'] for th in thresholds],
    }, index=[str(th) for th in thresholds])  

    # Plot settings
    New_Colors = ['darkblue', 'darkred']
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Verdana']

    # Creating the bar chart
    plotdata.plot(kind="bar", stacked=True, color=New_Colors, figsize=(10, 6))
    plt.xticks(rotation=45, horizontalalignment="right")
    plt.title('Volume Changes by Threshold [m³]')
    plt.xlabel('Threshold')
    plt.ylabel('Volume Changes [m³]')
    plt.grid(axis='y', linestyle='--')

    # Save plot
    plt.savefig(f'{final_path}/threshold_volume_changes.png', bbox_inches='tight')


## ----------------------------------------------------------------------------------------------------------##
# Save grid as raster
def save_grid_as_raster(filename, grid, geo_transform, epsg_code):
    driver = gdal.GetDriverByName("GTiff")
    rows, cols = grid.shape
    dataset = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    dataset.SetGeoTransform(geo_transform)
    #epsg_code = 32633 #anpassen
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg_code))
    projection = srs.ExportToWkt()
    dataset.SetProjection(projection)
    band = dataset.GetRasterBand(1)
    band.WriteArray(grid)
    band.SetNoDataValue(np.nan)
    band.FlushCache()
    dataset = None  # Close the file


## ----------------------------------------------------------------------------------------------------------##
# Visualize in QGIS
def visualize_qgis(thresholds, final_path):
    # Zugriff auf die Root-Ebene des Layer-Baums
    root = QgsProject.instance().layerTreeRoot()
    
    # Access to the root level of the layer tree
    change_detection_group = root.findGroup("ChangeDetection")
    if not change_detection_group:
        # Create group if it does not exist
        change_detection_group = root.addGroup("ChangeDetection")

    for threshold in thresholds:
        if threshold == "UC":
            threshold_label = 'uc'
        elif threshold == 0.0: 
            threshold_label = 'raw'
        else:
            threshold_value = float(threshold)
            threshold_label = f'{threshold_value}m'

        raster_path = os.path.join(final_path, f'corrected_diff_{threshold_label}.tif')

        raster_layer = QgsRasterLayer(raster_path, str(threshold_label))

        if not raster_layer.isValid():
            print("Layer konnte nicht geladen werden!")
        else:
            QgsProject.instance().addMapLayer(raster_layer, False)
            change_detection_group.addLayer(raster_layer)
            
            # Min/max values for automatic classification
            stats = raster_layer.dataProvider().bandStatistics(1, QgsRasterBandStats.All)  
            min_value = stats.minimumValue
            max_value = stats.maximumValue

            fcn = QgsColorRampShader()
            fcn.setColorRampType(QgsColorRampShader.Interpolated)

            colDic = {min_value: '#7a0615', -3: '#d7191c', -2: '#eb6640', -1: '#feb165', -0.3: '#ffdc96',
                      0.3: '#9cd3a7', 1: '#5ea7b1', 2: '#2b83ba', 3: '#1e5c83', max_value: '#05005b'}
            colorRampItemList = [QgsColorRampShader.ColorRampItem(value, QColor(color)) for value, color in colDic.items()]
            fcn.setColorRampItemList(colorRampItemList)

            shader = QgsRasterShader()
            shader.setRasterShaderFunction(fcn)

            renderer = QgsSingleBandPseudoColorRenderer(raster_layer.dataProvider(), 1, shader)
            raster_layer.setRenderer(renderer)

            # Transparency settings
            transparency = QgsRasterTransparency()
            transparentPixel = QgsRasterTransparency.TransparentSingleValuePixel()
            transparentPixel.min = -0.1
            transparentPixel.max = 0.1
            transparentPixel.percentTransparent = 100
            transparency.setTransparentSingleValuePixelList([transparentPixel])
            raster_layer.renderer().setRasterTransparency(transparency)

            raster_layer.triggerRepaint()

## ----------------------------------------------------------------------------------------------------------##
# Open Raster
def open_rasters(diff_path, uc_path):
    """Opens the difference and uncertainty rasters and reads their data."""
    try:
        diff_rds = gdal.Open(diff_path)
        uc_rds = gdal.Open(uc_path)
        if not diff_rds or not uc_rds:
            raise IOError("Could not open raster files.")
        pix_diff = diff_rds.GetRasterBand(1).ReadAsArray()
        pix_uc = uc_rds.GetRasterBand(1).ReadAsArray()
        return diff_rds, uc_rds, pix_diff, pix_uc
    except Exception as e:
        print(f"Error opening rasters: {e}")
        raise



## ----------------------------------------------------------------------------------------------------------##
## ----------------------------------------------------------------------------------------------------------##
## ----------------------------------------------------------------------------------------------------------##
class ChangeDetectionProcessingAlgorithm(QgsProcessingAlgorithm):
    """
Python Script to detect and analyse surface and volumetric changes in multi-temporal digital terrain models
Input- & Output-Parameter
    """
    INPUT_Shape = 'INPUT_SHP'
    INPUT_ALS_new = 'DTM1'
    INPUT_ALS_old = 'DTM2'
    SELECTION = 'Selection'
    INPUT_UC = 'UC'
    INPUT_threshold = 'Threshold'
    TEMP = 'TEMP'
    
    SELECTIONS = ['UncertaintyModel', 'Threshold']
    
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return ChangeDetectionProcessingAlgorithm()

    def name(self):
        return 'changedetectiontool'

    def displayName(self):
        return self.tr('QCD_ChangeDetection_v.2.0.0-beta')

    def group(self):
        return self.tr('pyAlpineRisk')

    def groupId(self):
        return 'pyAlpineRiskScripts'

    def shortHelpString(self):
        return self.tr("Detection and analysis of surface and volumetric changes in multi-temporal digital terrain models (Nicole Kamp, 2024) \nPUBLICATIONS\nKamp N., Krenn P., Avian M., Sass O. (2023). Comparability of multi‐temporal DTMs derived from different LiDAR platforms: Error sources and uncertainties in the application of geomorphic impact studies. In: Earth Surface Processes and Landforms 48(6): 1152–1175. \nKrenn P., Kamp N., Peßenteiner S., Sass O. (in review). Analysing geomorphic impacts of an extreme precipitation event on a torrential catchment in Upper Styria (Austria) using Unmanned-Aerial-Vehicle borne laser scanning (ULS). In: Earth Surface Processes and Landforms.")
        
        
    def initAlgorithm(self, config=None):
        self.selections = [self.tr('UncertaintyModel'),
                        self.tr('Threshold')]
                        
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_Shape,
                self.tr('Study Area'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_ALS_new,
                self.tr('New DTM Survey')
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_ALS_old,
                self.tr('Old DTM survey')
            )
        )
 
        self.addParameter(
            QgsProcessingParameterEnum(
                self.SELECTION,
                self.tr('Use calculated Uncertainty Model [m] or the threshold [m] of the minimum surface change'),
                self.selections,
                defaultValue=0
            )
        )
        
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_UC,
                self.tr('Uncertainty Model'),
                #extension='tiff',
                #fileFilter="tiff (*.tif)",
                optional=True
            )
        )
        
        self.addParameter(QgsProcessingParameterNumber(
            self.INPUT_threshold, 
            self.tr('Threshold'),
            QgsProcessingParameterNumber.Double,
            0.3,
            optional=True
            )
        )

        self.addParameter(
            QgsProcessingParameterFolderDestination(
                self.TEMP, 
                self.tr('Output Folder'), 
                defaultValue='C:/temp/CD'
            )
        )


    ## ----------------------------------------------------------------------------------------------------------##
    ## ----------------------------------------------------------------------------------------------------------##
    def processAlgorithm(self, parameters, context, feedback):
        feedback = QgsProcessingMultiStepFeedback(1, feedback)
        results = {}
        final_results = {}
        outputs = {}
        
        ## ----------------------------------------------------------------------------------------------------------##
        # Paths + Timestamp
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        temp_path = str(parameters[self.TEMP])+'/temp_'+str(timestamp)
        final_path = str(parameters[self.TEMP])+'/final_'+str(timestamp)
        
        
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
        if not os.path.exists(final_path):
            os.makedirs(final_path)
                
        ## Cell Size and EPSG-Code
        raster_cs = gdal.Open(str(parameters[self.INPUT_ALS_new]))
        gt_cs =raster_cs.GetGeoTransform()
        proj = osr.SpatialReference(wkt=raster_cs.GetProjection()) 
        cs = gt_cs[1]
        epsg = proj.GetAttrValue('AUTHORITY',1)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg))
        
        #Uncertainty Model - yes or no
        sel = self.SELECTIONS[self.parameterAsEnum(parameters, self.SELECTION, context)]

        ## ----------------------------------------------------------------------------------------------------------##
        ## Buffer
        buffer = out = temp_path + '/' + 'extent.shp'
        alg_params = {
            'DISSOLVE': False,
            'DISTANCE': -0.1,
            'END_CAP_STYLE': 2,
            'INPUT': str(parameters[self.INPUT_Shape]),
            'JOIN_STYLE': 1,
            'MITER_LIMIT': 2,
            'SEGMENTS': 1,
            'OUTPUT': buffer
        }
        processing.run('native:buffer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)

        ## ----------------------------------------------------------------------------------------------------------##
        ## Clip Rasters
        out_new = temp_path + '/' + 'dhm_new.tif'
        out_old = temp_path + '/' + 'dhm_old.tif'
        alg_params = {
            'ALPHA_BAND': False,
            'CROP_TO_CUTLINE': True,
            'DATA_TYPE': 0,
            'EXTRA': '',
            'INPUT': str(parameters[self.INPUT_ALS_new]),
            'KEEP_RESOLUTION': False,
            'MASK': str(buffer),
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'SET_RESOLUTION': False,
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg)),
            'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg)),
            'X_RESOLUTION': None,
            'Y_RESOLUTION': None,
            'OUTPUT': out_new
        }
        processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
        
        alg_params = {
            'ALPHA_BAND': False,
            'CROP_TO_CUTLINE': True,
            'DATA_TYPE': 0,
            'EXTRA': '',
            'INPUT': str(parameters[self.INPUT_ALS_old]),
            'KEEP_RESOLUTION': False,
            'MASK': str(buffer),
            'MULTITHREADING': False,
            'NODATA': None,
            'OPTIONS': '',
            'SET_RESOLUTION': False,
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg)),
            'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg)),
            'X_RESOLUTION': None,
            'Y_RESOLUTION': None,
            'OUTPUT': out_old
        }
        processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)


        ## ----------------------------------------------------------------------------------------------------------##
        ## ----------------------------------------------------------------------------------------------------------##
        ## Calculate DEM of Difference (DoD)
        format = "GTiff"
        driver = gdal.GetDriverByName("GTiff")
        if not driver:
            raise ValueError("GDAL TIFF driver not available.")
        
        # Open the input datasets and read data
        out2_old = temp_path + '/' + 'diff.tif'
        out2 = final_path + '/' + 'diff.tif'
        
        # DTM 2 - new DTM
        new_rds = gdal.Open(out_new)
        band_new = new_rds.GetRasterBand(1)
        pix_new = band_new.ReadAsArray()

        # DTM 1 - old DTM 
        old_rds = gdal.Open(out_old)
        band_old = old_rds.GetRasterBand(1)
        pix_old = band_old.ReadAsArray()
        
        # Get geo-transform and other properties from the old dataset
        gt = old_rds.GetGeoTransform()
        width = new_rds.RasterXSize
        height = new_rds.RasterYSize
#        minx = gt[0]
#        miny = gt[3] + width*gt[4] + height*gt[5] 
#        maxx = gt[0] + width*gt[1] + height*gt[2]
#        maxy = gt[3]
#        extent=str(minx)+","+str(maxx)+","+str(miny)+","+str(maxy)
#        extent2=(minx,maxx,miny,maxy)

        # Create new array and calculate difference
        grid_diff = np.full((height, width), np.nan, dtype=np.float32)
        nodata_value = -9999.0  # Define your actual NoData value
        valid_mask = (pix_new > nodata_value) & (pix_old > nodata_value)
        grid_diff[valid_mask] = pix_new[valid_mask] - pix_old[valid_mask]

        # Create a new raster for the difference
        out_diff_raster = driver.Create(out2_old, width, height, 1, gdal.GDT_Float32)
        out_diff_raster.SetGeoTransform(gt)  # Set the same geo-transform as the input
        out_diff_raster.SetProjection(new_rds.GetProjection())  # Set the same projection as the new_rds

        # Write the difference array to the new raster
        out_diff_band = out_diff_raster.GetRasterBand(1)
        out_diff_band.WriteArray(grid_diff)
        out_diff_band.SetNoDataValue(np.nan)  # Set NoData value to NaN or another appropriate value
        out_diff_band.FlushCache()

        # Copy the temporary raster to the final location
        driver.CreateCopy(out2, out_diff_raster, 0)

        # Cleanup
        out_diff_band = None
        out_diff_raster = None
        dst_ds = None
        new_rds = None
        old_rds = None
        
        
        ## ----------------------------------------------------------------------------------------------------------##
        ## ----------------------------------------------------------------------------------------------------------##
        ## ChangeDetection (Uncertainty Model)
        if str(sel) == str('UncertaintyModel'):
            feedback.pushInfo(str(sel))
            out_uc = temp_path + '/' + 'dhm_uc.tif'
            alg_params = {
                'ALPHA_BAND': False,
                'CROP_TO_CUTLINE': True,
                'DATA_TYPE': 0,
                'EXTRA': '',
                'INPUT': str(parameters[self.INPUT_UC]),
                'KEEP_RESOLUTION': False,
                'MASK': str(buffer),
                'MULTITHREADING': False,
                'NODATA': None,
                'OPTIONS': '',
                'SET_RESOLUTION': False,
                'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg)),
                'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:'+str(epsg)),
                'X_RESOLUTION': None,
                'Y_RESOLUTION': None,
                'OUTPUT': out_uc
            }
            processing.run('gdal:cliprasterbymasklayer', alg_params, context=context, feedback=feedback, is_child_algorithm=True)
            
            # Opening output files
            analyse = open(final_path + '/' + 'change_detection_UC.txt', "w")
            analyse_csv = open(final_path + '/' + 'change_detection_UC.csv', "w")
            analyse_csv.write(" ; AOI [m2]; Detectable Change [m2]; Detectable Change [%]; Surface Lowering [m2]; Surface Raising [m2]; Surface Lowering [m3]; Surface Raising [m3]; Volume of Difference [m3]; Net Volume of Difference [m3]" + "\n")

            # Rasters to array
            diff_array = raster_to_array(out2)
            uc_array = raster_to_array(out_uc)
            
            #Count AOI for raw DoD
            countAOI = 0
            for row in range(height):
                for col in range(width):
                    diff = diff_array[row, col]
                    if diff > -9999.0 and diff < 100:
                        countAOI += cs * cs 
            
            #Analyse
            thresholds = [0, 0.1, 0.3, "UC"]
            for threshold in thresholds:
                results[str(threshold)] = process_threshold(threshold, diff_array, uc_array, cs, width, height, countAOI, analyse, analyse_csv, final_path, gt_cs, epsg)

            # Close all files
            analyse.close()
            analyse_csv.close()
            
            #Plot Results
            plot_data(results, thresholds, final_path)


        ## ----------------------------------------------------------------------------------------------------------##
        ## ChangeDetection (Threshold)
        if str(sel) == str('Threshold'):
            feedback.pushInfo(str(sel))

            # Opening output files
            analyse = open(final_path + '/' + 'change_detection_UC.txt', "w")
            analyse_csv = open(final_path + '/' + 'change_detection_UC.csv', "w")
            analyse_csv.write(" ; AOI [m2]; Detectable Change [m2]; Detectable Change [%]; Surface Lowering [m2]; Surface Raising [m2]; Surface Lowering [m3]; Surface Raising [m3]; Volume of Difference [m3]; Net Volume of Difference [m3]" + "\n")

            # Rasters to array
            diff_array = raster_to_array(out2)
            
            #Count AOI for raw DoD
            countAOI = 0
            for row in range(height):
                for col in range(width):
                    diff = diff_array[row, col]
                    if diff > -9999.0 and diff < 100:
                        countAOI += cs * cs 
            
            #Analyse
            val = (parameters[self.INPUT_threshold])
            uc_array = np.array([])
            thresholds = [0, 0.1, 0.3, val]
            for threshold in thresholds:
                results[str(threshold)] = process_threshold(threshold, diff_array, uc_array, cs, width, height, countAOI, analyse, analyse_csv, final_path, gt_cs, epsg)

            # Close all files
            analyse.close()
            analyse_csv.close()
        
            #Plot Results
            plot_data(results, thresholds, final_path)

        ## ----------------------------------------------------------------------------------------------------------##
        ## ----------------------------------------------------------------------------------------------------------##
        # Visualize in QGIS
        visualize_qgis(thresholds, final_path)


        ## ----------------------------------------------------------------------------------------------------------##
        ## ----------------------------------------------------------------------------------------------------------##
        ## Process completed 
        #the_end = 'The End'
        #outputs['LastStep'] = the_end
        final_results['ChangeDetection completed: '] = results
        return final_results


