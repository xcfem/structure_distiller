# -*- coding: utf-8 -*-
''' Utilities for library testing.'''

import ifcopenshell
from structure_distiller import structure_distiller
from structure_distiller import ifc_utils

def openTestFile(pth, ifcFileName):
    ''' Open a test file from the data directory.

    :param ifcFileName: name of the IFC file to open.
    '''
    if(not pth):
      pth= "."
    fName= pth+"/../data/"+ifcFileName
    return ifcopenshell.open(fName)

def extractDatum(ifcModel, ifcOutputFileName, projectName, modelName= None):
    ''' Extract mid-surfraces from the IFC products of the input model.

    :param ifcModel: IFC input model.
    :param ifcOutputFile: name of the file for the resulting IFC model.
    :param projectName: human readable project name.
    :param modelName: human readable model name.
    '''
    # The geometric elements in an IFC file are the IfcProduct elements. So
    # these are stored in product_shapes.
    productShapes= ifc_utils.getProductShapes(ifcModel)
    midSurfaces, axes= ifc_utils.computeShapesDatum(productShapes)

    # Write ouput
    outputModel= structure_distiller.StructureDistiller(outputFileName= ifcOutputFileName, projectName= projectName, modelName= modelName)
    outputModel.setup()
    outputModel.dumpSurfaces(midSurfaces)
    outputModel.dumpLines(axes)
    outputModel.write()
