# -*- coding: utf-8 -*-
''' Compute wall mid-plane for use in structural analytical model.'''

__author__= "Luis C. Pérez Tato (LCPT) and Ana Ortega (AOO)"
__copyright__= "Copyright 2021, LCPT and AOO"
__license__= "GPL"
__version__= "3.0"
__email__= "l.pereztato@ciccp.es"

import ifcopenshell
from structure_distiller import structure_distiller
from structure_distiller import ifc_utils
            
        
        
# Open the IFC file using IfcOpenShell
import os
pth= os.path.dirname(__file__)
# print("pth= ", pth)
if(not pth):
  pth= "."
fName= pth+"/../data/revit2011_wall1.ifc"
inputFile= ifcopenshell.open(fName)


# The geometric elements in an IFC file are the IfcProduct elements. So
# these are stored in product_shapes.

productShapes= ifc_utils.getProductShapes(inputFile)

midSurfaces= ifc_utils.computeMidSurfaces(productShapes) 

# Write ouput
outputModel= structure_distiller.StructureDistiller(outputFileName= '/tmp/test.ifc', modelName= 'My Model')

outputModel.setup()
outputModel.dumpSurfaces(midSurfaces)
outputModel.write()


