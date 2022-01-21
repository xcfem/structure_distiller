# -*- coding: utf-8 -*-
''' Compute wall mid-plane for use in structural analytical model.'''

__author__= "Luis C. Pérez Tato (LCPT) and Ana Ortega (AOO)"
__copyright__= "Copyright 2021, LCPT and AOO"
__license__= "GPL"
__version__= "3.0"
__email__= "l.pereztato@ciccp.es"

import os
from structure_distiller import test_utils
                    
inputFileName= '/../data/south_roof.ifc'
outputFileName= '/tmp/test.ifc'
projectName= os.path.basename(inputFileName).replace('.ifc','')

pth= os.path.dirname(__file__)
ifcModel= test_utils.openTestFile(pth, inputFileName)
test_utils.extractDatum(ifcModel, outputFileName, projectName)
        
        
