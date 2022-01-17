# -*- coding: utf-8 -*-
''' Tools to create an structural analysis model from an IFC structural model.'''

__author__= "Luis C. Pérez Tato (LCPT) and Ana Ortega (AOO)"
__copyright__= "Copyright 2021, LCPT and AOO"
__license__= "GPL"
__version__= "3.0"
__email__= "l.pereztato@ciccp.es"

import os
import pwd
from datetime import datetime

import ifcopenshell

class StructureDistiller(object):
    ''' Tools to create an structural analysis model from an IFC 
        structural model. Inspired on ca2ifc module from IfcOpenShell-0.6.0.

    :ivar ifcModel: IFC structural analysis model.
    :ivar outputFileName: name of the output file.
    '''
    def __init__(self, outputFileName: str, modelName: str):
        ''' Constructor. 

        :param outputFileName: name of the output file.
        :param modelName: human readable model name.
        '''
        self.ifcModel= None
        self.outputFileName= outputFileName
        self.modelName= modelName
        self.reps= dict()
        self.ownerHistory= None
        self.globalAxes= None
        self.localPlacement= None
        self.building= None
        self.analysisModel= None
        self.ifcElements= list()

    def createHeader(self):
        self.ifcModel.wrapped_data.header.file_name.name = os.path.basename(self.outputFileName)

    def getUserName(self):
        ''' Return the name of the user that is running this code.'''
        return pwd.getpwuid( os.getuid() )[0]
        
    def createGlobalAxes(self):
        ''' Create global axes for the structural analysis model.'''
        xAxis= self.ifcModel.createIfcDirection((1.0, 0.0, 0.0))
        yAxis= self.ifcModel.createIfcDirection((0.0, 1.0, 0.0))
        zAxis= self.ifcModel.createIfcDirection((0.0, 0.0, 1.0))
        origin= self.ifcModel.createIfcCartesianPoint((0.0, 0.0, 0.0))
        self.globalAxes= self.ifcModel.createIfcAxis2Placement3D(origin, zAxis, xAxis)

    def createLocalPlacement(self, relative_to= None):
        if(self.globalAxes==None):
            self.createGlobalAxes()
        return self.ifcModel.createIfcLocalPlacement(relative_to,self.globalAxes)

    def createOwnerHistory(self):
        ''' Create IFC owner history. IfcOwnerHistory is used to identify 
        the creating and owning application and user for the associated 
        object, as well as capture the last modifying application and user.

        '''
        actor = self.ifcModel.createIfcActorRole("ENGINEER", None, None)
        userName= self.getUserName()
        person = self.ifcModel.createIfcPerson(userName, None, '', None, None, None, (actor,))
        organization = self.ifcModel.createIfcOrganization(
            None,
            "IfcOpenShell",
            "IfcOpenShell, an open source (LGPL) software library that helps users and software developers to work with the IFC file format.",
        )
        p_o = self.ifcModel.createIfcPersonAndOrganization(person, organization)
        application = self.ifcModel.createIfcApplication(organization, "v0.0.x", "IFC2CA", "IFC2CA")
        timestamp = int(datetime.now().timestamp())
        return self.ifcModel.createIfcOwnerHistory(p_o, application, "READWRITE", None, None, None, None, timestamp)
    
    def createReferenceSubrep(self):
        modelRep = self.ifcModel.createIfcGeometricRepresentationContext(None, "Model", 3, 1.0e-05, self.globalAxes, None)
        bodySubRep = self.ifcModel.createIfcGeometricRepresentationSubContext("Body", "Model", None, None, None, None, modelRep, None, "MODEL_VIEW", None)
        refSubRep = self.ifcModel.createIfcGeometricRepresentationSubContext("Reference", "Model", None, None, None, None, modelRep, None, "GRAPH_VIEW", None)

        return {"model": modelRep, "body": bodySubRep, "reference": refSubRep}
    
    def guid(self):
        return ifcopenshell.guid.new()

    def setup(self):
        ''' IFC structural analysis model setup.'''
        self.ifcModel= ifcopenshell.file() # initiate ifc file
        self.createHeader() # create header
        
        # create global axes
        self.localPlacement= self.createLocalPlacement()
        
        # create units
        lengthUnit= self.ifcModel.createIfcSIUnit(None, "LENGTHUNIT", None, "METRE")
        unitAssignment= self.ifcModel.createIfcUnitAssignment((lengthUnit,))

        self.ownerHistory= self.createOwnerHistory() # create owner history

        self.reps = self.createReferenceSubrep() # create representations and subrepresentations
        
        # IFC hierarcy creation.
        
        # create project
        project = self.ifcModel.createIfcProject(
            self.guid(), self.ownerHistory, "A Project", None, None, None, None, (self.reps["model"],), unitAssignment
        )

        sitePlacement= self.createLocalPlacement()
        site= self.ifcModel.createIfcSite(self.guid(), self.ownerHistory, "Site", None, None, sitePlacement, None, None, "ELEMENT", None, None, None, None, None)
        buildingPlacement= self.createLocalPlacement()
        self.building= self.ifcModel.createIfcBuilding(self.guid(), self.ownerHistory, 'Building', None, None, buildingPlacement, None, None, "ELEMENT", None, None, None)

        containerSite = self.ifcModel.createIfcRelAggregates(self.guid(), self.ownerHistory, "Site Container", None, site, [self.building])
        containerProject = self.ifcModel.createIfcRelAggregates(self.guid(), self.ownerHistory, "Project Container", None, project, [site])
        

        
        # create model
        self.analysisModel= self.ifcModel.createIfcStructuralAnalysisModel(
            self.guid(),
            self.ownerHistory,
            self.modelName,
            None,
            None,
            "NOTDEFINED",
            self.globalAxes,
            None,
            None,
            self.localPlacement,
        )
        self.ifcModel.createIfcRelDeclares(self.guid(), self.ownerHistory, None, None, project, (self.analysisModel,))
        

    def createSurfaceGeometry(self, wires):
        ''' Create an IfcFaceSurface object from the wires argument.

        :param wires: polygons defining the surface and its holes (if any).
        '''
        sz= len(wires)
        faceBounds= list()
        if(sz>0): # Outer bound.
            outerBoundWire= wires[0]
            ifcPoints= list()
            for p in outerBoundWire[:-1]:
                ifcPoints.append(self.ifcModel.createIfcCartesianPoint((p.x, p.y, p.z)))
            polyLoop= self.ifcModel.createIfcPolyLoop(ifcPoints)
            outerBound= self.ifcModel.createIfcFaceOuterBound(Bound=polyLoop, Orientation=True)  # orientation of vertices is CCW
            faceBounds.append(outerBound)
        if(sz>1): # Inner bounds.
            for interiorWire in wires[1:]:
                ifcPoints= list()
                for p in interiorWire[:-1]:
                    ifcPoints.append(self.ifcModel.createIfcCartesianPoint((p.x, p.y, p.z)))
                polyLoop= self.ifcModel.createIfcPolyLoop(ifcPoints)
                innerBound= self.ifcModel.createIfcFaceBound(Bound=polyLoop, Orientation=False) # orientation of vertices is CW
                faceBounds.append(innerBound)
        return self.ifcModel.createIfcFaceSurface(Bounds= faceBounds)

    def createDefinitionShape(self, wires):
        ''' Create the product definition shape.

        :surfaceData: 
        '''
        ifcFaceSurface= self.createSurfaceGeometry(wires)
        topologyRep = self.ifcModel.createIfcTopologyRepresentation(self.reps["reference"], "Reference", "Face", (ifcFaceSurface,))
        return self.ifcModel.createIfcProductDefinitionShape(None, None, (topologyRep,))

    def createMaterial(self, materialData):
        ''' Create IFC material from the data argument.

        :param materialData: dictionary containing material names and 
                             thicknesses for each layer. 
        '''
        materialLayers= list()
        for matData in materialData:
            name= matData['materialName']
            if(name==None):
                name= 'air'
            thickness= matData['thickness']
            print('name= ', name)
            material= self.ifcModel.createIfcMaterial(name)
            materialLayer= self.ifcModel.createIfcMaterialLayer(material, thickness, None)
            materialLayers.append(materialLayer)
            print(materialLayers)
            
        materialLayerSet= self.ifcModel.createIfcMaterialLayerSet(materialLayers, None)
        materialLayerSetUsage= self.ifcModel.createIfcMaterialLayerSetUsage(materialLayerSet, "AXIS2", "POSITIVE", 0.0)
        return materialLayerSetUsage
        


    def dumpSurfaces(self, surfaces):
        ''' Create the IfcStructuralSurfaceMember objects corresponding to
            the surface arguments.

        :param surfaces: dictionary containing surface data.
        '''
        for shapeKey in surfaces:
            shapeData= surfaces[shapeKey]
            wires= shapeData['midSectionWires']
            prodDefShape = self.createDefinitionShape(wires)
            label= str(shapeKey)
            thickness= shapeData['thickness']
            materialLayerSetUsage= self.createMaterial(shapeData['materials'])
            surface= self.ifcModel.createIfcStructuralSurfaceMember(self.guid(), self.ownerHistory, label, None, None, self.localPlacement, prodDefShape, "SHELL", thickness)
            self.ifcElements.append(surface)
            self.ifcModel.createIfcRelAssociatesMaterial(self.guid(), self.ownerHistory, RelatedObjects=[surface], RelatingMaterial= materialLayerSetUsage)
        
    def write(self):
        ''' Writes the IFC file.'''
        # assign elements
        self.ifcModel.createIfcRelAssignsToGroup(self.guid(), self.ownerHistory, None, None, tuple(self.ifcElements), None, self.analysisModel)
        self.ifcModel.write(self.outputFileName)