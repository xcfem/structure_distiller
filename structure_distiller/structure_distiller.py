# -*- coding: utf-8 -*-
''' Tools to create an structural analysis model from an IFC structural model.'''

__author__= "Luis C. Pérez Tato (LCPT) and Ana Ortega (AOO)"
__copyright__= "Copyright 2021, LCPT and AOO"
__license__= "GPL"
__version__= "3.0"
__email__= "l.pereztato@ciccp.es"

import xc_base
import geom

import os
import sys
import pwd
from datetime import datetime


import ifcopenshell
from misc_utils import log_messages as lmsg

class StructureDistiller(object):
    ''' Tools to create an structural analysis model from an IFC 
        structural model. Inspired on ca2ifc module from IfcOpenShell-0.6.0.

    :ivar ifcModel: IFC structural analysis model.
    :ivar outputFileName: name of the output file.
    '''
    applicationName= 'XCStructureDistiller'
    def __init__(self, outputFileName: str, projectName: str, modelName: str= None):
        ''' Constructor. 

        :param outputFileName: name of the output file.
        :param projectName: human readable project name.
        :param modelName: human readable model name.
        '''
        self.ifcModel= None
        self.outputFileName= outputFileName
        self.projectName= projectName
        if(modelName==None):
            self.modelName= projectName
        else:
            self.modelName= modelName
        self.reps= dict()
        self.ownerHistory= None
        self.globalAxes= None
        self.localPlacement= None
        self.building= None
        self.analysisModel= None
        self.ifcElements= list()
        self.ifcMaterials= list()
        self.ifcProfiles= dict()

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
        application = self.ifcModel.createIfcApplication(organization, "v0.0.x", self.applicationName, self.applicationName)
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
            self.guid(), self.ownerHistory, self.projectName, None, None, None, None, (self.reps["model"],), unitAssignment
        )

        # sitePlacement= self.createLocalPlacement()
        # site= self.ifcModel.createIfcSite(self.guid(), self.ownerHistory, "Site", None, None, sitePlacement, None, None, "ELEMENT", None, None, None, None, None)
        # buildingPlacement= self.createLocalPlacement()
        # self.building= self.ifcModel.createIfcBuilding(self.guid(), self.ownerHistory, 'Building', None, None, buildingPlacement, None, None, "ELEMENT", None, None, None)

        # containerSite = self.ifcModel.createIfcRelAggregates(self.guid(), self.ownerHistory, "Site Container", None, site, [self.building])
        # containerProject = self.ifcModel.createIfcRelAggregates(self.guid(), self.ownerHistory, "Project Container", None, project, [site])
        

        
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
        

    # def createSurfaceGeometry(self, wires):
    #     ''' Create an IfcFaceSurface object from the wires argument.

    #     :param wires: polygons defining the surface and its holes (if any).
    #     '''
    #     sz= len(wires)
    #     faceBounds= list()
    #     if(sz>0): # Outer bound.
    #         outerBoundWire= wires[0]
    #         ifcPoints= list()
    #         for p in outerBoundWire[:-1]:
    #             ifcPoints.append(self.ifcModel.createIfcCartesianPoint((p.x, p.y, p.z)))
    #         polyLoop= self.ifcModel.createIfcPolyLoop(ifcPoints)
    #         outerBound= self.ifcModel.createIfcFaceOuterBound(Bound=polyLoop, Orientation=True)  # orientation of vertices is CCW
    #         faceBounds.append(outerBound)
    #     if(sz>1): # Inner bounds.
    #         for interiorWire in wires[1:]:
    #             ifcPoints= list()
    #             for p in interiorWire[:-1]:
    #                 ifcPoints.append(self.ifcModel.createIfcCartesianPoint((p.x, p.y, p.z)))
    #             polyLoop= self.ifcModel.createIfcPolyLoop(ifcPoints)
    #             innerBound= self.ifcModel.createIfcFaceBound(Bound=polyLoop, Orientation=False) # orientation of vertices is CW
    #             faceBounds.append(innerBound)
    #     return self.ifcModel.createIfcFaceSurface(Bounds= faceBounds)

    def createFaceBound(self, plg, outerBound):
        ''' Create an IfcFaceOuterBound or an IfcBound object from 
            the polygon argument.

        :param plg: point sequence.
        :param outerBound: if true the polygon defines an outer bound 
                           of the face.
        '''
        retval= None
        # Create vertices.
        ifcVertices= list()
        for p in plg[:-1]:
            pt= self.ifcModel.createIfcCartesianPoint((p.x, p.y, p.z))
            vertex= self.ifcModel.createIfcVertexPoint(pt)
            ifcVertices.append(vertex)
        # Create edges.
        orientedEdges= list()
        for i, v in enumerate(ifcVertices):
            v2Index= (i + 1) if i < len(ifcVertices) - 1 else 0
            edge= self.ifcModel.createIfcEdge(v, ifcVertices[v2Index])
            orientedEdges.append(self.ifcModel.createIfcOrientedEdge(None, None, edge, True))
        # Create edge loop.
        edgeLoop= self.ifcModel.createIfcEdgeLoop(tuple(orientedEdges))
        if(outerBound):
            retval= self.ifcModel.createIfcFaceOuterBound(Bound= edgeLoop, Orientation= True)  # orientation of vertices is CCW
        else:
            retval= self.ifcModel.createIfcFaceBound(Bound= edgeLoop, Orientation= False) # orientation of vertices is CW
        return retval

    def createPlane(self, plg):
        ''' Create an IfcPlane object from the polygon argument.

        :param plg: point sequence.
        '''
        points= list()
        for p in plg[:-1]:
            points.append(p)
        plane= geom.Plane3d(points)
        xAxis= plane.getBase1()
        ifcXAxis= self.ifcModel.createIfcDirection((xAxis.x,xAxis.y,xAxis.z))
        zAxis= plane.getNormal()
        ifcZAxis= self.ifcModel.createIfcDirection((zAxis.x,zAxis.y,zAxis.z))
        origin= plg.getCenterOfMass()
        ifcOrigin= self.ifcModel.createIfcCartesianPoint((origin.x, origin.y, origin.z))
        localAxes= self.ifcModel.createIfcAxis2Placement3D(ifcOrigin, ifcZAxis, ifcXAxis)
        return self.ifcModel.createIfcPlane(localAxes)
        
                
    def createSurfaceGeometry(self, wires):
        ''' Create an IfcFaceSurface object from the wires argument.

        :param wires: polygons defining the surface and its holes (if any).
        '''
        sz= len(wires)
        retval= None
        if(sz>0): # Outer bound.
            faceBounds= list()
            outerBoundWire= wires[0]
            faceBounds.append(self.createFaceBound(plg= outerBoundWire, outerBound= True))
            ifcPlane= self.createPlane(plg= outerBoundWire)
            
            if(sz>1): # Ther are inner bounds.
                for interiorWire in wires[1:]:
                    faceBounds.append(self.createFaceBound(plg= interiorWire, outerBound= False))
            retval= self.ifcModel.createIfcFaceSurface(faceBounds, ifcPlane, True)
        else:
            className= type(self).__name__
            methodName= sys._getframe(0).f_code.co_name
            lmsg.error(className+'.'+methodName+': no input polygons.')
        return retval
                
    def createSurfaceDefinitionShape(self, wires):
        ''' Create the surface product definition shape.

        :param wires: polygons defining the surface and its holes (if any).
        '''
        ifcFaceSurface= self.createSurfaceGeometry(wires)
        topologyRep = self.ifcModel.createIfcTopologyRepresentation(self.reps["reference"], "Reference", "Face", (ifcFaceSurface,))
        return self.ifcModel.createIfcProductDefinitionShape(None, None, (topologyRep,))

    def createLineGeometry(self, segment):
        ''' Create the IfcEdge object that represents the segment argument.

        :param segment: 3D segment defining the axis.
        '''
        p0= segment.getFromPoint()
        p1= segment.getToPoint()
        startPoint = self.ifcModel.createIfcCartesianPoint((p0.x, p0.y, p0.z))
        startVertex = self.ifcModel.createIfcVertexPoint(startPoint)
        endPoint = self.ifcModel.createIfcCartesianPoint((p1.x, p1.y, p1.z))
        endVertex = self.ifcModel.createIfcVertexPoint(endPoint)
        return self.ifcModel.createIfcEdge(startVertex, endVertex)
        
    def createLineDefinitionShape(self, segment):
        ''' Create the linear product definition shape.

        :param segment: 3D segment defining the axis.
        '''
        ifcEdge= self.createLineGeometry(segment)
        topologyRep = self.ifcModel.createIfcTopologyRepresentation(self.reps["reference"], "Reference", "Edge", (ifcEdge,))
        return self.ifcModel.createIfcProductDefinitionShape(None, None, (topologyRep,))

    def createLocalZAxis(self, refSys):
        ''' Create local z axis from the system of reference argument.

        :param refSys: reference system that define the measurement directions. 
        '''
        zVector= refSys.getKVector()
        return self.ifcModel.createIfcDirection((zVector.x, zVector.y, zVector.z))

    def createMaterial(self, materialName):
        ''' Create an IfcMaterial object.

        :param MaterialName: name of the material.
        '''
        retval= None
        if(materialName in self.ifcMaterials):
            retval= self.ifcMaterials[profileName]
        else:
            retval= self.ifcModel.createIfcMaterial(Name= materialName)
            self.ifcProfiles[materialName]= retval
        return retval        
 
    def createMaterialLayers(self, materialData):
        ''' Create IFC material from the data argument.

        :param materialData: dictionary containing material names and 
                             thicknesses for each layer. 
        '''
        materialLayers= list()
        retval= None
        if(len(materialData)>0):       
            for matData in materialData:
                name= matData['materialName']
                if(name==None):
                    name= 'air'
                thickness= matData['thickness']
                material= self.createMaterial(name)
                materialLayer= self.ifcModel.createIfcMaterialLayer(material, thickness, None)
                materialLayers.append(materialLayer)
            materialLayerSet= self.ifcModel.createIfcMaterialLayerSet(materialLayers, None)
            retval= self.ifcModel.createIfcMaterialLayerSetUsage(materialLayerSet, "AXIS2", "POSITIVE", 0.0)
        return retval

    def createProfile(self, profileName):
        ''' Creates an IfcMaterialProfileSet object.

        :param profileName: name of the profile.
        '''
        retval= None
        if(profileName in self.ifcProfiles):
            retval= self.ifcProfiles[profileName]
        else:
            retval= self.createMaterial(profileName)
            # matProf= self.ifcModel.createIfcMaterialProfile(profileName, None, material)
            # profileSet= self.ifcModel.createIfcMaterialProfileSet(None, None, (matProf,))
            # retval= self.ifcModel.createIfcMaterialProfileSetUsage(profileSet)
            self.ifcProfiles[profileName]= retval
        return retval

    def dumpSurfaces(self, surfaces):
        ''' Create the IfcStructuralSurfaceMember objects corresponding to
            the surface arguments.

        :param surfaces: dictionary containing surface data.
        '''
        for shapeKey in surfaces:
            shapeData= surfaces[shapeKey]
            wires= shapeData['midSectionWires']
            prodDefShape = self.createSurfaceDefinitionShape(wires)
            label= shapeData['product'].Name
            if(label==None):
                label= str(shapeKey)
            thickness= shapeData['thickness']
            surface= self.ifcModel.createIfcStructuralSurfaceMember(self.guid(), self.ownerHistory, label, None, None, self.localPlacement, prodDefShape, "SHELL", thickness)
            self.ifcElements.append(surface)
            # Set material.
            materialLayerSetUsage= self.createMaterialLayers(shapeData['materials'])
            if(materialLayerSetUsage): # There are materials to link to the surface.
                self.ifcModel.createIfcRelAssociatesMaterial(self.guid(), self.ownerHistory, RelatedObjects=[surface], RelatingMaterial= materialLayerSetUsage)
                
    def dumpLines(self, lines):
        ''' Create the IfcStructuralSurfaceMember objects corresponding to
            the surface arguments.

        :param lines: dictionary containing line data.
        '''
        for shapeKey in lines:
            shapeData= lines[shapeKey]
            axis= shapeData['axis']
            prodDefShape= self.createLineDefinitionShape(axis)
            label= shapeData['product'].Name
            if(label==None):
                label= str(shapeKey)
            length= shapeData['length']
            footprint= shapeData['footprint']
            localZAxis= self.createLocalZAxis(footprint.refSys)
            line= self.ifcModel.createIfcStructuralCurveMember(self.guid(), self.ownerHistory, label, None, None, self.localPlacement, prodDefShape, "RIGID_JOINED_MEMBER", localZAxis)
            self.ifcElements.append(line)
            # Set material.
            profileName= shapeData['product'].Description
            materials= shapeData['materials']
            if(len(materials)>0):
                if(isinstance(materials, list)):
                    profileName= materials[0]
                else: # dict
                    profileName= materials[0]['materialName']
            profile= self.createProfile(profileName)
            if(profile): # There is a profile to link to the surface.
                self.ifcModel.createIfcRelAssociatesMaterial(self.guid(), self.ownerHistory, RelatedObjects=[line], RelatingMaterial= profile)
        
    def write(self):
        ''' Writes the IFC file.'''
        # assign elements
        self.ifcModel.createIfcRelAssignsToGroup(self.guid(), self.ownerHistory, None, None, tuple(self.ifcElements), None, self.analysisModel)
        self.ifcModel.write(self.outputFileName)
