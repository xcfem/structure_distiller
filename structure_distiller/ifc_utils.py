# -*- coding: utf-8 -*-
''' Compute shell mid-surface and beam-column axes for use in structural 
    analytical model.'''

__author__= "Luis C. Pérez Tato (LCPT) and Ana Ortega (AOO)"
__copyright__= "Copyright 2021, LCPT and AOO"
__license__= "GPL"
__version__= "3.0"
__email__= "l.pereztato@ciccp.es"

import ifcopenshell
import ifcopenshell.geom

import OCC.Core.gp
import OCC.Core.BRepBuilderAPI
import OCC.Core.BRepAlgoAPI
import OCC.Core.TopExp
import OCC.Core.ShapeAnalysis
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop_VolumeProperties
from OCC.Core.BRepBndLib import brepbndlib_Add
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.Core.BRep import BRep_Tool

import xc_base
import geom
import numpy as np

def getBoundingBox(shape, tol:float =1e-6) -> tuple:
    """ return the bounding box of the TopoDS_Shape `shape`


    :param shape : TopoDS_Shape or a subclass such as TopoDS_Face
                   the shape to compute the bounding box from
    :param tol: tolerance of the computed boundingbox
    """
    bbox= Bnd_Box()
    bbox.SetGap(tol)
    brepbndlib_Add(shape, bbox, False)

    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    return ((xmin, ymin, zmin), (xmax, ymax, zmax))

def getLocalReferenceSystem(shape):
    ''' Return a local reference system for the shape argument using
        the principal axis of inertia.

    :param shape : TopoDS_Shape to compute the reference system for.
    '''
    # Compute inertia properties
    props= GProp_GProps()
    brepgprop_VolumeProperties(shape, props)

    # Get centroid
    cog = props.CentreOfMass()
    cog_x, cog_y, cog_z = cog.Coord()
    centroidPos= geom.Pos3d(cog_x, cog_y, cog_z)

    # Get inertia tensor.
    _m = props.MatrixOfInertia()
    matrix_of_inertia = [
         (_m.Row(row + 1).X(), _m.Row(row + 1).Y(), _m.Row(row + 1).Z())
         for row in range(3)
     ]
    eigenvalues, eigenvectors = np.linalg.eig(np.array(matrix_of_inertia))
    v0= geom.Vector3d(eigenvectors[0][0], eigenvectors[0][1], eigenvectors[0][2])
    v1= geom.Vector3d(eigenvectors[1][0], eigenvectors[1][1], eigenvectors[1][2])
    v2= geom.Vector3d(eigenvectors[2][0], eigenvectors[2][1], eigenvectors[2][2])

    return geom.Ref3d3d(centroidPos, centroidPos+v0, centroidPos+v1)

def computeShapeDimensions(shape, refSys):
    ''' Compute the dimensions of the shape arguments measures along the axes
        of the reference system argument.

    :param shape : TopoDS_Shape to compute the reference system for.
    :param refSys: reference system that define the measurement directions. 
    '''
    box= getBoundingBox(shape)
    pMin= geom.Pos3d(box[0][0], box[0][1], box[0][2])
    pMax= geom.Pos3d(box[1][0], box[1][1], box[1][2])

    localMin= refSys.getLocalPosition(pMin)
    localMax= refSys.getLocalPosition(pMax)
    return [abs(localMax.x-localMin.x), abs(localMax.y-localMin.y), abs(localMax.z-localMin.z)]

def computeDimensionality(shapeDimensions, thresholdFactor= 0.2):
    ''' Compute the dimensionality from the shape dimensions.

    :param shapeDimensions: dimensions of the shape along its local axes.
    :param thresholdFactor: factor to consider a dimension as of a lower order
                            (if dim/maxDim<thresholdFactor => lower order dimension). 
    '''
    retval= 0
    maxDimension= max(shapeDimensions) # Compute maximum dimension.
    normalizedDimensions= [shapeDimensions[0]/maxDimension, shapeDimensions[1]/maxDimension, shapeDimensions[2]/maxDimension]
    for d in normalizedDimensions:
        if(d>thresholdFactor):
            retval+= 1
    return retval

class GeometryFootprint(object):
    ''' geometric parameters that will allow to get 
        the mid-surface or the axis of the shape.

    :ivar refSys: local reference system for the shape computed using
                  its principal axis of inertia.
    :ivar shapeDimensions: dimensions of the shape along its local axes.
    :ivar dimensionality: dimensionality of the shape (0, 1, 2 or 3).
    '''

    def __init__(self, shape, thresholdFactor= 0.2):
        ''' Compute some geometric parameters that will allow to get 
            the mid-surface or the axis of the shape.

        :param shape : TopoDS_Shape to compute the reference system for.
        :param thresholdFactor: factor to consider a dimension as of a lower order
                                (if dim/maxDim<thresholdFactor => lower order dimension). 
        '''
        self.refSys= getLocalReferenceSystem(shape)
        self.shapeDimensions= computeShapeDimensions(shape, self.refSys)
        self.dimensionality= computeDimensionality(self.shapeDimensions)

    def getCentroid(self):
        ''' Return the origin of the local coordinate system.'''
        return self.refSys.getOrg()

    def getAxisDirections(self):
        ''' Return the unary vectors of the local axes.'''
        return [self.refSys.getIVector(), self.refSys.getJVector(), self.refSys.getKVector()]
    

def getMidPlane(shape, tol= 1e-6):
    ''' Return the mid-plane of the shape argument.

    :param shape : TopoDS_Shape or a subclass such as TopoDS_Face
                   the shape to compute the bounding box from.
    :param tol: tolerance of the computed boundingbox.
    '''
    shapeFootprint= GeometryFootprint(shape)
    centroidPos= shapeFootprint.getCentroid()
    shapeVectors= shapeFootprint.getAxisDirections()
    shapeDimensions= shapeFootprint.shapeDimensions

    # Compute index of minimum dimension
    thickness, idx = min((thickness, idx) for (idx, thickness) in enumerate(shapeDimensions))
    extension= max(shapeDimensions) # Compute maximum dimension.
    normalVector= shapeVectors[idx]
    midPlane= OCC.Core.gp.gp_Pln( OCC.Core.gp.gp_Pnt(centroidPos.x, centroidPos.y, centroidPos.z), OCC.Core.gp.gp_Dir(normalVector.x, normalVector.y, normalVector.z) )
    return midPlane, thickness, extension

def extractEdges(section):
    ''' Return a edges of the section argument.'''
    exp = OCC.Core.TopExp.TopExp_Explorer(section, OCC.Core.TopAbs.TopAbs_EDGE)
    retval= list()
    while exp.More():
        retval.append(OCC.Core.TopoDS.topods.Edge(exp.Current()))
        exp.Next()
    return retval

def extractWires(section):
    ''' Return a wire from the vertices of the section argument.'''
    retval= list()
    section_edges= extractEdges(section)
    if len(section_edges) > 0:
        segments= list()
        for edge in section_edges:
            vertices= extractVertices(edge)
            x0, y0, z0= BRep_Tool().Pnt(vertices[0]).Coord()
            pt0= geom.Pos3d(x0, y0, z0)
            x1, y1, z1= BRep_Tool().Pnt(vertices[1]).Coord()
            pt1= geom.Pos3d(x1, y1, z1)
            sg= geom.Segment3d(pt0, pt1)
            segments.append(sg)

        # A wire is formed by connecting the edges
        polylines= geom.get_3d_polylines(segments,.001)
        retval= polylines
    return retval

def extractVertices(edge):
    ''' Return the vertices of the edge argument.'''
    exp = OCC.Core.TopExp.TopExp_Explorer(edge, OCC.Core.TopAbs.TopAbs_VERTEX)
    retval= list()
    while exp.More():
        retval.append(OCC.Core.TopoDS.topods.Vertex(exp.Current()))
        exp.Next()
    return retval

def getMaterials(ifcProduct):
    ''' Return the materials associated with an IFC product.'''
    retval= list()
    if ifcProduct.HasAssociations:
        for i in ifcProduct.HasAssociations:
             if i.RelatingMaterial.is_a('IfcMaterial'):
                retval.append(i.RelatingMaterial.Name)
             if i.RelatingMaterial.is_a('IfcMaterialList'):
                 for materials in i.RelatingMaterial.Materials:
                     retval.append(materials.Name)
             if i.RelatingMaterial.is_a('IfcMaterialLayerSetUsage'):
                 for materials in i.RelatingMaterial.ForLayerSet.MaterialLayers:
                     material= materials.Material
                     materialName= None
                     if(material):
                         materialName= material.Name
                     thickness= materials.LayerThickness
                     retval.append({'materialName':materialName, 'thickness':thickness})
    return retval
    
# Specify to return pythonOCC shapes from ifcopenshell.geom.create_shape()
settings= ifcopenshell.geom.settings()
settings.set(settings.USE_PYTHON_OPENCASCADE, True)

def getProductShapes(ifcModel):
    ''' For every product a shape is created if the shape has a Representation.

    :param ifcModel: IFC model to get the shapes from. 
    '''
    products= ifcModel.by_type("IfcProduct")
    retval= dict()

    # For every product a shape is created if the shape has a Representation.
    for product in products:
        if product.is_a("IfcOpeningElement") or product.is_a("IfcSite"):
            continue
        if product.Representation is not None:
            print(product.Name, product.is_a())
            shape= ifcopenshell.geom.create_shape(settings, product).geometry
            materialList= getMaterials(product)
            retval[product.id()]= {'product':product, 'shape':shape, 'materials': materialList}
    return retval

shellIFCTypes= ['IfcCurtainWall', 'IfcFooting', 'IfcOpeningElement', 'IfcRailing', 'IfcRamp', 'IfcRampFlight', 'IfcReinforcingMesh', 'IfcReinforcingMesh', 'IfcSlab', 'IfcSlab', 'IfcStair', 'IfcStairFlight', 'IfcWall', 'IfcWallStandardCase']

beamColumnIFCTypes= ['IfcBeam', 'IfcColumn', 'IfcPile', 'IfcReinforcingBar', 'IfcTendon']

def computeMidSurfaces(productShapes):
    ''' For each thin walled shape in product shapes compute the contour of 
        its mid-surface.

    :param productShapes: dictionary containing the shapes to process.
    '''
    retval= dict()
    for key in productShapes:
        productData= productShapes[key]
        product= productData['product']
        ifcType= product.is_a()
        shape= productData['shape']
        materialList= productData['materials']
        midPlane, thickness, extension= getMidPlane(shape)
        midSectionFace= OCC.Core.BRepBuilderAPI.BRepBuilderAPI_MakeFace(midPlane, -extension, extension, -extension, extension).Face()
        shapeSection= OCC.Core.BRepAlgoAPI.BRepAlgoAPI_Section(midSectionFace, shape).Shape()
        midSectionWires= extractWires(shapeSection)
        retval[product.id()]= {'product':product, 'shape':shape, 'midPlane':midPlane, 'midSection':midSectionFace, 'midSectionWires':midSectionWires, 'thickness':thickness,'materials':materialList, 'ifcType': ifcType}

    return retval



