# -*- coding: utf-8 -*-
''' Compute shell mid-surface and beam-column axes for use in structural 
    analytical model.'''

__author__= "Luis C. Pérez Tato (LCPT) and Ana Ortega (AOO)"
__copyright__= "Copyright 2021, LCPT and AOO"
__license__= "GPL"
__version__= "3.0"
__email__= "l.pereztato@ciccp.es"

import sys
import math

import ifcopenshell
import ifcopenshell.geom

from misc_utils import log_messages as lmsg

import OCC.Core.gp
import OCC.Core.BRepBuilderAPI
import OCC.Core.BRepAlgoAPI
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakePrism
import OCC.Core.TopExp
import OCC.Core.ShapeAnalysis
import OCC.Core.TopoDS
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import topods_Face
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop_VolumeProperties
from OCC.Core.BRepBndLib import brepbndlib_Add, brepbndlib_AddOptimal, brepbndlib_AddOBB
from OCC.Core.Bnd import Bnd_Box, Bnd_OBB
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.Core.BRep import BRep_Tool

import xc_base
import geom
import numpy as np

def getBoundingBox(shape, tol:float =1e-6) -> tuple:
    ''' return the bounding box of the TopoDS_Shape `shape`


    :param shape : TopoDS_Shape or a subclass such as TopoDS_Face
                   the shape to compute the bounding box from
    :param tol: tolerance of the computed boundingbox
    '''
    bbox= Bnd_Box()
    bbox.SetGap(tol)
    brepbndlib_Add(shape, bbox, False)

    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    return ((xmin, ymin, zmin), (xmax, ymax, zmax))

def getOBB(shape, optimal_OBB=True):
    '''return the oriented bounding box of the TopoDS_Shape `shape`

    Parameters
    ----------

    shape : TopoDS_Shape or a subclass such as TopoDS_Face
        the shape to compute the bounding box from
    optimal_OBB : bool, True by default. If set to True, compute the
        optimal (i.e. the smallest oriented bounding box). Optimal OBB is
        a bit longer.
    Returns
    -------
        the oriented bounding box.
    '''
    obb = Bnd_OBB()
    if optimal_OBB:
        is_triangulation_used = True
        is_optimal = True
        is_shape_tolerance_used = False
        brepbndlib_AddOBB(
            shape, obb, is_triangulation_used, is_optimal, is_shape_tolerance_used
        )
    else:
        brepbndlib_AddOBB(shape, obb)

    return obb

def getOrientation(shape, optimal_OBB=True):
    '''return the vectors orientation vectors of the oriented bounding box of the TopoDS_Shape `shape`

    Parameters
    ----------

    shape : TopoDS_Shape or a subclass such as TopoDS_Face
        the shape to compute the bounding box from
    optimal_OBB : bool, True by default. If set to True, compute the
        optimal (i.e. the smallest oriented bounding box). Optimal OBB is
        a bit longer.
    Returns
    -------
        the baricentre and three vectors.
    '''
    obb = getOBB(shape, optimal_OBB)

    # extract the required values
    cog = obb.Center()
    x_direction = obb.XDirection()
    y_direction = obb.YDirection()
    z_direction = obb.ZDirection()
    vx= geom.Vector3d(x_direction.X(), x_direction.Y(), x_direction.Z())
    vy= geom.Vector3d(y_direction.X(), y_direction.Y(), y_direction.Z())
    vz= geom.Vector3d(z_direction.X(), z_direction.Y(), z_direction.Z())
    bary_center= geom.Pos3d(cog.X(), cog.Y(), cog.Z())
    return bary_center, [vx, vy, vz]
    

def getOrientedBoundingbox(shape, optimal_OBB=True):
    '''return the oriented bounding box of the TopoDS_Shape `shape`

    Parameters
    ----------

    shape : TopoDS_Shape or a subclass such as TopoDS_Face
        the shape to compute the bounding box from
    optimal_OBB : bool, True by default. If set to True, compute the
        optimal (i.e. the smallest oriented bounding box). Optimal OBB is
        a bit longer.
    Returns
    -------
        a list with center, x, y and z sizes

        a shape
    '''
    obb= getOBB(shape, optimal_OBB)

    # converts the bounding box to a shape
    bary_center = obb.Center()
    x_direction = obb.XDirection()
    y_direction = obb.YDirection()
    z_direction = obb.ZDirection()
    a_half_x = obb.XHSize()
    a_half_y = obb.YHSize()
    a_half_z = obb.ZHSize()

    ax = OCC.Core.gp.gp_XYZ(x_direction.X(), x_direction.Y(), x_direction.Z())
    ay = OCC.Core.gp.gp_XYZ(y_direction.X(), y_direction.Y(), y_direction.Z())
    az = OCC.Core.gp.gp_XYZ(z_direction.X(), z_direction.Y(), z_direction.Z())
    p = OCC.Core.gp.gp_Pnt(bary_center.X(), bary_center.Y(), bary_center.Z())
    an_axe = OCC.Core.gp.gp_Ax2(p, OCC.Core.gp.gp_Dir(z_direction), OCC.Core.gp.gp_Dir(x_direction))
    an_axe.SetLocation(OCC.Core.gp.gp_Pnt(p.XYZ() - ax * a_half_x - ay * a_half_y - az * a_half_z))
    a_box = BRepPrimAPI_MakeBox(
        an_axe, 2.0 * a_half_x, 2.0 * a_half_y, 2.0 * a_half_z
    ).Shape()
    return bary_center, [a_half_x, a_half_y, a_half_z], a_box


def getShapeFaces(shape):
    ''' Return the faces of the shape argument.

    :param shape : TopoDS_Shape or a subclass such as TopoDS_Face
                   the shape to compute the bounding box from
    '''
    retval= list()
    topExp = TopExp_Explorer()
    topExp.Init(shape, TopAbs_FACE)

    while topExp.More():
        fc = topods_Face(topExp.Current())
        retval.append(fc)
        topExp.Next()

    return retval

def getAngle(vA, vB, tol= .1):
    retval= vA.getAngle(vB)
    if(abs(retval-math.pi)<tol):
        retval= 0.0
    if(abs(retval-2*math.pi)<tol):
        retval= 0.0
    return retval

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

    # Get inertia tensor.
    _m = props.MatrixOfInertia()
    matrix_of_inertia = [
         (_m.Row(row + 1).X(), _m.Row(row + 1).Y(), _m.Row(row + 1).Z())
         for row in range(3)
     ]
    eigenvalues, eigenvectors = np.linalg.eig(np.array(matrix_of_inertia))
    sortedIndexes= np.argsort(eigenvalues)
    i0= sortedIndexes[0]; i1= sortedIndexes[1]; i2= sortedIndexes[2]
    v0= geom.Vector3d(eigenvectors[i0][0], eigenvectors[i0][1], eigenvectors[i0][2])
    v1= geom.Vector3d(eigenvectors[i1][0], eigenvectors[i1][1], eigenvectors[i1][2])
    v2= geom.Vector3d(eigenvectors[i2][0], eigenvectors[i2][1], eigenvectors[i2][2])

    # Get oriented bounding box.
    bary_center, o_vectors= getOrientation(shape)
    [vx, vy, vz]= o_vectors
    ang_v0= [getAngle(v0, vx), getAngle(v0, vy), getAngle(v0, vz)]
    index_v0= np.argmin(ang_v0)
    ang_v1= [getAngle(v1, vx), getAngle(v1, vy), getAngle(v1, vz)]
    index_v1= np.argmin(ang_v1)
    ang_v2= [getAngle(v2, vx), getAngle(v2, vy), getAngle(v2, vz)]
    index_v2= np.argmin(ang_v2)
    if(index_v0==index_v2):
        print(o_vectors)
        print('v0: ', v0, o_vectors[index_v0], ang_v0, index_v0)
        print('v2: ', v2, o_vectors[index_v2], ang_v2, index_v2)
        index_v2= index_v1
    v0= o_vectors[index_v0]
    v2= o_vectors[index_v2]
    centroidPos= bary_center
    
    return geom.Ref3d3d(centroidPos, centroidPos+v0, centroidPos+v2)

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
    ''' characteristic geometric parameters obtained from a shape
        local axes, dimensions and dimensionality.

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
    

def getShapeMidPlane(shape, shapeFootprint, tol= 1e-6):
    ''' Return the mid-plane of the shape argument.

    :param shape : TopoDS_Shape or a subclass such as TopoDS_Face
                   the shape to compute the bounding box from.
    :param shapeFootprint: characteristic geometric parameters obtained 
                           the shape arguments.
    :param tol: tolerance of the computed boundingbox.
    '''
    centroidPos= shapeFootprint.getCentroid()
    shapeVectors= shapeFootprint.getAxisDirections()
    shapeDimensions= shapeFootprint.shapeDimensions

    # Compute index of minimum dimension
    thickness, idx = min((thickness, idx) for (idx, thickness) in enumerate(shapeDimensions))
    extension= max(shapeDimensions) # Compute maximum dimension.
    #normalVector= shapeVectors[idx]
    normalVector= shapeVectors[1]
    midPlane= OCC.Core.gp.gp_Pln( OCC.Core.gp.gp_Pnt(centroidPos.x, centroidPos.y, centroidPos.z), OCC.Core.gp.gp_Dir(normalVector.x, normalVector.y, normalVector.z) )
    planeNormal= midPlane.Axis().Direction().XYZ()
    return midPlane, thickness, extension

def getShapeAxis(shape, shapeFootprint):
    ''' Return the axis of the shape argument.

    :param shape : TopoDS_Shape or a subclass such as TopoDS_Face
                   the shape to compute the bounding box from.
    :param tol: tolerance of the computed boundingbox.
    '''
    centroidPos= shapeFootprint.getCentroid()
    shapeVectors= shapeFootprint.getAxisDirections()
    shapeDimensions= shapeFootprint.shapeDimensions
    length, idx = max((length, idx) for (idx, length) in enumerate(shapeDimensions))
    auxVector= shapeVectors[idx]*(length/2.0)
    fromPoint= centroidPos-auxVector
    toPoint= centroidPos+auxVector
    axis= geom.Segment3d(fromPoint, toPoint)
    return axis, length

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
            shape= ifcopenshell.geom.create_shape(settings, product).geometry
            materialList= getMaterials(product)
            retval[product.id()]= {'product':product, 'shape':shape, 'materials': materialList}
    return retval

shellIFCTypes= ['IfcCurtainWall', 'IfcFooting', 'IfcOpeningElement', 'IfcRailing', 'IfcRamp', 'IfcRampFlight', 'IfcReinforcingMesh', 'IfcReinforcingMesh', 'IfcSlab', 'IfcPlate', 'IfcStair', 'IfcStairFlight', 'IfcWall', 'IfcWallStandardCase']

beamColumnIFCTypes= ['IfcBeam','IfcBeamStandardCase', 'IfcColumn', 'IfcPile', 'IfcReinforcingBar', 'IfcTendon']

def computeShapesDatum(productShapes):
    ''' For each thin walled shape in product shapes compute the contour of 
        its mid-surface and for each linear shape computes its axis.

    :param productShapes: dictionary containing the shapes to process.
    '''
    surfaces= dict()
    lines= dict()
    for key in productShapes:
        productData= productShapes[key]
        product= productData['product']
        ifcType= product.is_a()
        shape= productData['shape']
        materialList= productData['materials']
        shapeFootprint= GeometryFootprint(shape)
        if(ifcType in shellIFCTypes):
            midPlane, thickness, extension= getShapeMidPlane(shape, shapeFootprint)
            planeNormal= midPlane.Axis().Direction().XYZ()
            midSectionFace= OCC.Core.BRepBuilderAPI.BRepBuilderAPI_MakeFace(midPlane, -extension, extension, -extension, extension).Face()
            shapeSection= OCC.Core.BRepAlgoAPI.BRepAlgoAPI_Section(midSectionFace, shape).Shape()
            midSectionWires= extractWires(shapeSection)
            surfaces[product.id()]= {'product':product, 'shape':shape, 'footprint':shapeFootprint, 'materials':materialList, 'ifcType': ifcType, 'midPlane':midPlane, 'midSection':midSectionFace, 'midSectionWires':midSectionWires, 'thickness':thickness}
        elif(ifcType in beamColumnIFCTypes):
            axis, length= getShapeAxis(shape, shapeFootprint)
            lines[product.id()]= {'product':product, 'shape':shape, 'footprint':shapeFootprint, 'materials':materialList, 'ifcType': ifcType, 'axis':axis, 'length':length}
        else:
            funcName= sys._getframe(0).f_code.co_name
            lmsg.error(funcName+'; IFC type: '+ifcType+' ignored.')

    return surfaces, lines



