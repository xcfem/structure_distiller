# -*- coding: utf-8 -*-
''' Compute wall mid-plane and beam-column axes for use in structural 
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

def getMidPlane(shape, tol= 1e-6):
    ''' Return the mid-plane of the shape argument.

    :param shape : TopoDS_Shape or a subclass such as TopoDS_Face
                   the shape to compute the bounding box from.
    :param tol: tolerance of the computed boundingbox.
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
    v0= geom.Vector3d(matrix_of_inertia[0][0], matrix_of_inertia[0][1], matrix_of_inertia[0][2])
    v1= geom.Vector3d(matrix_of_inertia[1][0], matrix_of_inertia[1][1], matrix_of_inertia[1][2])
    v2= geom.Vector3d(matrix_of_inertia[2][0], matrix_of_inertia[2][1], matrix_of_inertia[2][2])
    shapeVectors= [v0, v1, v2]

    ref= geom.Ref3d3d(centroidPos, centroidPos+v0, centroidPos+v1)

    box= getBoundingBox(shape)
    pMin= geom.Pos3d(box[0][0], box[0][1], box[0][2])
    pMax= geom.Pos3d(box[1][0], box[1][1], box[1][2])

    localMin= ref.getLocalPosition(pMin)
    localMax= ref.getLocalPosition(pMax)
    shapeDimensions= [localMax.x-localMin.x, localMax.y-localMin.y, localMax.z-localMin.z]

    # Compute index of minimum dimension
    thickness, idx = min((thickness, idx) for (idx, thickness) in enumerate(shapeDimensions))
    normalVector= shapeVectors[idx]
    midPlane= OCC.Core.gp.gp_Pln( OCC.Core.gp.gp_Pnt(centroidPos.x, centroidPos.y, centroidPos.z), OCC.Core.gp.gp_Dir(normalVector.x, normalVector.y, normalVector.z) )
    return midPlane, thickness

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

def extractVertices(section):
    ''' Return the vertices of the section argument.'''
    exp = OCC.Core.TopExp.TopExp_Explorer(section, OCC.Core.TopAbs.TopAbs_VERTEX)
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
                     print('materials: ', materials)
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

def computeMidSurfaces(productShapes):
    ''' For each thin walled shape in product shapes compute the contour of 
        its mid-surface.

    :param productShapes: dictionary containing the shapes to process.
    '''
    retval= dict()
    for key in productShapes:
        productData= productShapes[key]
        product= productData['product']
        shape= productData['shape']
        materialList= productData['materials']
        midPlane, thickness= getMidPlane(shape)
        midSection= OCC.Core.BRepBuilderAPI.BRepBuilderAPI_MakeFace(midPlane, -10, 10, -10, 10).Face()
        midSectionWires= extractWires(midSection)
        retval[product.id()]= {'product':product, 'shape':shape, 'midPlane':midPlane, 'midSection':midSection, 'midSectionWires':midSectionWires, 'thickness':thickness,'materials':materialList}

    return retval



