# -*- coding: utf-8 -*-
''' Compute wall mid-plane for use in structural analytical model.'''

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

def getBoundingBox(shape, tol=1e-6) -> tuple:
    """ return the bounding box of the TopoDS_Shape `shape`
    Parameters
    ----------
    shape : TopoDS_Shape or a subclass such as TopoDS_Face
        the shape to compute the bounding box from
    tol: float
        tolerance of the computed boundingbox
    """
    bbox= Bnd_Box()
    bbox.SetGap(tol)
    brepbndlib_Add(shape, bbox, False)

    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    return ((xmin, ymin, zmin), (xmax, ymax, zmax))

def getMidPlane(shape, tol= 1e-6):
    ''' Return the mid-plane of the shape argument.

    shape : TopoDS_Shape or a subclass such as TopoDS_Face
            the shape to compute the bounding box from
    tol: tolerance of the computed boundingbox
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
    val, idx = min((val, idx) for (idx, val) in enumerate(shapeDimensions))
    normalVector= shapeVectors[idx]
    midPlane= OCC.Core.gp.gp_Pln( OCC.Core.gp.gp_Pnt(centroidPos.x, centroidPos.y, centroidPos.z), OCC.Core.gp.gp_Dir(normalVector.x, normalVector.y, normalVector.z) )

    print('shape vectors: ', shapeVectors)
    print('shape dimensions: ', shapeDimensions)
    print('minimum dimension: ', val, ' at index:', idx)
    print('normal vector: ', shapeVectors[idx])
    return midPlane

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
            print(pt0, pt1)
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

def createStructuralContext(ifcfile):

    """Creates an additional geometry context for structural objects. Returns the new context"""

    contexts = ifcfile.by_type("IfcGeometricRepresentationContext")
    # filter out subcontexts
    contexts = [c for c in contexts if c.is_a() == "IfcGeometricRepresentationContext"]
    geomContext = contexts[0] # arbitrarily take the first one...
    structContext = ifcfile.createIfcGeometricRepresentationSubContext(
        'Analysis', 'Axis', None, None, None, None, geomContext, None, "GRAPH_VIEW", None)
    return structContext

def createFaceSurface(ifcModel, wires):
    ''' Create an IfcFaceSurface object from the wires argument.

    :param ifcModel: IFC model where the IfcFaceSurface will be created.
    :param wires: polygons defining the surface and its holes (if any).
    '''
    sz= len(wires)
    faceBounds= list()
    if(sz>0): # Outer bound.
        outerBoundWire= wires[0]
        ifcPoints= list()
        for p in outerBoundWire[:-1]:
            ifcPoints.append(ifcModel.createIfcCartesianPoint((p.x, p.y, p.z)))
        polyLoop= ifcModel.createIfcPolyLoop(ifcPoints)
        outerBound= ifcModel.createIfcFaceOuterBound(Bound=polyLoop, Orientation=True)  # orientation of vertices is CCW
        faceBounds.append(outerBound)
    if(sz>1): # Inner bounds.
        for interiorWire in wires[1:]:
            ifcPoints= list()
            for p in interiorWire[:-1]:
                ifcPoints.append(ifcModel.createIfcCartesianPoint((p.x, p.y, p.z)))
            polyLoop= ifcModel.createIfcPolyLoop(ifcPoints)
            innerBound= ifcModel.createIfcFaceBound(Bound=polyLoop, Orientation=False) # orientation of vertices is CW
            faceBounds.append(innerBound)
    return ifcModel.createIfcFaceSurface(Bounds= faceBounds)
            

def createGlobalAxes(ifcModel):
    xAxis= ifcModel.createIfcDirection((1.0, 0.0, 0.0))
    yAxis= ifcModel.createIfcDirection((0.0, 1.0, 0.0))
    zAxis= ifcModel.createIfcDirection((0.0, 0.0, 1.0))
    origin= ifcModel.createIfcCartesianPoint((0.0, 0.0, 0.0))
    axes= ifcModel.createIfcAxis2Placement3D(origin, zAxis, xAxis)
    return axes
        
# Open the IFC file using IfcOpenShell
import os
pth= os.path.dirname(__file__)
# print("pth= ", pth)
if(not pth):
  pth= "."
fName= pth+"/../data/Wall.ifc"
inputFile= ifcopenshell.open(fName)

# Specify to return pythonOCC shapes from ifcopenshell.geom.create_shape()
settings= ifcopenshell.geom.settings()
settings.set(settings.USE_PYTHON_OPENCASCADE, True)

# The geometric elements in an IFC file are the IfcProduct elements. So
# these are stored in product_shapes.

products= inputFile.by_type("IfcProduct")
product_shapes= dict()

# For every product a shape is created if the shape has a Representation.
for product in products:
    if product.is_a("IfcOpeningElement") or product.is_a("IfcSite"): continue
    if product.Representation is not None:
        shape= ifcopenshell.geom.create_shape(settings, product).geometry
        midPlane= getMidPlane(shape)
        midSection= OCC.Core.BRepBuilderAPI.BRepBuilderAPI_MakeFace(midPlane, -10, 10, -10, 10).Face()
        midSectionWires= extractWires(midSection)
        product_shapes[product.id()]= {'product':product, 'shape':shape, 'midPlane':midPlane, 'midSection':midSection, 'midSectionWires':midSectionWires}

# Write ouput
outputFile= ifcopenshell.file(schema=inputFile.schema)
ifcProject= inputFile.by_type("IfcProject")[0]
outputFile.add(ifcProject)
uid = ifcopenshell.guid.new
owhList= outputFile.by_type("IfcOwnerHistory")
if(len(owhList)>0):
    owh= owhList[0]
else:
    owh= outputFile.createIfcOwnerHistory()
globalAxes= createGlobalAxes(outputFile)
localPlacement= outputFile.createIfcLocalPlacement(None, globalAxes)
structContext = createStructuralContext(outputFile)
mod = outputFile.createIfcStructuralAnalysisModel(uid(),owh,"Structural Analysis Model",None,None,"NOTDEFINED",None,None,None,localPlacement)

for shapeKey in product_shapes:
    print(shapeKey)
    shapeData= product_shapes[shapeKey]
    wires= shapeData['midSectionWires']
    ifcFaceSurface= createFaceSurface(outputFile, wires)
    uid = ifcopenshell.guid.new()
    label= str(shapeKey)
    topologyRep = outputFile.createIfcTopologyRepresentation(structContext, "Analysis", "Face", (ifcFaceSurface,))
    prodDefShape = outputFile.createIfcProductDefinitionShape(None, None, (topologyRep,))
    thickness= 0.1
    outputFile.createIfcStructuralSurfaceMember(uid, owh, label, None, None, localPlacement, prodDefShape, "SHELL", thickness)

outputFileName= '/tmp/test.ifc'
outputFile.write(outputFileName)



