module JOcTree

include("sparse3.jl")

# using jInv.Mesh.AbstractMesh
# using jInv.Mesh.ndgrid
importall jInv.Mesh
using jInv.Utils
export OcTreeMesh
abstract OcTreeMesh <: AbstractMesh


include("OcTreeMeshFV.jl")
include("OcTreeMeshFEM.jl")

include("findNonRegularBlocks.jl")
include("findBlocks.jl")
include("getCellNumbering.jl")
include("getCellCenteredGrid.jl")
include("getCurlMatrixRec.jl")
include("getDivergenceMatrixRec.jl")
include("getEdgeToCellCenteredMatrix.jl")
include("getEdgeMassMatrix.jl")
include("getEdgeMassMatrixAnisotropic.jl")
include("getEdgeNumbering.jl")
include("getEdgeGrids.jl")
include("getEdgeSize.jl")
include("getEdgeInterpolationMatrix.jl")
include("getFaceToCellCenteredMatrix.jl")
include("getFaceMassMatrix.jl")
include("getFaceGrids.jl")
include("getFaceSize.jl")
include("getFaceNumbering.jl")
include("getFaceInterpolationMatrix.jl")
include("getInterpolationMatrix.jl")
include("getNodalConstraints.jl")
include("getFaceConstraints.jl")
include("getNodalMassMatrix.jl")
include("getNodalInterpolationMatrix.jl")
include("getNodalNumbering.jl")
include("getNodalGradientRec.jl")
include("getNodalGrid.jl")
include("getNodesToCellCenteredMatrix.jl")
include("getNumberOfNeighbors.jl")
include("regularizeOcTree.jl")
include("refineOcTree.jl")
include("getVolume.jl")
include("getLength.jl")
include("uniteOcTrees.jl")

include("createOcTreeFromBox.jl")
include("getEdgeIntegralOfPolygonalChain.jl")
include("createOcTreeFromImage.jl")
include("initCoarseOcTree.jl")

include("display.jl")
include("sparse3.jl")

include("getEdgeConstraints.jl")
include("createOcTreeMesh.jl")
include("OctreeBoxPolygon.jl")

#include("createSmallMeshFromTX.jl")
#include("createOcTreeFromTRX.jl")

include("regularizeOcTree2.jl")
include("splitCells.jl")
include("getInterfaceWeights.jl")

#include("getLocalElementMatrices.jl") 
#include("getMassMatrixFEM.jl") 
#include("getDiffMassMatrixFEM.jl")

# include("plot.jl")


end