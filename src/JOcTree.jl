module JOcTree

    importall jInv.Mesh
    using jInv.Utils
    export OcTreeMesh
    abstract type OcTreeMesh <: AbstractMesh end

    include("sparse3.jl")
    include("OcTreeMeshFV.jl")
    include("display.jl")

    include("getGrids.jl")
    include("getSizeNumbering.jl")
    include("getSizeOfNeighbors.jl")

    include("getCurlMatrixRec.jl")
    include("getDivergenceMatrixRec.jl")
    include("getNodalGradientRec.jl")
    include("getCellCenterGradientMatrix.jl")

    include("getEdgeMassMatrix.jl")
    include("getFaceMassMatrix.jl")
    include("getNodalMassMatrix.jl")

    include("getEdgeInterpolationMatrix.jl")
    include("getFaceInterpolationMatrix.jl")
    include("getNodalInterpolationMatrix.jl")

    include("getInterpolationMatrix.jl")
    include("getEdgeToCellCenteredMatrix.jl")
    include("getFaceToCellCenteredMatrix.jl")
    include("getNodesToCellCenteredMatrix.jl")

    include("findNonRegularBlocks.jl")
    include("findBlocks.jl")

    include("getNodalConstraints.jl")
    include("getEdgeConstraints.jl")
    include("getFaceConstraints.jl")

    include("regularizeOcTree.jl")
    include("getNumberOfNeighbors.jl")
    include("refineOcTree.jl")
    include("splitCells.jl")
    include("uniteOcTrees.jl")

    include("Utils.jl")

    include("getEdgeIntegralOfPolygonalChain.jl")

    include("getInterfaceWeights.jl")

    # Mesh Creation
    include("createOcTree/createOcTreeFromBox.jl")
    include("createOcTree/createOcTreeFromImage.jl")
    include("createOcTree/initializeOctree.jl")

    include("IO/importUBC.jl")
    include("IO/exportUBC.jl")

end
