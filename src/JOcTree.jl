module JOcTree

    importall jInv.Mesh
    using jInv.Utils
    export OcTreeMesh
    abstract OcTreeMesh <: AbstractMesh

    include("sparse3.jl")
    include("OcTreeMeshFV.jl")
    include("display.jl")

    include("getNumbering.jl")
    include("getGrids.jl")
    include("getEdgeSize.jl")
    include("getFaceSize.jl")

    include("getCurlMatrixRec.jl")
    include("getDivergenceMatrixRec.jl")
    include("getNodalGradientRec.jl")

    include("getEdgeMassMatrix.jl")
    include("getEdgeMassMatrixAnisotropic.jl")
    include("getEdgeMassMatrixAnisotropicNoWeight.jl")
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
    include("regularizeOcTree2.jl")

    include("getNumberOfNeighbors.jl")
    include("getVolume.jl")
    include("getLength.jl")

    include("getEdgeIntegralOfPolygonalChain.jl")
    
    include("refineOcTree.jl")    
    include("splitCells.jl")
    include("uniteOcTrees.jl")

    include("getInterfaceWeights.jl")

    # Mesh Creation
    include("createOcTree/initCoarseOcTree.jl")
    include("createOcTree/createOcTreeFromBox.jl")
    include("createOcTree/createOcTreeFromImage.jl")
    # include("createOcTree/createOcTreeFromPoints.jl")
    # include("createOcTree/createOcTreeFromTopo.jl")
    include("createOcTree/createOcTreeMesh.jl")
    include("createOcTree/OctreeBoxPolygon.jl")

    include("IO/importUBC.jl")
    include("IO/exportUBC.jl")

end