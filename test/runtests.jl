using JOcTree
using Base.Test

include("randomOctreeMesh.jl") 

@testset "JOcTree" begin
include("testIO.jl")
include("testMatrices.jl")
include("testMassMatrices.jl")
include("testMassMatrixDerivatives.jl")
include("testInterpolationMatrix.jl")
end