using JOcTree
using Base.Test

include("randomOctreeMesh.jl")
intTypes = (Int32,Int64)
@testset "JOcTree" begin
    include("testIO.jl")
    include("testMatrices.jl")
    include("testMassMatrices.jl")
    include("testMassMatrixDerivatives.jl")
    include("testInterpolationMatrix.jl")
    include("testConstraints.jl")
end
