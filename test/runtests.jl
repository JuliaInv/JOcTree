
using JOcTree
using Base.Test

@testset "JOcTree" begin
println("==== test input & output ====")
include("testIO.jl")
include("testMatrices.jl")
include("testMassMatrixDerivatives.jl")
end