export plot

if Pkg.installed("MATLAB") == Nothing()
	
	function missingPackage()
		
		warn("Install Julia package MATLAB to enable OcTree mesh plotting.")
		return
		
	end
	
	plot(mesh::OcTreeMesh) = missingPackage()
	plot(mesh::OcTreeMesh, u::Array{Int64,1}) = missingPackage()
	plot(mesh::OcTreeMesh, u::Array{Float64,1}) = missingPackage()
	
else
	
	using MATLAB
	
	function plot(mesh::OcTreeMesh; f = 1)
		
		u = iround(log2(nonzeros(mesh.S)) .+ 1)
		
		plot(mesh, u; f = f)
		
		return
	
	end

	function plot(mesh::OcTreeMesh, u::Array{Int64,1}; f = 1)
	
		# to avoid incompatibilities between Julia's and Matlab's sparse3 we pass the indices and values
		i,j,k,bsz = find3(mesh.S)
		n         = mesh.n
		h         = mesh.h
		x0        = mesh.x0
	
		# set path to Matlab function plotOcTree.m which resides in the same directory like this Julia function
		(dname,fname) = splitdir(@__FILE__())
		mxcall(:addpath, 0, dname)
	
		# plot
		mxcall(:plotOcTreeMesh, 0, f, i, j, k, bsz, n, h, x0, u)
	
		return
	
	end

	function plot(mesh::OcTreeMesh, u::Array{Float64,1}; f = 1)
	
		# to avoid incompatibilities between Julia's and Matlab's sparse3 we pass the indices and values
		i,j,k,bsz = find3(mesh.S)
		n         = mesh.n
		h         = mesh.h
		x0        = mesh.x0
	
		# set path to Matlab function plotOcTree.m which resides in the same directory like this Julia function
		(dname,fname) = splitdir(@__FILE__())
		mxcall(:addpath, 0, dname)
	
		# plot
		mxcall(:plotOcTreeMesh, 0, f, i, j, k, bsz, n, h, x0, u)
	
		return
	
	end
	
end