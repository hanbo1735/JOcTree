#module JOcTree


export OcTreeMeshFV, getOcTreeMeshFV


"""
    struct MassMatrix

Internal storage for mass matrix integration
 - `M = getNodalMassMatrix(mesh, sigma)`
 - `M = getEdgeMassMatrix(mesh, sigma)`
 - `M = getFaceMassMatrix(mesh, sigma)`

Fields
 - `n::Int64`: size of mass matrix `M` (number of nodes, edges, faces in `mesh`)
 - `A::SparseMatrixCSC{Float64,Int64}`: map coefficient vector `sigma` to
      nonzero entries in mass matrix `M` (numerical integration)
 - `rowval::Array{Int64,1}`: row indices of nonzero entries in `M`
 - `colptr::Array{Int64,1}`: CSC format column pointers in `M`
 - `colval::Array{Int64,1}`: column indices of nonzero entries in `M`
"""
struct MassMatrix
    n::Int64
    A::SparseMatrixCSC{Float64,Int64}
    rowval::Array{Int64,1}
    colptr::Array{Int64,1}
    colval::Array{Int64,1}
end


"""
    MassMatrix()

Default constructor for `MassMatrix`
"""
function MassMatrix()
    F = Array{Float64}(0)
    I = Array{Int64}(0)
    S = SparseMatrixCSC(0, 0, [1], I, F)
    return MassMatrix(0, S, I, I, I)
end


mutable struct OcTreeMeshFV <: OcTreeMesh
	S::SparseArray3D    # i,j,k, bsz
	h::Vector{Float64}  # (3) underlying cell size
	x0::Vector{Float64} # coordinates of corner of mesh
	n::Vector{Int64} # underlying mesh
	nc::Int          # number of cells
	nf::Vector{Int}  # (3) number of faces
	ne::Vector{Int}  # (3) number of edges
  nn::Int          # number of nodes
	Div::SparseMatrixCSC
	Grad::SparseMatrixCSC
	Curl::SparseMatrixCSC
	Pf::Dict{Int64,MassMatrix} # face mass matrix integration storage
	Pe::Dict{Int64,MassMatrix} # edge mass matrix integration storage
	Pn::Dict{Int64,MassMatrix} # nodal mass matrix integration storage
	Af::SparseMatrixCSC # Face to cell-centre matrix
	Ae::SparseMatrixCSC # Edge to cell-centre matrix
	An::SparseMatrixCSC # Node to cell-centre matrix
	V::SparseMatrixCSC # cell volume
	L::SparseMatrixCSC # edge lengths
	Ne::SparseMatrixCSC # Edge nullspace matrix
	Qe::SparseMatrixCSC # Edge projection matrix
	activeEdges::Vector{Int64}   # lookup table for new edge enumeration
	activeFaces::Vector{Int64}   # lookup table for new face enumeration
	activeNodes::Vector{Int64}   # lookup table for new node enumeration
	Nn::SparseMatrixCSC # Node nullspace matrix
	Qn::SparseMatrixCSC # Node projection matrix
	Nf::SparseMatrixCSC # Face nullspace matrix
	Qf::SparseMatrixCSC # Face projection matrix
	FX::SparseArray3D  # X face size
	FY::SparseArray3D  # Y face size
	FZ::SparseArray3D  # Z face size
	EX::SparseArray3D  # X edge size
	EY::SparseArray3D  # Y edge size
	EZ::SparseArray3D  # Z edge size
  NC::SparseArray3D  # CellNumbering
	NFX::SparseArray3D # X FaceNumbering
	NFY::SparseArray3D # Y FaceNumbering
	NFZ::SparseArray3D # Z FaceNumbering
	NEX::SparseArray3D # X EdgeNumbering
	NEY::SparseArray3D # Y EdgeNumbering
	NEZ::SparseArray3D # Z EdgeNumbering
  NN::SparseArray3D  # NodalNumbering
	dim::Int           # Mesh dimension
end # type OcTreeMeshFV


function getOcTreeMeshFV(S,h;x0=zeros(3))

    # get number of cells
    NC = getCellNumbering(S)
    nc = nnz(NC)

		# get number of faces
    FX,FY,FZ, NFX, NFY, NFZ = getFaceSizeNumbering(S)
		nf = [nnz(FX), nnz(FY), nnz(FZ)]

		# get number of edges
    EX,EY,EZ, NEX, NEY, NEZ = getEdgeSizeNumbering(S)
		ne = [nnz(EX), nnz(EY), nnz(EZ)]

		# get number of nodes
    NN = getNodalNumbering(S)
		nn = nnz(NN)

		empt = spzeros(0,0)
		empt3 = sparse3([size(S,1),size(S,2),size(S,3)])

		return OcTreeMeshFV(S, h, x0, S.sz,
                        nc,nf,ne,nn,
                        empt,empt,empt,       # no Div, Grad, Curl
                        Dict{Int64,MassMatrix}(), # no Pf
                        Dict{Int64,MassMatrix}(), # no Pe
                        Dict{Int64,MassMatrix}(), # no Pn
                        empt,empt,empt,   # no Af,Ae,An
                        empt,empt,empt,empt, # no V,L,Ne,Qe,
                        Int64[],Int64[],Int64[],  # no active edges, active faces, active nodes
                        empt,empt,empt,empt, #no Nn,Qn,Nf,Qf
                        FX,FY,FZ,
                        EX,EY,EZ,
                        NC,
                        NFX, NFY, NFZ,
                        NEX, NEY, NEZ,
                        NN,
                        3)
end  # function getOcTreeMeshFV

import Base.clear!
function clear!(M::OcTreeMeshFV)
  M.S    = sparse3([0,0,0])
  M.h    = zeros(3)
  M.x0   = zeros(3)
  M.n    = [0,0,0]
  M.nc   = 0
  M.nf   = [0,0,0]
  M.ne   = [0,0,0]
  M.nn   = 0
  M.Div  = spzeros(0,0)
  M.Grad = spzeros(0,0)
  M.Curl = spzeros(0,0)
  M.Pf   = Dict{Int64,MassMatrix}()
  M.Pe   = Dict{Int64,MassMatrix}()
  M.Pn   = Dict{Int64,MassMatrix}()
  M.Af   = spzeros(0,0)
  M.Ae   = spzeros(0,0)
  M.An   = spzeros(0,0)
  M.V    = spzeros(0,0)
  M.L    = spzeros(0,0)
  M.Ne   = spzeros(0,0)
  M.Qe   = spzeros(0,0)
  M.activeEdges = Array{Int64}(0)
  M.activeFaces = Array{Int64}(0)
  M.activeNodes = Array{Int64}(0)
  M.Nn   = spzeros(0,0)
  M.Qn   = spzeros(0,0)
  M.Nf   = spzeros(0,0)
  M.Qf   = spzeros(0,0)
  M.FX   = sparse3([0,0,0])
  M.FY   = sparse3([0,0,0])
  M.FZ   = sparse3([0,0,0])
  M.EX   = sparse3([0,0,0])
  M.EY   = sparse3([0,0,0])
  M.EZ   = sparse3([0,0,0])
  M.NC   = sparse3([0,0,0])
  M.NFX  = sparse3([0,0,0])
  M.NFY  = sparse3([0,0,0])
  M.NFZ  = sparse3([0,0,0])
  M.NEX  = sparse3([0,0,0])
  M.NEY  = sparse3([0,0,0])
  M.NEZ  = sparse3([0,0,0])
  M.NN   = sparse3([0,0,0])
  return
end  # function clear

import Base.==
function ==(M1::OcTreeMeshFV,M2::OcTreeMeshFV)
	isEqual = trues(21)

	# check mandatory fields
	isEqual[1] =  M1.S==M2.S
	isEqual[2] =  M1.h==M2.h
	isEqual[3] =  (M1.x0    == M2.x0)
	isEqual[5] =  (M1.n     == M2.n)
	isEqual[6] =  (M1.nc    == M2.nc)
	isEqual[7] =  (M1.nf    == M2.nf)
	isEqual[8] =  (M1.ne    == M2.ne)

	# check fields that might be empty
	if !(isempty(M1.Div)) && !(isempty(M2.Div))
		isEqual[9] = (M1.Div == M2.Div)
	end
	if !(isempty(M1.Grad)) && !(isempty(M2.Grad))
		isEqual[10] = (M1.Grad == M2.Grad)
	end
	if !(isempty(M1.Curl)) && !(isempty(M2.Curl))
		isEqual[11] = (M1.Curl == M2.Curl)
	end
	if !(isempty(M1.Af)) && !(isempty(M2.Af))
		isEqual[12] = (M1.Af == M2.Af)
	end
	if !(isempty(M1.Ae)) && !(isempty(M2.Ae))
		isEqual[13] = (M1.Ae == M2.Ae)
	end
	if !(isempty(M1.V)) && !(isempty(M2.V))
		isEqual[15] = (M1.V == M2.V)
	end
	if !(isempty(M1.Ne)) && !(isempty(M2.Ne))
		isEqual[20] = (M1.Ne == M2.Ne)
	end
	if !(isempty(M1.Nn)) && !(isempty(M2.Nn))
		isEqual[21] = (M1.Nn == M2.Nn)
	end

	return all(isEqual)
end  # function ==
