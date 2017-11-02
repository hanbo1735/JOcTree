#module JOcTree


export OcTreeMeshFV, getOcTreeMeshFV


"""
    struct MassMatrix

Internal storage for mass matrix integration
 - `M = getNodalMassMatrix(mesh, sigma)`
 - `M = getEdgeMassMatrix(mesh, sigma)`
 - `M = getFaceMassMatrix(mesh, sigma)`

Fields
 - `n::Integer`: size of mass matrix `M` (number of nodes, edges, faces in `mesh`)
 - `A::SparseMatrixCSC{Float64,Integer}`: map coefficient vector `sigma` to
      nonzero entries in mass matrix `M` (numerical integration)
 - `rowval::Array{Integer,1}`: row indices of nonzero entries in `M`
 - `colptr::Array{Integer,1}`: CSC format column pointers in `M`
 - `colval::Array{Integer,1}`: column indices of nonzero entries in `M`
"""
struct MassMatrix{T<:Real,N<:Integer}
    n::Int
    A::SparseMatrixCSC{T,N}
    rowval::Array{N,1}
    colptr::Array{N,1}
    colval::Array{N,1}
end


"""
    MassMatrix()

Default constructor for `MassMatrix`
"""
function MassMatrix(T::Type{t},N::Type{n}) where t <: Real where n <: Integer
    F = Array{T}(0)
    I = Array{T}(0)
    S = SparseMatrixCSC(0, 0, N[1], I, F)
    return MassMatrix(0, S, I, I, I)
end


mutable struct OcTreeMeshFV{T<:Real,N<:Integer,N2<:Integer} <: OcTreeMesh
	S::SparseArray3D{N,N2}    # i,j,k, bsz
	h::Vector{T}  # (3) underlying cell size
	x0::Vector{T} # coordinates of corner of mesh
	n::Vector{N2} # underlying mesh
	nc::N          # number of cells
	nf::Vector{N}  # (3) number of faces
	ne::Vector{N}  # (3) number of edges
    nn::N          # number of nodes
	Div::SparseMatrixCSC{T,N}
	Grad::SparseMatrixCSC{T,N}
	Curl::SparseMatrixCSC{T,N}
	Pf::Dict{Int64,MassMatrix{T,N}} # face mass matrix integration storage
	Pe::Dict{Int64,MassMatrix{T,N}} # edge mass matrix integration storage
	Pn::Dict{Int64,MassMatrix{T,N}} # nodal mass matrix integration storage
	Af::SparseMatrixCSC{T,N} # Face to cell-centre matrix
	Ae::SparseMatrixCSC{T,N} # Edge to cell-centre matrix
	An::SparseMatrixCSC{T,N} # Node to cell-centre matrix
	V::SparseMatrixCSC{T,N} # cell volume
	L::SparseMatrixCSC{T,N} # edge lengths
	Ne::SparseMatrixCSC{T,N} # Edge nullspace matrix
	Qe::SparseMatrixCSC{T,N} # Edge projection matrix
	activeEdges::Vector{N}   # lookup table for new edge enumeration
	activeFaces::Vector{N}   # lookup table for new face enumeration
	activeNodes::Vector{N}   # lookup table for new node enumeration
	Nn::SparseMatrixCSC{T,N} # Node nullspace matrix
	Qn::SparseMatrixCSC{T,N} # Node projection matrix
	Nf::SparseMatrixCSC{T,N} # Face nullspace matrix
	Qf::SparseMatrixCSC{T,N} # Face projection matrix
	FX::SparseArray3D{N,N2}  # X face size
	FY::SparseArray3D{N,N2}  # Y face size
	FZ::SparseArray3D{N,N2}  # Z face size
	EX::SparseArray3D{N,N2}  # X edge size
	EY::SparseArray3D{N,N2}  # Y edge size
	EZ::SparseArray3D{N,N2}  # Z edge size
    NC::SparseArray3D{N,N2}  # CellNumbering
	NFX::SparseArray3D{N,N2} # X FaceNumbering
	NFY::SparseArray3D{N,N2} # Y FaceNumbering
	NFZ::SparseArray3D{N,N2} # Z FaceNumbering
	NEX::SparseArray3D{N,N2} # X EdgeNumbering
	NEY::SparseArray3D{N,N2} # Y EdgeNumbering
	NEZ::SparseArray3D{N,N2} # Z EdgeNumbering
    NN::SparseArray3D{N,N2}  # NodalNumbering
	dim::N           # Mesh dimension
end # type OcTreeMeshFV


function getOcTreeMeshFV{T<:Real,N<:Integer,N2<:Integer}(S::SparseArray3D{N,N2},
                         h::Vector{T};x0::Vector{T}=zeros(T,3))

    # get number of cells
    NC = getCellNumbering(S)
    nc = N(nnz(NC))

	# get number of faces
    FX,FY,FZ, NFX, NFY, NFZ = getFaceSizeNumbering(S)
	nf = convert(Vector{N},[nnz(FX), nnz(FY), nnz(FZ)])

	# get number of edges
    EX,EY,EZ, NEX, NEY, NEZ = getEdgeSizeNumbering(S)
	ne = convert(Vector{N},[nnz(EX), nnz(EY), nnz(EZ)])

	# get number of nodes
    NN = getNodalNumbering(S)
	nn = N(nnz(NN))

	empt  = spzeros(T,N,0,0)
    sz    = [size(S,1),size(S,2),size(S,3)]
	empt3 = SparseArray3D(spzeros(N,N2,prod(sz)),sz)

		return OcTreeMeshFV(S, h, x0, S.sz,
                        nc,nf,ne,nn,
                        empt,empt,empt,       # no Div, Grad, Curl
                        Dict{Int64,MassMatrix{T,N}}(), # no Pf
                        Dict{Int64,MassMatrix{T,N}}(), # no Pe
                        Dict{Int64,MassMatrix{T,N}}(), # no Pn
                        empt,empt,empt,   # no Af,Ae,An
                        empt,empt,empt,empt, # no V,L,Ne,Qe,
                        N[],N[],N[],  # no active edges, active faces, active nodes
                        empt,empt,empt,empt, #no Nn,Qn,Nf,Qf
                        FX,FY,FZ,
                        EX,EY,EZ,
                        NC,
                        NFX, NFY, NFZ,
                        NEX, NEY, NEZ,
                        NN,
                        N(3))
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
