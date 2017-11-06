export importUBCOcTreeMesh, importUBCOcTreeModel

"""
    mesh = importUBCOcTreeMesh(meshfile)

    Reads an OcTree mesh in UBC format from disk.

    Input:

        meshfile::AbstractString - File to read

    Output:

        mesh::OcTreeMesh - The mesh

"""
function importUBCOcTreeMesh(meshfile::AbstractString;Tn::Type{N}=Int64,Tn2::Type{N2}=Int64) where N <: Integer where N2 <: Integer

    # open file (throws error if file doesn't exist)
    f    = open(meshfile,"r")

    # number of cells of underlying tensor mesh along dimension 1, 2, 3
    line = split(readline(f))
    m1 = parse(Tn2,line[1])
    m2 = parse(Tn2,line[2])
    m3 = parse(Tn2,line[3])

    # top corner coordinates
    line = split(readline(f))
    x1 = parse(Float64,line[1])
    x2 = parse(Float64,line[2])
    x3 = parse(Float64,line[3])

    # cell size
    line = split(readline(f))
    h1 = parse(Float64,line[1])
    h2 = parse(Float64,line[2])
    h3 = parse(Float64,line[3])

    # number of OcTree cells
    line = split(readline(f))
    n = parse(Tn,line[1])

    # read rest of file at ones
    lines = readlines(f)

    # close file
    close(f)

    # check correct number of lines read
    if n != length(lines)
        error("Invalid number of (i,j,k,bsz) lines in file $meshfile.")
    end

    # allocate space
    i1  = zeros(Tn2, n)
    i2  = zeros(Tn2, n)
    i3  = zeros(Tn2, n)
    bsz = zeros(Tn , n)

    # convert string array to numbers
    for i = 1:n
        line   = split(lines[i])

        i1[i]  = parse(Tn2,line[1])
        i2[i]  = parse(Tn2,line[2])
        i3[i]  = parse(Tn2,line[3])
        bsz[i] = parse(Tn ,line[4])
    end

    # Roman's code starts the OcTree at the top corner. Change to bottom
    # corner.
    i3 = m3 + 2 .- i3 - bsz
    x3 = x3 - m3 * h3

    #S   = sortrows([i3 i2 i1 bsz])
    S = sub2ind( (m1,m2,m3), i1,i2,i3 )
    p = sortpermFast(S)[1]

    i1  = i1[p] #S[:,3]
    i2  = i2[p] #S[:,2]
    i3  = i3[p] #S[:,1]
    bsz = bsz[p] #S[:,4]

    # create mesh object
    S = sparse3(i1,i2,i3,bsz,[m1,m2,m3])
    mesh = getOcTreeMeshFV(S,[h1,h2,h3];x0=[x1,x2,x3])
    return mesh
end

"""
    model = importUBCOcTreeModel(modelfile, mesh, DataType=Float64)

    Reads an OcTree model in UBC format from disk.

    Input:

        modelfile::AbstractString - File to read
        mesh::OcTreeMeshFV        - The corresponding mesh
        T::DataType               - Data type of model (Float64, Int64, Bool, ...)

    Output:

        model::Array{Float64,1} - The model
        model::Array{Float64,2}

"""
function importUBCOcTreeModel(modelfile::AbstractString, mesh::OcTreeMesh, T::DataType=Float64)

    # open file (throws error if file doesn't exist)
    f = open(modelfile,"r")

    # read everything
    s = readlines(f)

    # close
    close(f)

    # check if we have the correct number of cell values
    n = mesh.nc
    if length(s) != n
        error("Incorrect number of cell values")
    end

    # Roman's code starts the OcTree at the top corner. Here, we start with the bottom corner. Therefore, we need to permute the cells values.
    m1,m2,m3 = mesh.n
    i1,i2,i3,bsz = find3(mesh.S)
    i3 = m3 + 2 .- i3 - bsz
    S = sub2ind( (m1,m2,m3), i1,i2,i3 )
    p = sortpermFast(S)[1]

    # check for multiple values per cell
  	d = split(s[1])
  	ncol = length(d)

  	# convert to numbers
    if ncol == 1
        model = Array{T}(n)
    else
        model = Array{T}(n,ncol)
    end
    for i = 1:n
        idx = p[i]
        d = split(s[i])
        for j = 1:ncol
            model[idx,j] = parse(T,d[j])
        end
    end

    return model
end
