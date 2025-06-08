

"""
    vec2ws(v::Vector{Matrix{Float64}})

Convert a vector of matrices into a Workspace object.

Each matrix in the input vector `v` is assigned to a key in the output Workspace.
The keys are symbols representing the 1-based index of the matrix in the vector
(e.g., the first matrix `v[1]` is stored under the key `:1`, `v[2]` under `:2`, and so on).

# Arguments
- `v::Vector{Matrix{Float64}}`: A vector of `Matrix{Float64}`.

# Returns
- `ws::Workspace`: A Workspace object where each key `Symbol(i)` corresponds to the matrix `v[i]`.
"""
function vec2ws(v::Vector{Matrix{Float64}})
    ws = Workspace()
    for (i,mat) ∈ enumerate(v)
        ws[Symbol(i)] = mat
    end
    return ws
end

"""
    vec2mat(v::Vector{Matrix{Float64}})

Convert a vector of square matrices into a single matrix containing their diagonals.

This function takes a vector of matrices, extracts the diagonal from each matrix,
and stores these diagonals as rows in a new matrix. It's assumed all matrices
in the input vector are square and have the same dimensions.

# Arguments
- `v::Vector{Matrix{Float64}}`: A vector of square matrices.

# Returns
- `mat::Matrix{Float64}`: A matrix where the i-th row is the diagonal of the i-th matrix from the input vector `v`. The number of rows is `length(v)` and the number of columns is `size(v[1],1)`.
"""
function vec2mat(v::Vector{Matrix{Float64}})
    mat = zeros(length(v), size(v[1],1))
    for (i,_mat) ∈ enumerate(v)
        mat[i,:] = diag(_mat)
    end
    return mat
end

"""
    save_chain_to_daec(path, ws)

Save estimation chain results from a Workspace to a DAEC database file.

This function preprocesses certain fields in the input Workspace `ws` before saving.
Specifically, `ws.𝚲draws` (a vector of matrices) is converted to a Workspace using `vec2ws`,
and `ws.𝛙draws` and `ws.Rdraws` (vectors of matrices) are converted to matrices of their
diagonals using `vec2mat`. Other fields are deep-copied, except for a predefined list
of fields to exclude (`:𝚲draws`, `:𝛙draws`, `:Rdraws`, `:S`, `:Λ_constraints`, `:SR`, `:fixed_parameters`).
The processed Workspace is then written to the specified DAEC file path.

# Arguments
- `path::String`: The file path (including filename) where the DAEC database will be saved.
- `ws::Workspace`: The Workspace object containing the estimation chain results.
                  Expected to have fields like `𝚲draws`, `𝛙draws`, `Rdraws`, and other
                  relevant results.

# Returns
- Nothing. The function writes to a file.

# Dependencies
- Requires the `DE` module (presumably DataFramesECO.jl or similar) for `DE.writedb`.
- Uses `vec2ws` and `vec2mat` for transformations.
"""
function save_chain_to_daec(path,ws)
    # convert nested matrices 
    _ws = Workspace()
    _ws.𝚲 = vec2ws(ws.𝚲draws)
    _ws.𝛙 = vec2mat(ws.𝛙draws)
    _ws.R = vec2mat(ws.Rdraws)

    for key ∈ keys(ws)
        if key ∈ (:𝚲draws, :𝛙draws, :Rdraws, :S, :Λ_constraints, :SR, :fixed_parameters)
            continue
        end
        _ws[key] = deepcopy(ws[key])        
    end
    
    DE.writedb(path,_ws)
end

"""
    store_results(res::Workspace, m::MacroModelling.ℳ, name = "", burnin=50000, compress=true)

Store estimation results Workspace to a DAEC file with a standardized naming convention.

This function serves as a wrapper around `save_chain_to_daec`. It constructs a filename
based on the model name (or a custom `name`) and the current timestamp.
The default name is `m.model_name`, but if `m.model_name` is "SW", it defaults to "gelfer".

# Arguments
- `res::Workspace`: The Workspace object containing the estimation results to be saved.
- `m::MacroModelling.ℳ`: The model object, used to derive a default filename.
- `name = ""`: An optional custom name for the chain. If empty, a name is derived from `m.model_name`.
- `burnin=50000`: This argument is present but not explicitly used by `store_results` itself or `save_chain_to_daec`. It might be a remnant or intended for other uses.
- `compress=true`: This argument is present but not explicitly used by `store_results` itself or `save_chain_to_daec`. It might be a remnant or intended for other uses.

# Returns
- Nothing. The function calls `save_chain_to_daec` which writes to a file.

# File Naming
The output file will be saved in the "data/results/" directory with a name like
`[derived_name]_[timestamp].daec`, e.g., "gelfer_2023-10-27T1530.daec".
"""
function store_results(res::Workspace, m::MacroModelling.ℳ, name = "")
    if name == ""
        name = m.model_name
        if m.model_name == "SW"
            name = "gelfer"
        end
    end
    save_chain_to_daec("""data/results/$(name)_$(Dates.format(now(),"yyyy-mm-ddTHHMM")).daec""",res)
end

"""
    load_results(path, full=false)

Load estimation results from a DAEC database file into a Workspace.

This function can load either the full database or a specific subset of commonly
used fields (`𝛉`, `𝚲`, `𝛙`, `R`). The `𝚲` field, if loaded as a Workspace
(from `vec2ws` during save), is converted back into a vector of matrices.

# Arguments
- `path::String`: The file path of the DAEC database to load.
- `full=false`: If `true`, loads the entire database content into the Workspace.
                If `false` (default), loads only the `𝛉`, `𝚲`, `𝛙`, and `R` fields.

# Returns
- `res::Workspace`: A Workspace object containing the loaded results. If `𝚲` was loaded,
                  it's converted from a Workspace of matrices (keyed by index symbols)
                  back to a `Vector{Matrix{Float64}}`.

# Dependencies
- Requires the `DE` module (presumably DataFramesECO.jl or similar) for `DE.readdb`.
"""
function load_results(path, full=false)
    res = Workspace()
    if full
        res = DE.readdb(path)
    else
        res.𝛉 = DE.readdb(path, "𝛉")
        res.𝚲 = DE.readdb(path, "𝚲")
        res.𝛙 = DE.readdb(path, "𝛙")
        res.R = DE.readdb(path, "R")
    end
    𝚲 = repeat([zeros(size(res.𝚲[Symbol(1)]))], length(keys(res.𝚲)))
    for (i,sym) ∈ enumerate(keys(res.𝚲))
        𝚲[i] = res.𝚲[sym]
    end
    res.𝚲 = 𝚲
    
    return res 
end