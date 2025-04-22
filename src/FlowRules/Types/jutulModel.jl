export jutulModel, CartesianMesh

struct jutulModel{D, T}
    n::NTuple{D, Int64}
    d::NTuple{D, T}
    ϕ::Vector{T}
    K::Union{Matrix{T}, T}
    h::T    # depth of the top grid
end

# Constructor from scalar ϕ
function jutulModel(n::NTuple{D, Int64}, d::NTuple{D, T}, ϕ::T, K::Union{Matrix{T}, T}; h::T=T(0), pad::Bool=true) where {D, T}
    ϕ_full = ϕ * ones(T, n)
    if pad
        ϕ_full[1,:,:] .= 1e8
        ϕ_full[end,:,:] .= 1e8
        ϕ_full[:,:,1] .= ϕ * (h/d[end] + T(1))
        ϕ_full[:,:,end] .= 1e8
    end
    ϕ_vec = vec(ϕ_full)
    return jutulModel(n, d, ϕ_vec, K, h)
end

# Constructor from vector ϕ with given depth
jutulModel(n::NTuple{D, Int64}, d::NTuple{D, T}, ϕ::Vector{T}, K::Union{Matrix{T}, T}; h::T=T(0)) where {D, T} =
    jutulModel{D, T}(n, d, ϕ, K, h)

# Display override
display(M::jutulModel{D, T}) where {D, T} =
    println("$(D)D jutulModel with size $(M.n) and grid spacing $(M.d) at depth $(M.h) m")

# Create Cartesian mesh
CartesianMesh(M::jutulModel{D, T}) where {D, T} = CartesianMesh(M.n, M.d .* M.n)

# Equality
==(A::jutulModel{D, T}, B::jutulModel{D, T}) where {D, T} = (A.n == B.n && A.d == B.d && A.ϕ == B.ϕ && A.K == B.K && A.h == B.h)
