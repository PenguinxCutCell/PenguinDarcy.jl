abstract type AbstractMobilityTensor end

struct ScalarMobility{M} <: AbstractMobilityTensor
    data::M
end

struct DiagonalMobility{M} <: AbstractMobilityTensor
    data::M
end

struct FullTensorMobility{M} <: AbstractMobilityTensor
    data::M
end

_wrap_mobility(λ::AbstractMobilityTensor, N::Int) = λ
_wrap_mobility(λ::Number, N::Int) = ScalarMobility(λ)
_wrap_mobility(λ::Function, N::Int) = ScalarMobility(λ)

function _wrap_mobility(λ::NTuple{N,<:Any}, ::Int) where {N}
    return DiagonalMobility(λ)
end

function _wrap_mobility(λ::AbstractVector, N::Int)
    length(λ) == N || throw(DimensionMismatch("diagonal mobility vector length must be $N"))
    return DiagonalMobility(λ)
end

function _wrap_mobility(λ::AbstractMatrix, N::Int)
    size(λ, 1) == N && size(λ, 2) == N || throw(DimensionMismatch("mobility tensor matrix must be $N×$N"))
    return FullTensorMobility(λ)
end

function _mobility_matrix(λ::Function, x::SVector{N,T}, t::T) where {N,T}
    if applicable(λ, x..., t)
        M = λ(x..., t)
    elseif applicable(λ, x...)
        M = λ(x...)
    elseif applicable(λ, t)
        M = λ(t)
    elseif applicable(λ)
        M = λ()
    else
        throw(ArgumentError("full tensor mobility callback must accept (x..., t), (x...), (t), or ()"))
    end
    size(M, 1) == N && size(M, 2) == N || throw(DimensionMismatch("full tensor mobility callback must return $N×$N matrix"))
    return M
end

function _mobility_entry(mob::ScalarMobility, i::Int, j::Int, x::SVector{N,T}, t::T, row_lin::Int) where {N,T}
    i == j || return zero(T)
    return convert(T, eval_coeff(mob.data, x, t, row_lin))
end

function _mobility_entry(mob::DiagonalMobility, i::Int, j::Int, x::SVector{N,T}, t::T, row_lin::Int) where {N,T}
    i == j || return zero(T)
    di = mob.data[i]
    return di isa Number ? convert(T, di) : _eval_fun_or_const(di, x, t)
end

function _mobility_entry(mob::FullTensorMobility, i::Int, j::Int, x::SVector{N,T}, t::T, row_lin::Int) where {N,T}
    M = if mob.data isa AbstractMatrix
        mob.data
    elseif mob.data isa Function
        _mobility_matrix(mob.data, x, t)
    else
        throw(ArgumentError("unsupported full tensor mobility payload $(typeof(mob.data))"))
    end
    return convert(T, M[i, j])
end

function _phase_mobility(model::DarcyModelMono{N,T}, phase::Int=1) where {N,T}
    return _wrap_mobility(model.λ, N)
end

function _phase_mobility(model::DarcyModelDiph{N,T}, phase::Int) where {N,T}
    λ = phase == 1 ? model.λ1 : model.λ2
    return _wrap_mobility(λ, N)
end

function _face_coordinate(cap::AssembledCapacity{N,T}, I, d::Int, LI) where {N,T}
    lin = LI[I]
    cω = cap.C_ω
    xlow = convert(T, cap.xyz[d][1])
    if I[d] == 1
        return SVector{N,T}(ntuple(k -> (k == d ? xlow : cω[lin][k]), N))
    end
    Iminus = CartesianIndex(ntuple(k -> (k == d ? I[k] - 1 : I[k]), N))
    linm = LI[Iminus]
    return SVector{N,T}(ntuple(k -> (cω[lin][k] + cω[linm][k]) / T(2), N))
end

function face_tensor_values(
    model,
    ops::DiffusionOps{N,T},
    cap::AssembledCapacity{N,T},
    i::Int,
    j::Int;
    t::T,
    phase::Int=1,
) where {N,T}
    1 <= i <= N || throw(ArgumentError("i must satisfy 1 <= i <= $N"))
    1 <= j <= N || throw(ArgumentError("j must satisfy 1 <= j <= $N"))

    nt = cap.ntotal
    vals = Vector{T}(undef, nt)
    LI = LinearIndices(cap.nnodes)
    mob = _phase_mobility(model, phase)

    @inbounds for I in CartesianIndices(cap.nnodes)
        lin = LI[I]
        if any(k -> I[k] == cap.nnodes[k], 1:N)
            vals[lin] = zero(T)
            continue
        end
        xface = _face_coordinate(cap, I, i, LI)
        vals[lin] = _mobility_entry(mob, i, j, xface, t, lin)
    end
    return vals
end

function _add_block!(A::SparseMatrixCSC{T,Int}, rows::UnitRange{Int}, cols::UnitRange{Int}, B::SparseMatrixCSC{T,Int}) where {T}
    @inbounds for j in 1:size(B, 2)
        for p in nzrange(B, j)
            i = B.rowval[p]
            A[rows[i], cols[j]] += B.nzval[p]
        end
    end
    return A
end

function tensor_flux_operator(
    model,
    ops::DiffusionOps{N,T},
    cap::AssembledCapacity{N,T};
    t::T,
    phase::Int=1,
) where {N,T}
    nt = cap.ntotal
    nfaces = N * nt
    Tf = spzeros(T, nfaces, nfaces)
    pflags = periodic_flags(model.bc_border, N)

    for i in 1:N
        ri = ((i - 1) * nt + 1):(i * nt)
        for j in 1:N
            cj = ((j - 1) * nt + 1):(j * nt)
            Dij = spdiagm(0 => face_tensor_values(model, ops, cap, i, j; t=t, phase=phase))
            Iji = cross_face_interp(cap, j, i; periodic_flags=pflags)
            _add_block!(Tf, ri, cj, Dij * Iji)
        end
    end

    return Tf
end
