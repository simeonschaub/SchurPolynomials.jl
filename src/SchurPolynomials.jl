module SchurPolynomials

using AbstractAlgebra, Transducers, DynamicPolynomials
using AbstractAlgebra: Partition
using Transducers: next, complete, Foldable
import MutableArithmetics as MA

const RowsType = Vector{SubArray{Int, 1, Vector{Int}, Tuple{UnitRange{Int}}, true}}

"""
Iterates over all semi-standard Young tableaux of shape `λ`.
"""
struct SemiStandardYoungTableaux <: Foldable
    T::Generic.YoungTableau{Int}
    rows::RowsType
    n::Int
    function SemiStandardYoungTableaux(λ, n)
        fill = Vector{Int}(undef, sum(λ))
        rows = map(λ, Iterators.accumulate(+, λ)) do len, stop
            view(fill, stop-len+1:stop)
        end
        return new(YoungTableau(λ, fill), rows, n)
    end
end

"""
Iterates over all configurations of row `i` that are allowed by the previous row.
Iterator returns `nothing` but mutates `row`.
"""
struct RowConfigs <: Foldable
    rows::RowsType
    i::Int
    n::Int
end
function Transducers.__foldl__(rf, acc, (; rows, i, n)::RowConfigs)
    row = @inbounds rows[i]
    acc = _row!(i, 1, row, rows, n, rf, acc)
    return complete(rf, acc)
end
Base.eltype(::RowConfigs) = Nothing

# row i, column j
function _row!(i, j, row, rows, n, rf, acc)
    j > length(row) && return next(rf, acc, nothing)
    above = i == 1 ? 0 : rows[i-1][j]
    left = get(row, j-1, 1)
    @inbounds for row[j] in max(above + 1, left):n
        acc = _row!(i, j + 1, row, rows, n, rf, acc)
        acc isa Reduced && return acc
    end
    acc
end

function _worker!(i, rf, acc, ssyt)
    (; T, rows, n) = ssyt
    i > length(rows) && return next(rf, acc, T)
    foldl(RowConfigs(rows, i, n); init=acc) do acc, _
        _worker!(i + 1, rf, acc, ssyt)
    end
end

function Transducers.__foldl__(rf, acc, ssyt::SemiStandardYoungTableaux)
    acc = _worker!(1, rf, acc, ssyt)
    return complete(rf, acc)
end

function Base.length((; T, n)::SemiStandardYoungTableaux)
    return Int(prod(
        (n-i+j) // hooklength(T, i, j)
        for (i, j_max) in enumerate(T.part)
        for j in 1:j_max
    ))
end
Base.eltype(::SemiStandardYoungTableaux) = Generic.YoungTableau

function schur(n, λ)
    SSYT = SemiStandardYoungTableaux(λ, n)
    k = length(SSYT)
    @polyvar x[1:n]
    monomials = [Vector{Int}(undef, n) for _ in 1:k]
    foreach(Enumerate(), SSYT) do (j, T)
        map!(monomials[j], 1:n) do i
            count(==(i), T.fill)
        end
        nothing
    end
    return DynamicPolynomials.polynomialclean(x, fill(1, k), monomials)
end

end
