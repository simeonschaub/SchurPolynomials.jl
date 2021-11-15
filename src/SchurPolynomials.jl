module SchurPolynomials

export schur, semi_standard_young_tableaux

using AbstractAlgebra, Transducers, Dictionaries
using YAPP: Monomial, Polynomial
using BangBang
using AbstractAlgebra: Partition
using Transducers: next, complete, Foldable, Map

const RowsType = Vector{SubArray{Int, 1, Vector{Int}, Tuple{UnitRange{Int}}, true}}

"""
Iterates over all semi-standard Young tableaux of shape `λ`.
"""
struct SemiStandardYoungTableaux <: Foldable
    n::Int
    T::Generic.YoungTableau{Int}
    rows::RowsType
    function SemiStandardYoungTableaux(n, λ)
        fill = Vector{Int}(undef, sum(λ))
        rows = map(λ, Iterators.accumulate(+, λ)) do len, stop
            view(fill, stop-len+1:stop)
        end
        return new(n, YoungTableau(λ, fill), rows)
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
Base.eltype(::RowConfigs) = Nothing
function Transducers.__foldl__(rf::RF, acc, (; rows, i, n)::RowConfigs) where {RF}
    row = @inbounds rows[i]
    acc = _row!(i, 1, row, rows, n, rf, acc)
    return complete(rf, acc)
end
# row i, column j
function _row!(i, j, row, rows, n, rf::RF, acc) where {RF}
    j > length(row) && return next(rf, acc, nothing)
    above = i == 1 ? 0 : rows[i-1][j]
    left = get(row, j-1, 1)
    @inbounds for row[j] in max(above + 1, left):n
        acc = _row!(i, j + 1, row, rows, n, rf, acc)
        acc isa Reduced && return acc
    end
    acc
end

function _worker!(i, rf::RF, acc, ssyt) where {RF}
    (; n, T, rows) = ssyt
    i > length(rows) && return next(rf, acc, T)
    foldl(RowConfigs(rows, i, n); init=acc) do acc, _
        _worker!(i + 1, rf, acc, ssyt)
    end
end

function Transducers.__foldl__(rf::RF, acc, ssyt::SemiStandardYoungTableaux) where {RF}
    acc = _worker!(1, rf, acc, ssyt)
    return complete(rf, acc)
end

function Base.length((; n, T)::SemiStandardYoungTableaux)
    return Int(prod(
        (n-i+j) // hooklength(T, i, j)
        for (i, j_max) in enumerate(T.part)
        for j in 1:j_max
    ))
end
Base.eltype(::SemiStandardYoungTableaux) = Generic.YoungTableau

semi_standard_young_tableaux(n, λ) = SemiStandardYoungTableaux(n, λ) |> Map() do T
    YoungTableau(T.part, copy(T.fill))
end

function schur(n, λ)
    SSYT = SemiStandardYoungTableaux(n, λ)
    foldl(SSYT; init=Polynomial{Int}()) do p, T
        exponents = map(1:n) do i
            count(==(i), T.fill)
        end
        add!!(p, Monomial(exponents))
    end
end

end
