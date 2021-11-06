module SchurPolynomials

using AbstractAlgebra, Transducers, DynamicPolynomials
using AbstractAlgebra: Partition
using Transducers: next, @next, complete
import MutableArithmetics as MA

struct SemiStandardYoungTableaux <: Foldable
    T::Generic.YoungTableau{Int}
    rows::Vector{SubArray{Int64, 1, Vector{Int}, Tuple{UnitRange{Int}}, true}}
    n::Int
    function SemiStandardYoungTableaux(λ, n)
        fill = Vector{Int}(undef, sum(λ)) #[i for (i, j) in enumerate(λ) for _ in 1:j]
        rows = map(λ, Iterators.accumulate(+, λ)) do len, stop
            view(fill, stop-len+1:stop)
        end
        return new(YoungTableau(λ, fill), rows, n)
    end
end

function row_configs!(rows, i, n)
    Transducers.AdHocFoldable((; rows, i, n)) do rf, acc, (; rows, i, n)
        row = @inbounds rows[i]
        # row i, column j
        function _row!(j, acc)
            j > length(row) && return next(rf, acc, row)
            above = i == 1 ? 0 : rows[i-1][j]
            left = get(row, j-1, 1)
            @inbounds for row[j] in max(above + 1, left):n
                acc = _row!(j + 1, acc)
                acc isa Reduced && return acc
            end
            acc
        end
        acc = _row!(1, acc)
        return complete(rf, acc)
    end
end

function Transducers.__foldl__(rf, acc, (; T, rows, n)::SemiStandardYoungTableaux)
    function worker!(i, acc)
        i > length(rows) && return next(rf, acc, T)
        foldl(row_configs!(rows, i, n); init=acc) do acc, _
            worker!(i + 1, acc)
        end
    end
    acc = worker!(1, acc)
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
    end
    return DynamicPolynomials.polynomialclean(x, fill(1, k), monomials)
end

end
