module SchurPolynomials

using AbstractAlgebra, FGenerators, Transducers, DynamicPolynomials
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

@fgenerator function row_configs!(rows, i, n)
	row = @inbounds rows[i]
	# row i, column j
	function _row!(j)
		j > length(row) && (@yield(row); return)
		above = i == 1 ? 0 : rows[i-1][j]
		left = get(row, j-1, 1)
		@inbounds for row[j] in max(above + 1, left):n
			_row!(j + 1)
		end
	end
	_row!(1)
end

@fgenerator(ssyt::SemiStandardYoungTableaux) do
    (; T, rows, n) = ssyt
	function worker!(i)
		i > length(rows) && (@yield(T); return)
		foreach(row_configs!(rows, i, n)) do _
			worker!(i + 1)
		end
	end
	worker!(1)
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
