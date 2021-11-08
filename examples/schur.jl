### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ 76f9e6bd-6c96-4441-8f80-dec874ea5a11
begin
	import Pkg
	Pkg.activate(".")
	Pkg.develop(; path="..")
	Pkg.add(["AbstractAlgebra", "DynamicPolynomials", "Transducers", "Dictionaries"])
end

# ╔═╡ ee3c9bca-3fa0-11ec-3fc4-55bf8b43afb4
using SchurPolynomials, AbstractAlgebra, DynamicPolynomials, Transducers, Dictionaries

# ╔═╡ 3049ac04-180c-45a1-9571-58afdb69dbbc
using Transducers: Map

# ╔═╡ 4a3fb69b-13c7-4fdf-b12f-3bd1e2abdaba
using AbstractAlgebra: Partition

# ╔═╡ 480c3e86-ef27-437d-872d-612792f259f9
λ = Partition([2, 1, 1])

# ╔═╡ fd2f0cdb-0674-4546-bae8-15dcacc5a18c
collect(semi_standard_young_tableaux(4, λ))

# ╔═╡ 53eee919-2c7b-4c0f-bd86-9e715d2b0129
n = 5

# ╔═╡ 97356e16-b186-4d7a-8dc4-0f062dc0c91f
S = Dict(λ.part => schur(n, λ) for λ in Generic.partitions(n))

# ╔═╡ 5ce04038-e4a1-4ab2-bb67-b2e379747e56
powerfun(k::Int, x) = sum(xᵢ -> xᵢ^k, x)

# ╔═╡ b626a35f-6732-41db-a3d9-c53a906bc95a
powerfun(ν::Generic.Partition, x) = prod(k -> powerfun(k, x), ν)

# ╔═╡ c6903995-493e-40d6-bf40-4b1aef56f224
function z(ν::Generic.Partition)
	row_counts = foldl(right, GroupBy(identity, (x, _) -> x + 1, 0), ν)
	prod(pairs(row_counts)) do (k, rₖ)
		factorial(rₖ) * k^rₖ
	end
end

# ╔═╡ d407ded7-6613-4068-b55c-c1869e9178fe
@polyvar x[1:n]

# ╔═╡ d4b75503-c43a-4c58-880b-628f5e9e0972
S′ = Dict(
	λ.part => convert(Polynomial{true, Int}, sum(Generic.partitions(n)) do ν
		character(λ, ν) // z(ν) * powerfun(ν, x)
	end)
	for λ in Generic.partitions(n)
)

# ╔═╡ ceebd22f-447c-4310-8a52-b22ec98ce876
S

# ╔═╡ 8a08ba51-f000-4771-854a-301133c7c405
let x = rand(n)
	(Ref(x) .|> Dictionary(S)) .- (Ref(x) .|> Dictionary(S′))
end

# ╔═╡ Cell order:
# ╠═76f9e6bd-6c96-4441-8f80-dec874ea5a11
# ╠═ee3c9bca-3fa0-11ec-3fc4-55bf8b43afb4
# ╠═3049ac04-180c-45a1-9571-58afdb69dbbc
# ╠═4a3fb69b-13c7-4fdf-b12f-3bd1e2abdaba
# ╠═480c3e86-ef27-437d-872d-612792f259f9
# ╠═fd2f0cdb-0674-4546-bae8-15dcacc5a18c
# ╠═53eee919-2c7b-4c0f-bd86-9e715d2b0129
# ╠═97356e16-b186-4d7a-8dc4-0f062dc0c91f
# ╠═5ce04038-e4a1-4ab2-bb67-b2e379747e56
# ╠═b626a35f-6732-41db-a3d9-c53a906bc95a
# ╠═c6903995-493e-40d6-bf40-4b1aef56f224
# ╠═d407ded7-6613-4068-b55c-c1869e9178fe
# ╠═d4b75503-c43a-4c58-880b-628f5e9e0972
# ╠═ceebd22f-447c-4310-8a52-b22ec98ce876
# ╠═8a08ba51-f000-4771-854a-301133c7c405