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
P(k::Int, x) = sum(i -> x[i]^k, 1:length(x))

# ╔═╡ a7542146-8fb4-470e-a63d-9e39b9286d5d
function P(ρ::Generic.Partition, x)
	row_counts = foldl(right, GroupBy(identity, (x, _) -> x + 1, 0), ρ)
	prod(pairs(row_counts)) do (k, rₖ)
		1//(factorial(rₖ) * k^rₖ) * P(k, x)^rₖ
	end
end

# ╔═╡ d407ded7-6613-4068-b55c-c1869e9178fe
@polyvar x[1:n]

# ╔═╡ fa0daebb-15c6-4b37-ab75-bf16be9bb6b0
S′ = Dict(
	λ.part => convert(Polynomial{true, Int}, sum(Generic.partitions(n)) do ρ
		character(λ, ρ) * P(ρ, x)
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
# ╠═a7542146-8fb4-470e-a63d-9e39b9286d5d
# ╠═d407ded7-6613-4068-b55c-c1869e9178fe
# ╠═fa0daebb-15c6-4b37-ab75-bf16be9bb6b0
# ╠═ceebd22f-447c-4310-8a52-b22ec98ce876
# ╠═8a08ba51-f000-4771-854a-301133c7c405
