### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 76f9e6bd-6c96-4441-8f80-dec874ea5a11
begin
	import Pkg
	Pkg.activate(".")
	#Pkg.add(; url="git@github.com:simeonschaub/YAPP.jl", rev="main")
	Pkg.develop(; path="/home/simeon/.julia/dev/YAPP")
	Pkg.develop(; path="..")
	Pkg.add(["AbstractAlgebra", "Transducers", "Dictionaries", "BenchmarkTools", "PlutoUI", "Plots"]) 
end

# ╔═╡ ee3c9bca-3fa0-11ec-3fc4-55bf8b43afb4
begin
	using SchurPolynomials, AbstractAlgebra, YAPP, Transducers, Dictionaries, LinearAlgebra
	using YAPP: Polynomial, Monomial
	using SchurPolynomials: schur
	using Transducers: Map
	using AbstractAlgebra: Partition
end

# ╔═╡ 2bf97bef-38c3-4140-b44d-b8bb9745673f
using Printf

# ╔═╡ 9e830e0c-2cab-4247-8487-ece4b6f82854
using BenchmarkTools

# ╔═╡ 555e1339-0d12-402b-957c-f0d8d70cbc42
using Plots

# ╔═╡ 480c3e86-ef27-437d-872d-612792f259f9
λ = Partition([2, 1, 1])

# ╔═╡ 5987e81e-b1d7-49cf-8bf1-83a43eece501
md"""
$$\begin{align}
s_\lambda(x_1, \dots, x_n) = \sum_{T \in \mathrm{SSYT}(\lambda)} x^T
= \sum_{T \in \mathrm{SSYT}(\lambda)} x_1^{T_1} \cdots x_n^{T_n}
\end{align}$$
"""

# ╔═╡ 53eee919-2c7b-4c0f-bd86-9e715d2b0129
n = 5

# ╔═╡ 97356e16-b186-4d7a-8dc4-0f062dc0c91f
S = Dict(YoungTableau(λ.part) => schur(n, λ) for λ in Generic.partitions(n))

# ╔═╡ d65b2b91-95cd-4bc5-a50c-60dfb994647a
md"""
$$\begin{align}
s_\lambda(x_1, \dots, x_n) = \sum_{\nu \vdash n} \frac{\chi_\nu^\lambda}{z_\nu} \cdot p_\nu(x_1, \dots, x_n)
\end{align}$$
"""

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
x = Monomial{Vector{Int}}.(eachcol(I(n)))#@polyvar x[1:n]

# ╔═╡ b7c2ff8f-ff71-499b-84fe-d0cc77ee25d3
Base.copy(m::Monomial) = Monomial(copy(m.exponents))

# ╔═╡ d4b75503-c43a-4c58-880b-628f5e9e0972
S′ = Dict(
	YoungTableau(λ.part) => copy(sum(Generic.partitions(n)) do ν
		character(λ, ν) // z(ν) * powerfun(ν, x)
	end, Int)
	for λ in Generic.partitions(n)
)

# ╔═╡ ceebd22f-447c-4310-8a52-b22ec98ce876
S

# ╔═╡ 8a08ba51-f000-4771-854a-301133c7c405
let x = rand(n)
	(Ref(x) .|> Dictionary(S)) .- (Ref(x) .|> Dictionary(S′))
end

# ╔═╡ 9ee2519b-786a-428e-a517-7e599c11f3e4
let m=11, n=6
	x = Monomial{Vector{Int}}.(eachcol(I(n)))
	x_ = randn(n, m)
	p, s = (Matrix{Float64}(undef, length(Generic.partitions(n)), m) for _ in 1:2)
	for (i, λ) in enumerate(Generic.partitions(n))
		p[i, :] .= (eachcol(x_) .|> powerfun(λ, x[1:n])) ./ z(λ)
		s[i, :] .= (eachcol(x_) .|> schur(n, λ))
	end
	s / p .|> x -> round(Int, x)
	Text.(Printf.format.([Printf.format"%.2f"], s / p))
end

# ╔═╡ b3f710d2-5216-47a8-8427-3b0a5d0f2dd7
@benchmark schur(10, Partition([4, 3, 2, 1]))

# ╔═╡ f80625c9-0bc7-410c-a1aa-cbd80829520a
function next_partition!(p)
	i = rand(0:length(p))
	if i == 0
		push!(p, 1)
	elseif get(p, i-1, typemax(Int)) > p[i]
		p[i] += 1
	else
		next_partition!(p)
	end
	p
end

# ╔═╡ 3eef3ce3-c303-41ff-b608-6b5653fc1dc4
p = accumulate(1:12; init=Int[]) do p, _
	next_partition!(copy(p))
end

# ╔═╡ 766ca80a-5918-4d1b-bd6c-d0dba151208c
times = [@belapsed(schur($n, $(Partition(p)))) for (n, p) in enumerate(p)]

# ╔═╡ 5cbe2be9-1308-4bb0-8955-8666ad489110
scatter(times; xscale=:log10, yscale=:log10, legend=:none)

# ╔═╡ Cell order:
# ╠═76f9e6bd-6c96-4441-8f80-dec874ea5a11
# ╠═ee3c9bca-3fa0-11ec-3fc4-55bf8b43afb4
# ╠═2bf97bef-38c3-4140-b44d-b8bb9745673f
# ╠═480c3e86-ef27-437d-872d-612792f259f9
# ╟─5987e81e-b1d7-49cf-8bf1-83a43eece501
# ╠═53eee919-2c7b-4c0f-bd86-9e715d2b0129
# ╠═97356e16-b186-4d7a-8dc4-0f062dc0c91f
# ╟─d65b2b91-95cd-4bc5-a50c-60dfb994647a
# ╠═5ce04038-e4a1-4ab2-bb67-b2e379747e56
# ╠═b626a35f-6732-41db-a3d9-c53a906bc95a
# ╠═c6903995-493e-40d6-bf40-4b1aef56f224
# ╠═d407ded7-6613-4068-b55c-c1869e9178fe
# ╠═b7c2ff8f-ff71-499b-84fe-d0cc77ee25d3
# ╠═d4b75503-c43a-4c58-880b-628f5e9e0972
# ╠═ceebd22f-447c-4310-8a52-b22ec98ce876
# ╠═8a08ba51-f000-4771-854a-301133c7c405
# ╠═9ee2519b-786a-428e-a517-7e599c11f3e4
# ╠═9e830e0c-2cab-4247-8487-ece4b6f82854
# ╠═b3f710d2-5216-47a8-8427-3b0a5d0f2dd7
# ╠═f80625c9-0bc7-410c-a1aa-cbd80829520a
# ╠═3eef3ce3-c303-41ff-b608-6b5653fc1dc4
# ╠═766ca80a-5918-4d1b-bd6c-d0dba151208c
# ╠═555e1339-0d12-402b-957c-f0d8d70cbc42
# ╠═5cbe2be9-1308-4bb0-8955-8666ad489110
