### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 11111111-1111-1111-1111-111111111111
begin
	using Oscar: ZZ, QQ, matrix_space
	using Latexify
end

# ╔═╡ 2a2a2a2a-2a2a-2a2a-2a2a-2a2a2a2a2a2a
begin
	# Global Mathematical Domain Bindings
	Z = ZZ
	Q = QQ
	Num_Nodes = Z(3)
	
	Space_Z = matrix_space(Z, 3, 3)
	Space_Q = matrix_space(Q, 3, 3)
	
	# Adjacency Matrix representing direct walks
	Adj_Alg = Space_Z([
		0 1 1;
		1 0 1;
		1 0 0
	])
	
	# Transition Matrix for Markov Chain (probabilities as rationals)
	M_Alg = Space_Q([
		0 Q(1,2) Q(1,2);
		Q(1,2) 0 Q(1,2);
		1 0 0
	])
	
	# Initial State Vector (Vector Space over Q)
	V_Q = matrix_space(Q, 3, 1)
	s_0 = V_Q([Q(1,3), Q(1,3), Q(1,3)])
end

# ╔═╡ 22222222-2222-2222-2222-222222222222
md"""
# Lecture Seventeen: Ways to Walk—Matrices and Markov Chains

## Scope
- Focuses on **combinatorics of paths** and **Markov chains**.
- Explores how multiplying adjacency matrices reveals the number of walks in a graph.
- Introduces Markov chains for analyzing systems that change randomly over time.

---

## I. Matrices as walks
A. The number of $k$-step walks from vertex $i$ to vertex $j$ is the $(i,j)$ entry of $A^k$.
B. Problem: Find the number of 3-step walks between all nodes.
"""

# ╔═╡ 3a3a3a3a-3a3a-3a3a-3a3a-3a3a3a3a3a3a
begin
	# SECTION I - PARADIGMS
	# Calculate 3-step walks
	
	# 1. Pure FP
	walks_fp = let
		n = Int(Num_Nodes)
		a_raw = [Int(Adj_Alg[i,j]) for i in 1:n, j in 1:n]
		
		# Matrix multiplication recursively
		mat_mul = (X, Y) -> [sum(X[i,k]*Y[k,j] for k in 1:n) for i in 1:n, j in 1:n]
		
		a2 = mat_mul(a_raw, a_raw)
		mat_mul(a2, a_raw)
	end
	
	# 2. Abstract Algebra
	walks_aa = let
		A3 = Adj_Alg^3
		[Int(A3[i,j]) for i in 1:3, j in 1:3]
	end
	
	# 3. Geometric Focus
	walks_grp = let
		A3 = Adj_Alg * Adj_Alg * Adj_Alg
		[Int(A3[i,j]) for i in 1:3, j in 1:3]
	end
end

# ╔═╡ 3b3b3b3b-3b3b-3b3b-3b3b-3b3b3b3b3b3b
begin
	Markdown.parse("""
	### Aggregate Verification (Section I):
	* **Pure FP**: $(walks_fp)
	* **Abstract Algebra**: $(walks_aa)
	* **Geometric**: $(walks_grp)
	
	All paradigms agree: **$(walks_fp == walks_aa == walks_grp)**
	""")
end

# ╔═╡ 44444444-4444-4444-4444-444444444444
md"""
---

## II. Markov Chains
A. A Markov chain describes a sequence of possible events where the probability of each event depends only on the state attained in the previous event.
B. Problem: Find the state vector after 2 transitions.
"""

# ╔═╡ 4a4a4a4a-4a4a-4a4a-4a4a-4a4a4a4a4a4a
begin
	# SECTION II - PARADIGMS
	# Calculate the state after 2 steps (M^2 * s_0 for left-multiplication, but here rows are states? 
	# Wait, usually it is s_0^T * M^2. Let's compute s_2).
	
	# 1. Pure FP
	state_2_fp = let
		n = Int(Num_Nodes)
		m_raw = [M_Alg[i,j] for i in 1:n, j in 1:n]
		s_raw = [s_0[i,1] for i in 1:n]
		
		# s_1 = M^T * s_0
		s1 = [sum(m_raw[j,i] * s_raw[j] for j in 1:n) for i in 1:n]
		s2 = [sum(m_raw[j,i] * s1[j] for j in 1:n) for i in 1:n]
		s2
	end
	
	# 2. Abstract Algebra
	state_2_aa = let
		M2 = M_Alg^2
		s2_mat = transpose(M2) * s_0
		[s2_mat[i,1] for i in 1:3]
	end
	
	# 3. Geometric Focus
	state_2_grp = let
		s1_mat = transpose(M_Alg) * s_0
		s2_mat = transpose(M_Alg) * s1_mat
		[s2_mat[i,1] for i in 1:3]
	end
end

# ╔═╡ 4b4b4b4b-4b4b-4b4b-4b4b-4b4b4b4b4b4b
begin
	Markdown.parse("""
	### Aggregate Verification (Section II):
	* **Pure FP State 2**: $(state_2_fp)
	* **Abstract Algebra State 2**: $(state_2_aa)
	* **Geometric State 2**: $(state_2_grp)
	
	All paradigms agree: **$(state_2_fp == state_2_aa == state_2_grp)**
	""")
end

# ╔═╡ 55555555-5555-5555-5555-555555555555
md"""
---

## III. The Stable Matrix
A. A regular Markov chain eventually reaches a stable state (equilibrium) independent of the initial state.
B. Problem: Find the limit matrix $M^\infty$ representing the steady state. For this simple 3-node graph, we will compute $M^4$ as a proxy for approaching stability.
"""

# ╔═╡ 5a5a5a5a-5a5a-5a5a-5a5a-5a5a5a5a5a5a
begin
	# SECTION III - PARADIGMS
	# Calculate M^4
	
	# 1. Pure FP
	m4_fp = let
		n = Int(Num_Nodes)
		m_raw = [M_Alg[i,j] for i in 1:n, j in 1:n]
		mat_mul = (X, Y) -> [sum(X[i,k]*Y[k,j] for k in 1:n) for i in 1:n, j in 1:n]
		m2 = mat_mul(m_raw, m_raw)
		mat_mul(m2, m2)
	end
	
	# 2. Abstract Algebra
	m4_aa = let
		M4 = M_Alg^4
		[M4[i,j] for i in 1:3, j in 1:3]
	end
	
	# 3. Geometric Focus
	m4_grp = let
		M2 = M_Alg * M_Alg
		M4 = M2 * M2
		[M4[i,j] for i in 1:3, j in 1:3]
	end
end

# ╔═╡ 5b5b5b5b-5b5b-5b5b-5b5b-5b5b5b5b5b5b
begin
	Markdown.parse("""
	### Aggregate Verification (Section III):
	* **Pure FP**: $(m4_fp)
	* **Abstract Algebra**: $(m4_aa)
	* **Geometric**: $(m4_grp)
	
	All paradigms agree: **$(m4_fp == m4_aa == m4_grp)**
	""")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
Oscar = "f1435218-dba5-11e9-1e4d-f1a5fab5fc13"

[compat]
Latexify = "~0.16.10"
Oscar = "~1.7.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.6"
manifest_format = "2.0"
project_hash = "1ab06f7055e32512a4ef690bbead085f66024b1b"
"""

# ╔═╡ Cell order:
# ╠═11111111-1111-1111-1111-111111111111
# ╠═2a2a2a2a-2a2a-2a2a-2a2a-2a2a2a2a2a2a
# ╟─22222222-2222-2222-2222-222222222222
# ╠═3a3a3a3a-3a3a-3a3a-3a3a-3a3a3a3a3a3a
# ╠═3b3b3b3b-3b3b-3b3b-3b3b-3b3b3b3b3b3b
# ╟─44444444-4444-4444-4444-444444444444
# ╠═4a4a4a4a-4a4a-4a4a-4a4a-4a4a4a4a4a4a
# ╠═4b4b4b4b-4b4b-4b4b-4b4b-4b4b4b4b4b4b
# ╟─55555555-5555-5555-5555-555555555555
# ╠═5a5a5a5a-5a5a-5a5a-5a5a-5a5a5a5a5a5a
# ╠═5b5b5b5b-5b5b-5b5b-5b5b-5b5b5b5b5b5b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
