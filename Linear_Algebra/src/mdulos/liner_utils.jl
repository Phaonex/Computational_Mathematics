module linearUtils

using Symbolics
using Latexify
using LinearAlgebra


export extract_coeffs, get_coefficients

extract_coeffs(eq, var, idx) = Symbolics.coeff(Symbolics.arguments(Symbolics.value(eq))[1], var[idx])

get_coefficients(eqs, var) = map(v -> extract_coeffs(eqs, var, v), range(1,length(var)))

dot_prods(v_1, v_2) = sum(v_1 .* v_2)

end; # module linearUtils
