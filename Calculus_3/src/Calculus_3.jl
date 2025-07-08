module Calculus_3
    using Oscar
   function gen_uni_polynomials!()
        XU, x_u = polynomial_ring(RealField(), :x_u)
        RR, X = polynomial_ring(RealField(), :X)
        YU, y_u = polynomial_ring(RealField(), :y_u)
        ZU, z_u = polynomial_ring(RealField(), :z_u)
        TU, t_u = polynomial_ring(RealField(), :t_u)
        DXU, dx_u = polynomial_ring(RealField(), :dx_u)
        DYU, dy_u = polynomial_ring(RealField(), :dy_u)
        DXU, dx_u = polynomial_ring(RealField(), :dx_u)
        DTU, dt_u = polynomial_ring(RealField(), :dt_u)
        DWU, dw_u = polynomial_ring(RealField(), :dw_u)
        ∆XU, ∆x_u = polynomial_ring(RealField(), :∆x_u)
        ∆YU, ∆y_u = polynomial_ring(RealField(), :∆y_u)
        ∂XU, ∂x_u = polynomial_ring(RealField(), :∂x_u)
        ∂YU, ∂y_u = polynomial_ring(RealField(), :∂y_u)
        ∂YU, ∂z_u = polynomial_ring(RealField(), :∂z_u)
        ∂WU, ∂w_u = polynomial_ring(RealField(), :∂w_u)
        ∂ƒU, ∂ƒ_u = polynomial_ring(RealField(), :∂ƒ_u)
    end


    gen_multi_polynomials!() = U, (a, A, b, c, C, d, dt, F, G, m, M, r, R, I, L, v, u, J, K, x, y, z, ƒ, t, ∫, θ, ∆x, ∆y, ∂x, ∂y, ∂z, ∂ƒ, ∞) = universal_polynomial_ring(QQ, [:a, :A, :b, :c, :C, :d, :dt, :F, :G, :m, :M, :r, :R, :I, :L, :v, :u, :J, :K, :x, :y, :z, :ƒ, :t, :∫, :θ, :∆x, :∆y, :∂x, :∂y, :∂z, :∂ƒ, :∞])
    
    to_ringElem(variable, field::Ring, value=1 ) = field(value)
    
    to_number(variable, value=1 ) = variable isa Number ? variable : value

    to_rational_coefficients(expr::Union{RealPolyRingElem, RealFieldElem}) = filter(!iszero, map(c -> rationalize(Float64(c)), coefficients(expr)))

    precision_real(real_expr::Union{RealPolyRingElem, RealFieldElem}, max::Int=4) = parse(Float64, !occursin("[", "$(real_expr)") ? SubString("$(real_expr)", 1:max) : SubString("$(real_expr)", 2:max))

    get_integrals_vars(expr::Union{RealPolyRingElem, RealFieldElem}, vars="x", oper_vars="*") =  map(s -> s[1:3], filter(c -> startswith(c[1:3], vars), split("$(expr)", oper_vars)))

    get_integrals(expr::Union{RealPolyRingElem, RealFieldElem}, oper_vars="*") =  map(filter(isnumeric), split("$(expr)", oper_vars))

    get_integrals_terms(to_rational_coefficients, get_integrals_vars) =  
        map(c -> string(c[1]) * c[2] , zip(to_rational_coefficients, get_integrals_vars))

    precision_integral(get_integrals, variable, precision::Int = 2) = length(get_integrals) >= 2 ? Rational(BigFloat(SubString(get_integrals[1], 1:precision)))*variable^ZZ(get_integrals[2]) : Rational(BigFloat(SubString(get_integrals[1], 1:precision)))*variable^ZZ(1)

    U(variable, parent = RealField()) = change_base_ring(parent, to_univariate(variable))

    Us(vector, parent = RealField()) = map(v -> U(v, parent), vector)

    terms_to_univariate(expr) = map(t -> U(t), terms(expr))

    to_univariate_terms(expr, var="", operator="*") = filter(c -> startswith(c[1:3], var), split("$(expr)", operator));

end # module Calculus_3