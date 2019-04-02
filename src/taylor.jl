# perform the taylor expansion (w.r.t. central moments) on an expression
function taylor_expand(expression::SymPy.Sym, momentstruct::MomentStruct)
    eq = SymPy.Sym(0)
    for orders in keys(momentstruct.cen_moments)
        diff_args = reduce(append!, repeat([momentstruct.species[i]], orders[i]) for i in eachindex(orders))
        factorial_coeff = inv_factorial_calc(orders)
        eq += foldl(diff, (expression, diff_args...))*momentstruct.cen_moments[orders]*factorial_coeff
    end
    eq
end

# perform the taylor expansion (w.r.t. central moments), but separating different order terms into elements of an array
function taylor_expand_vec(expression::SymPy.Sym, momentstruct::MomentStruct)
    eq = repeat([SymPy.Sym(0)], momentstruct.order+1)
    for orders in keys(momentstruct.cen_moments)
        diff_args = reduce(append!, repeat([momentstruct.species[i]], orders[i]) for i in eachindex(orders))
        factorial_coeff = inv_factorial_calc(orders)
        eq[sum(orders)+1] += foldl(diff, (expression, diff_args...))*momentstruct.cen_moments[orders]*factorial_coeff
    end
    eq
end
