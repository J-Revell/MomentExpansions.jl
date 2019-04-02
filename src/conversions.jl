# convert raw moments to central moments
function raw_to_cen(momentstruct::MomentStruct, n_order::Vector{Int})
    eq = SymPy.Sym(0)
    for ord in Iterators.product([0:n_max for n_max in n_order]...)
        __n =collect(ord)
        eq += binomial_coeff(n_order, __n) * α_coeff(momentstruct.species, n_order, __n) * momentstruct.cen_moments[__n]
    end
    eq
end

# convert central moments to raw moments (currently unused)
function cen_to_raw(momentstruct::MomentStruct, n_order::Vector{Int})
	eq = SymPy.Sym(0)
	for ord in Iterators.product([0:n_max for n_max in n_order]...)
		__n =collect(ord)
		if sum(__n) == 1
			eq += binomial_coeff(n_order, __n) * sign_coeff(n_order, __n) * α_coeff(momentstruct.species, n_order, __n) * momentstruct.species[findall(isone,__n)[1]]
		else
			eq += binomial_coeff(n_order, __n) * sign_coeff(n_order, __n) * α_coeff(momentstruct.species, n_order, __n) * momentstruct.raw_moments[__n]
		end
	end
	eq
end
