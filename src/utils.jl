# calculate the stoichiometry coefficients
function stoich_coeff(species::Vector{T}, reaction::DiffEqBiological.ReactionStruct, e_order::Vector{Int}) where T <: AlgebraSet
    coeff = 1
    for ind in eachindex(e_order)
        coeff *= DiffEqBiological.get_stoch_diff(reaction, Symbol(species[ind]))^e_order[ind]
    end
    coeff
end

# calculate the sign coefficients
function sign_coeff(a_ords::Vector{Int}, b_ords::Vector{Int})
    *((-1).^(a_ords .- b_ords)...)
end

# calculate the exponentiated mean terms
function α_coeff(raw_means::Vector{T}, n_orders::Vector{Int}, k_orders::Vector{Int}) where T <: AlgebraSet
    *(raw_means.^(n_orders .- k_orders)...)
end

# calculate the bionomial coefficients
function binomial_coeff(a_ords::Vector{Int}, b_ords::Vector{Int})
    *(binomial.(a_ords, b_ords)...)
end

# calculate the 1 / n! coefficients
function inv_factorial_calc(orders::Vector{Int})
    1 / *(factorial.(orders)...)
end

# to rewrite the mass rate calculation of DiffEqBiological, so that it can be differentiated later on.
function mass_rate_MA(reaction::DiffEqBiological.ReactionStruct; algebra::Algebra=sympy())
    #println(SymPy.Sym(string(reaction.rate_org)))
    rate = algebra.Sym(string(reaction.rate_org))
    if reaction.is_pure_mass_action
        for sub in reaction.substrates
             for i in 1:sub.stoichiometry
                rate *= (algebra.Sym(sub.reactant) + 1 - i) / i
            end
        end
    end
    rate
end

# convert raw moments to central moments
function raw_to_cen(momentstruct::MomentStruct, n_order::Vector{Int})
    eq = zero(eltype(momentstruct.species))
    for ord in Iterators.product([0:n_max for n_max in n_order]...)
        __n =collect(ord)
        eq += binomial_coeff(n_order, __n) * α_coeff(momentstruct.species, n_order, __n) * momentstruct.cen_moments[__n]
    end
    eq
end

# convert central moments to raw moments (currently unused)
function cen_to_raw(momentstruct::MomentStruct, n_order::Vector{Int})
	eq = zero(eltype(momentstruct.species))
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

# compute taylor expansion
function taylor_expand(expr::AlgebraSet, momentstruct::MomentStruct; algebra::Algebra=sympy())
    eq = zero(algebra)
    for orders in keys(momentstruct.cen_moments)
        diff_args = reduce(append!, repeat([momentstruct.species[i]], orders[i]) for i in eachindex(orders))
        factorial_coeff = inv_factorial_calc(orders)
        eq += foldl(algebra.diff, (expr, diff_args...))*momentstruct.cen_moments[orders]*factorial_coeff
    end
    eq
end

# generate ode expressions and function from symbolic arrays
function get_ode(network::DiffEqBase.AbstractReactionNetwork, moments::MomentStruct, dμ::Vector{T}, dM::Vector{T}) where T <: AlgebraSet
    filtered_cen_moments = filter(p->sum(p.first)>1, moments.cen_moments)
    central_dict = OrderedDict{Symbol, Int}(Symbol(pair[2])=>i+moments.num_species for (i,pair) in enumerate(filtered_cen_moments))
    ode_expr = Vector{Expr}(undef, moments.num_species + length(central_dict))
    for i = 1:moments.num_species
        ex = dμ[i]
        ode_expr[i] = :(internal_var___du[$i] = $(Meta.parse(string(ex))))
    end
    for i = (moments.num_species+1):(moments.num_species+length(central_dict))
        ex = dM[i-moments.num_species]
        ode_expr[i] = :(internal_var___du[$i] = $(Meta.parse(string(ex))))
    end
    combined_dict = merge(speciesmap(network), central_dict)
    ode_func = DiffEqBiological.make_func(ode_expr, combined_dict, paramsmap(network))
    ode_expr, ODEFunction(eval(ode_func))
end
