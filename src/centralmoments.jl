# generate higher order functions of the propensity to be expanded
function gen_fs(network::DiffEqBase.AbstractReactionNetwork, momentstruct::MomentStruct, k_orders::Vector{Int}, e_orders::Vector{Int})
	propensities = mass_rate_MA.(network.reactions)
	propensities .* α_coeff(momentstruct.species, k_orders, e_orders)
end

# expand higher order functions of the propensity
function expand_fs(momentstruct::MomentStruct, fs::Vector{SymPy.Sym})
	e_of_fs = [taylor_expand(f, momentstruct) for f in fs]
end

# dα term
function dα_dt(momentstruct::MomentStruct, dμ::Vector{SymPy.Sym}, n_orders::Vector{Int}, k_orders::Vector{Int})
	dα = SymPy.Sym(0)
	for i in eachindex(n_orders)
		dα += (n_orders[i] - k_orders[i]) * (1/momentstruct.species[i]) * α_coeff(momentstruct.species, n_orders, k_orders) * dμ[i]
	end
	dα
end

# dβ term
function dβ_dt(network::DiffEqBase.AbstractReactionNetwork, momentstruct::MomentStruct, k_order::Vector{Int})
	eq = repeat([SymPy.Sym(0)], length(network.reactions))
	for e_order in Iterators.product([0:k_max for k_max in k_order]...)
		if sum(e_order) == 0
			continue
		end
		fs = gen_fs(network, momentstruct, k_order, collect(e_order))
		e_of_fs = expand_fs(momentstruct, fs)
		coeffs = [binomial_coeff(k_order, collect(e_order)) * stoich_coeff(momentstruct.species, reaction, collect(e_order)) for reaction in network.reactions]
		eq += coeffs .* e_of_fs
	end
	sum(eq)
end

# central moments
function dcen_dt(network::DiffEqBase.AbstractReactionNetwork, momentstruct::MomentStruct, dμ::Vector{SymPy.Sym})
	dcen = repeat([SymPy.Sym(0)], length(momentstruct.cen_moments) - momentstruct.num_species -1)
	counter = 0
	for n_order in keys(momentstruct.cen_moments)
		if sum(n_order) <= 1
			continue
		else
			counter += 1
			for k_order in Iterators.product([0:n_max for n_max in n_order]...)
				__k = collect(k_order)
				β = momentstruct.raw_moments[__k]
				α = α_coeff(momentstruct.species, n_order, __k)
				dβ = dβ_dt(network, momentstruct, __k)
				dα = dα_dt(momentstruct, dμ, n_order, __k)
				nck = binomial_coeff(n_order, __k)
				__sign = sign_coeff(n_order, __k)
				__β = raw_to_cen(momentstruct, __k)
				dcen[counter] += nck*__sign*(α*dβ+__β*dα)
			end
		end
	end
	SymPy.simplify.(dcen)
end
