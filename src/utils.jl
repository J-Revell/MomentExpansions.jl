# calculate the stoichiometry coefficients
function stoich_coeff(species::Vector{SymPy.Sym}, reaction::DiffEqBiological.ReactionStruct, e_order::Vector{Int})
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
function α_coeff(raw_means::Vector{SymPy.Sym}, n_orders::Vector{Int}, k_orders::Vector{Int})
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
function mass_rate_MA(reaction::DiffEqBiological.ReactionStruct)
    rate = SymPy.Sym(reaction.rate_org)
    if reaction.is_pure_mass_action
        for sub in reaction.substrates
             for i in 1:sub.stoichiometry
                rate *= (SymPy.Sym(sub.reactant) + 1 - i) / i
            end
        end
    end
    rate
end

# function to generate ode expressions and to create the simulation function
function get_ode(network::DiffEqBase.AbstractReactionNetwork, moments::MomentStruct, dμ::Vector{SymPy.Sym}, dM::Vector{SymPy.Sym})
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