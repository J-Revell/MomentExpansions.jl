# compute the dynamic change in mean over time
function dμ_dt(network::DiffEqBase.AbstractReactionNetwork, momentstruct::MomentStruct)
    species_dict = speciesmap(network)
    dμ = repeat([SymPy.Sym(0)], momentstruct.num_species)
    for reaction in network.reactions
        for reactant in union(getfield.(reaction.products, :reactant), getfield.(reaction.substrates, :reactant))
            propensity = mass_rate_MA(reaction)
            dμ_i = SymPy.simplify(taylor_expand(propensity, momentstruct))
            dμ_i *= DiffEqBiological.get_stoch_diff(reaction, reactant)
            dμ[species_dict[reactant]] += dμ_i
        end
    end
    dμ
end

# compute the dynamic change in mean over time,
# but expand orders into separate elements of an array
function dμ_dt_vec(network::DiffEqBase.AbstractReactionNetwork, momentstruct::MomentStruct)
    species_dict = speciesmap(network)
    dμ = repeat([repeat([SymPy.Sym(0)], momentstruct.order+1)], momentstruct.num_species)
    for reaction in network.reactions
        for reactant in union(getfield.(reaction.products, :reactant), getfield.(reaction.substrates, :reactant))
            propensity = mass_rate_MA(reaction)
            dμ_i = SymPy.simplify.(taylor_expand_vec(propensity, momentstruct))
            dμ_i .*= DiffEqBiological.get_stoch_diff(reaction, reactant)
            dμ[species_dict[reactant]] += dμ_i
        end
    end
    dμ
end
