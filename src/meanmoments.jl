# compute the dynamic change in mean over time
function dμ_dt(network::DiffEqBase.AbstractReactionNetwork, momentstruct::MomentStruct; algebra=sympy())
    species_dict = speciesmap(network)
    dμ = repeat([zero(algebra.Sym)], momentstruct.num_species)
    for reaction in network.reactions
        for reactant in union(getfield.(reaction.products, :reactant), getfield.(reaction.substrates, :reactant))
            propensity = mass_rate_MA(reaction; algebra=algebra)
            dμ_i = algebra.expand(taylor_expand(propensity, momentstruct; algebra=algebra))
            dμ_i *= DiffEqBiological.get_stoch_diff(reaction, reactant)
            dμ[species_dict[reactant]] += dμ_i
        end
    end
    dμ
end

