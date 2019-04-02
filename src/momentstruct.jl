"""
    MomentStruct

Structure storing:
-order::Int
-num_species::Int
-species::Vector{SymPy.Sym}
-raw_moments::OrderedDict{Vector{Int}, SymPy.Sym}
-cen_moments::OrderedDict{Vector{Int}, SymPy.Sym}
"""
struct MomentStruct
    order::Int
    num_species::Int
    species::Vector{SymPy.Sym}
    raw_moments::OrderedDict{Vector{Int}, SymPy.Sym}
    cen_moments::OrderedDict{Vector{Int}, SymPy.Sym}
end


"""
    MomentStruct

Given an AbstractReactionNetwork, and expansion order, generates a MomentStruct.
"""
function MomentStruct(network::DiffEqBase.AbstractReactionNetwork, order::Int)
    species_dict = speciesmap(network)
    species = SymPy.Sym.(collect(keys(species_dict)))
    num_species = length(species_dict)
    num_reactions = length(network.reactions)

    # zero order moments
    raw_moments = OrderedDict(zeros(Int, num_species) => SymPy.Sym(1))
    cen_moments = OrderedDict(zeros(Int, num_species) => SymPy.Sym(1))

    # map out first order moments
    first_order_keys = Vector{Int}[]
    for i in 1:num_species
        mkey = zeros(Int, num_species)
        mkey[i] = 1
        push!(first_order_keys, mkey)
    end
    first_order_syms = [SymPy.Sym(key) for key in keys(species_dict)]

    # first order central moments == 0
    cen_first_order_syms = [SymPy.Sym(0) for orders in first_order_keys]

    merge!(raw_moments, OrderedDict(first_order_keys .=> first_order_syms))
    merge!(cen_moments, OrderedDict(first_order_keys .=> cen_first_order_syms))

    # map out high order raw moments
    high_order_keys = collect.([ords for ords in Base.product(repeat([0:order], num_species)...) if 1<sum(ords)<=order])
    sort!(high_order_keys, by=sum)
    raw_high_order_syms = [SymPy.Sym("x___"*join([string(ords) for ords in orders])) for orders in high_order_keys]
    merge!(raw_moments, OrderedDict(high_order_keys .=> raw_high_order_syms))

    # map out high order central moments
    cen_high_order_syms = [SymPy.Sym("M___"*join([string(ords) for ords in orders])) for orders in high_order_keys]
    merge!(cen_moments, OrderedDict(high_order_keys .=> cen_high_order_syms))

    MomentStruct(order, num_species, species, raw_moments, cen_moments)
end
