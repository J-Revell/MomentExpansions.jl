const AlgebraSet = Union{SymEngine.Basic, SymPy.Sym}

# To facilitate the switching between SymPy and SymEngine function calls.
struct Algebra
    Sym::DataType
    diff::Function
    expand::Function
end
Base.zero(A::Algebra) = zero(A.Sym)
Base.one(A::Algebra) = one(A.Sym)

# Declaring the algebra functions
symengine() = Algebra(SymEngine.Basic, SymEngine.diff, SymEngine.expand)
sympy() = Algebra(SymPy.Sym, SymPy.diff, SymPy.expand)

# Type for storing dicts regarding the symbolic expressions for the moments.
struct MomentStruct{T <: AlgebraSet}
    order::Int
    num_species::Int
    species::Vector{T}
    raw_moments::OrderedDict{Vector{Int}, T}
    cen_moments::OrderedDict{Vector{Int}, T}
end
MomentStruct(order, num_species, species, raw_moments, cen_moments; algebra=sympy()) = MomentStruct{algebra.Sym}(order, num_species, species, raw_moments, cen_moments)

# function to generate the moment dict
function gen_moment_syms(network::DiffEqBase.AbstractReactionNetwork, order::Int; algebra=sympy())
    species_dict = speciesmap(network)
    species = algebra.Sym.(collect(keys(species_dict)))
    num_species = length(species_dict)
    num_reactions = length(network.reactions)

    # zero order moments
    raw_moments = OrderedDict(zeros(Int, num_species) => one(algebra.Sym))
    cen_moments = OrderedDict(zeros(Int, num_species) => one(algebra.Sym))

    # map out first order moments
    first_order_keys = Vector{Int}[]
    for i in 1:num_species
        mkey = zeros(Int, num_species)
        mkey[i] = 1
        push!(first_order_keys, mkey)
    end
    first_order_syms = [algebra.Sym(key) for key in keys(species_dict)]

    # first order central moments == 0
    cen_first_order_syms = [zero(algebra) for orders in first_order_keys]

    merge!(raw_moments, OrderedDict(first_order_keys .=> first_order_syms))
    merge!(cen_moments, OrderedDict(first_order_keys .=> cen_first_order_syms))

    # map out high order raw moments
    high_order_keys = collect.([ords for ords in Base.product(repeat([0:order], num_species)...) if 1<sum(ords)<=order])
    sort!(high_order_keys, by=sum)
    raw_high_order_syms = [algebra.Sym("x___"*join([string(ords) for ords in orders])) for orders in high_order_keys]
    merge!(raw_moments, OrderedDict(high_order_keys .=> raw_high_order_syms))

    # map out high order central moments
    cen_high_order_syms = [algebra.Sym("M___"*join([string(ords) for ords in orders])) for orders in high_order_keys]
    merge!(cen_moments, OrderedDict(high_order_keys .=> cen_high_order_syms))

    MomentStruct(order, num_species, species, raw_moments, cen_moments; algebra=algebra)
end

# function to store the moment expansion
struct MomentExpansion{T <: AlgebraSet}
    network::DiffEqBase.AbstractReactionNetwork
    moments::MomentStruct{T}
    dmean::Vector{T}
    dcen::Vector{T}
    ode_expr::Vector{Expr}
    ode::Function
end
MomentExpansion(network, moments, dmean, dcen, ode_expr, ode; algebra=sympy()) = MomentExpansion{algebra.Sym}(network, moments, dmean, dcen, ode_expr, ode)

Base.show(io::IO,M::MomentExpansions.MomentExpansion)=print("MomentExpansion, order: ", M.moments.order)
