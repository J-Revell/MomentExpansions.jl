module MomentExpansions
    using DataStructures
    import SymPy
    import SymEngine
    import Base.zero, Base.one, Base.show
    using Reexport
    @reexport using DiffEqBase
    @reexport using DiffEqBiological

    include("types.jl")
    include("utils.jl")
    include("meanmoments.jl")
    include("centralmoments.jl")

    """
		approx(network::AbstractReactionNetwork, order::Int; algebra::Algebra)

	Given an AbstractReactionNetwork, and an expansion order, returns the standard moment expansion.
    To switch between using SymPy or SymEngine as the backend, use the kwargs `algebra=sympy()` or `algebra=symengine()`.
	"""
    function approx(network::DiffEqBase.AbstractReactionNetwork; order::Int=2, algebra::Algebra=symengine())
        momentstruct = gen_moment_syms(network, order, algebra=algebra)
        dmean = dμ_dt(network, momentstruct, algebra=algebra)
        dcen = dcen_dt(network, momentstruct, dmean, algebra=algebra)
        ode_expr, ode_func = get_ode(network, momentstruct, dmean, dcen)
        MomentExpansion(network, momentstruct, dmean, dcen, ode_expr, ode_func)
    end

    function DiffEqBase.ODEProblem(momentexpansion::MomentExpansion, u0::Union{AbstractArray, Number}, args...; kwargs...)
        DiffEqBase.ODEProblem(momentexpansion.ode, u0, args...; kwargs...)
    end

    export approx, sympy, symengine
end
