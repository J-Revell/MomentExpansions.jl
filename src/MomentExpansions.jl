__precompile__()

module MomentExpansions

	using DataStructures
	import SymPy

	using Reexport
	@reexport using DiffEqBase
	@reexport using DiffEqBiological

	include("momentstruct.jl")
	include("utils.jl")
	include("taylor.jl")
	include("meanmoments.jl")
	include("conversions.jl")
	include("centralmoments.jl")

	struct MomentExpansion
		network::DiffEqBase.AbstractReactionNetwork
		moments::MomentStruct
		dμ::Vector{SymPy.Sym}
		dM::Vector{SymPy.Sym}
		ode_expr::Vector{Expr}
		ode::Function
	end

	"""
		MomentExpansion(network::AbstractReactionNetwork, order::Int)


	Given an AbstractReactionNetwork, and an expansion order, returns the standard moment expansion.
	"""
	function MomentExpansion(network::DiffEqBase.AbstractReactionNetwork, order::Int)
		moments = MomentStruct(network, order)
		dμ = dμ_dt(network, moments)
		dM = dcen_dt(network, moments, dμ)
		ode_expr, ode_func = get_ode(network, moments, dμ, dM)
		MomentExpansion(network, moments, dμ, dM, ode_expr, ode_func)
	end

	const MomentProblem = DiffEqBase.ODEProblem

	### create ODEProblem from MomentExpansion ###
	function DiffEqBase.ODEProblem(momentexpansion::MomentExpansion, u0::Union{AbstractArray, Number}, args...; kwargs...)
		DiffEqBase.ODEProblem(momentexpansion.ode, u0, args...; kwargs...)
	end

	### avoid confusion between diffeqbase and sympy ###
	#solve(args...; kwargs...) = DiffEqBase.solve(args...; kwargs...)

	export MomentExpansion, MomentStruct, ODEProblem#, solve

end # module
