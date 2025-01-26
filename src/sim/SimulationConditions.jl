# SimulationConditions

module SimulationConditions

#=
Packages
=#

using ..States: State
using ..DomainDescriptors: DomainDescriptor

#=
Public API
=#

export BoundTuple
export SimulationCondition, planeCouetteFlow, planePoiseuilleFlow
export init_flowfield!, apply_perturbation!

#=
Simulation Condition
=#

function check_option_validity(init_cond::Char, force_type::Char)

    # Admissible Set
    set_init_cond::Set{Char} = Set(['p', 'c', 'r'])
    set_force_type::Set{Char} = Set(['m', 'p'])

    # Error Toggle
    bool_init::Bool = false
    bool_force::Bool = false

    if !(init_cond in set_init_cond)
        bool_init = true
    end

    if !(force_type in set_force_type)
        bool_force = true
    end

    # Error Handling
    if (bool_init && bool_force)
        throw(
            ArgumentError(
                "init_cond argument ($(init_cond)) and force_type argument ($(force_type)) not in admissible option.",
            ),
        )
    end

    if (bool_init)
        throw(ArgumentError("init_cond argument ($(init_cond)) not in admissible option."))
    end

    if (bool_force)
        throw(
            ArgumentError("force_type argument ($(force_type)) not in admissible option."),
        )
    end

    return nothing

end

# BoundTuple NamedTuple
BoundTuple = @NamedTuple begin

    vel::Char # Velocity Component: 'u', 'v', 'w'
    wall::Char # Wall: upper wall 'u', lower wall 'l'

end

"""
SimulationCondition defines the simulation meta-variables (boundary conditions and forcing).
"""
struct SimulationCondition{T<:AbstractFloat}

    # Initial Condition
    init_cond::Char # 'p' for Poiseuille flow, 'c' for Couette flow, 'r' for restart
    init_kick::T

    # Pressure Gradient Characteristic
    force_constraint::Char # 'm' for mass flow rate, 'p' for pressure gradient
    force_magnitude::T

    # Boundary Condition
    bound_cond::Dict{BoundTuple,T} # Dict[(<velocity vector>,<wall>)] = velocity magnitude

    # Viscosity
    nu::T

    function SimulationCondition{T}(
        init_cond::Char,
        init_kick::T,
        force_type::Char,
        force_magnitude::T,
        bound_cond::Dict{BoundTuple,T},
        nu::T,
    ) where {T<:AbstractFloat}

        @info("--- Initialising SimulationCondition ---")

        # Validity Check
        check_option_validity(init_cond, force_type)

        @info("--- SimulationCondition initialised ---")
        @info(" ")

        return new{T}(init_cond, init_kick, force_type, force_magnitude, bound_cond, nu)

    end

end

#=
Canonical Configurations
=#

function planeCouetteFlow(init_cond::Char, init_kick::T) where {T<:AbstractFloat}

    # Initial Condition
    if (init_cond == 'p')
        throw(
            ArgumentError(
                "PlaneCouetteFlow flow configuration not compatible with init_cond argument.",
            ),
        )
    end

    init_kick::T = init_kick

    # Forcing Condition
    force_type::Char = 'p'
    force_magnitude::T = 0.0

    # Boundary Condition
    bound_cond::Dict{BoundTuple,T} = Dict()

    # u
    bound_cond[(vel = 'u', wall = 'l')] = -1.0
    bound_cond[(vel = 'u', wall = 'u')] = 1.0

    # v
    bound_cond[(vel = 'v', wall = 'l')] = 0.0
    bound_cond[(vel = 'v', wall = 'u')] = 0.0

    # w
    bound_cond[(vel = 'w', wall = 'l')] = 0.0
    bound_cond[(vel = 'w', wall = 'u')] = 0.0

    return SimulationCondition{T}(
        init_cond,
        init_kick,
        force_type,
        force_magnitude,
        bound_cond,
    )

end

function planePoiseuilleFlow(
    init_cond::Char,
    init_kick::T,
    force_type::Char,
    force_magnitude::T,
) where {T<:AbstractFloat}

    # Validity Check
    check_option_validity(init_cond::Char, force_type::Char)

    # Initial Condition
    if (init_cond == 'c')
        throw(
            ArgumentError(
                "PlanePoiseuilleFlow flow configuration not compatible with init_cond argument.",
            ),
        )
    end

    init_kick::T = init_kick

    # Boundary Condition
    bound_cond::Dict{BoundTuple,T} = Dict()

    # u
    bound_cond[(vel = 'u', wall = 'l')] = 0.0
    bound_cond[(vel = 'u', wall = 'u')] = 0.0

    # v
    bound_cond[(vel = 'v', wall = 'l')] = 0.0
    bound_cond[(vel = 'v', wall = 'u')] = 0.0

    # w
    bound_cond[(vel = 'w', wall = 'l')] = 0.0
    bound_cond[(vel = 'w', wall = 'u')] = 0.0

    return SimulationCondition{T}(
        init_cond,
        init_kick,
        force_type,
        force_magnitude,
        bound_cond,
    )

end

end
