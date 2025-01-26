module IntegrationSuite

#=
Package
=#

#=
Public API
=#

export integral_y

#=
Integration
=#

function integral_y(arr::Vector{Complex{T}}, dy::Vector{T})::T where {T<:AbstractFloat}

    # Integral
    integral::Complex{T} = 0.0

    # Loop
    for i = 1:length(dy)

        integral += dy[i] * (arr[i+1] + arr[i]) / 2.0

    end

    return real(integral)

end

end
