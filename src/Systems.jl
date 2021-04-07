module Systems

export PHPSystem,Tube,Evaporator,Condenser,Liquid,Vapor,Wall,PHPResult

# using ..Tools

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω0
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Tube
    L::Float64
    L2D::Float64
    angle::Float64
end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω0
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Evaporator
    He::Float64
    θe::Float64
    Xe::Array{Tuple{Float64,Float64},1}
end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω0
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Condenser
    Hc::Float64
    θc::Float64
    Xc::Array{Tuple{Float64,Float64},1}
end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω0
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Liquid
    γ::Float64
    ω0::Array{Float64,1}
    ℘::Array{Float64,1}
    Xp::Array{Tuple{Float64,Float64},1}
    dXdt::Array{Tuple{Float64,Float64},1}
    Xarrays::Array{Array{Float64,1},1}
    θarrays::Array{Array{Float64,1},1}
end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω0
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Vapor
    γ::Float64
    P::Array{Float64,1}
end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω0
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Wall
    α::Float64
    Xarray::Array{Float64,1}
    θarray::Array{Float64,1}
end



"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω0
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""
# add type names!
mutable struct PHPSystem
    tube::Tube
    evaporator::Evaporator
    condenser::Condenser
    liquid::Liquid
    vapor::Vapor
    wall::Wall
end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω0
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

struct PHPResult
    t
    Xp
    dXdt
    P
    θ
    M
end


end
