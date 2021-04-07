module Thermomodel

export zhang2002model!,dMdtzhang2002model

using ..Systems,..Tools

"""
    This function is required by "DifferentialEquation.jl" Package.
        du :: an empty state vector to be derived at the end of this function.
        u  :: the state vector input from the last time step.
        p  :: some important parameters which do not belong to u
        t  :: the current time

    This approach solves the same set of dimensionless governing equations as
    "https://doi.org/10.1016/S0017-9310(01)00348-9"(Zhang et al. 2002)
    Instead of solving them piece by piece. Our current approach solve them as a whole system.
"""

function zhang2002model!(du::Array{Float64,1},u::Array{Float64,1},p::PHPSystem,t::Float64)


    (Xp,dXdt0,M)=vectoXM(u)


    numofliquidslug =  Integer( (length(u) - 1)/5 )
    sys0 = p

    γ = sys0.liquidslug.γ
    ω0 = sys0.liquidslug.ω0
    ℘ = sys0.liquidslug.℘
    Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L)
#     Lliquidslug = XptoLliquidslug(Xp)
    height = getheight(Xp,sys0.tube.L2D,sys0.tube.angle)
    Xpvapor = getXpvapor(Xp,sys0.tube.L)



    # get P from M and γ
    # P = (M./Lvaporplug).^(γ)
    # P = zeros(size(M))
    # for i in length(M)
    #        P[i] = M[i] > 0 ? (M[i]./Lvaporplug[i]).^(γ) : -(-M[i]./Lvaporplug[i]).^(γ)
    #    end
    P = real.((M./Lvaporplug .+ 0im).^(γ))

    # get θ from P and γ
    # θ = zeros(size(P))
    # for i in length(P)
    #        θ[i] = P[i] > 0 ? P[i].^((γ-1)/γ) : -(-P[i]).^((γ-1)/γ)
    #    end
    θ = real.((P .+ 0im).^((γ-1)/γ))
    # θ = P.^((γ-1)/γ)


    for i = 1:numofliquidslug
        du[2*i-1] = u[2*numofliquidslug+2*i-1]
        du[2*i] = du[2*i-1]

        du[2*numofliquidslug + 2*i-1] = -32*u[2*numofliquidslug + 2*i-1] - (ω0[i]^2)*(0.5*(height[i][end]-height[i][1])) + ℘[i]*(P[i]-P[i+1])
        du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

    end


# not sure which one is better
        du[4*numofliquidslug+1:5*numofliquidslug+1] .= dMdtzhang2002model(Xpvapor,θ,sys0)
        # du[4*numofliquidslug+1:5*numofliquidslug+1] .= dMdtconstantH(Xpvapor,θ,sys0)

    return du

end

"""
    This function is a sub-function of zhang2002model. This function is the phase change part,
        Xpvapor ::  storing the location of two ends of each vapor
        θ       ::  the temperature in each vapor.
        sys0    ::  the struct that stores every needed initial conditions and boundary conditions

    In this case. He and Hc are preset constants.
"""

function dMdtzhang2002model(Xpvapor::Array{Tuple{Float64,Float64},1},θ::Array{Float64,1},sys0::PHPSystem)

    dMdt=zeros(length(Xpvapor))


    Xe = sys0.evaporator.Xe
    He = sys0.evaporator.He
    θe = sys0.evaporator.θe

    Xc = sys0.condenser.Xc
    Hc = sys0.condenser.Hc
    θc = sys0.condenser.θc

    Levapoverlap=XpvaportoLoverlap(Xpvapor,Xe)
    Lcondoverlap=XpvaportoLoverlap(Xpvapor,Xc)


    # May not be right for multi liquid flow
    for i = 1:length(Xpvapor)
        if Lcondoverlap[i] < 1e-8
            dMdt[i] = He * Levapoverlap[i] * (θe-θ[i])
        else
            dMdt[i] = -Hc * Lcondoverlap[i] * (θ[i]-θc)
        end
    end

    return dMdt

end

function dMdtconstantH(Xpvapor::Array{Tuple{Float64,Float64},1},θ::Array{Float64,1},sys0::PHPSystem)

    dMdt=zeros(length(Xpvapor))


    Xe = sys0.evaporator.Xe
    He = sys0.evaporator.He
    θe = sys0.evaporator.θe

    Xc = sys0.condenser.Xc
    Hc = sys0.condenser.Hc
    θc = sys0.condenser.θc

    Levapoverlap=XpvaportoLoverlap(Xpvapor,Xe)
    Lcondoverlap=XpvaportoLoverlap(Xpvapor,Xc)

    for i = 1:length(Xpvapor)
        dMdt[i] = He * Levapoverlap[i] * (θe-θ[i]) -Hc * Lcondoverlap[i] * (θ[i]-θc)
    end
    return dMdt

end

end
