module Thermomodel

export zhang2002model!,dMdtzhang2002model

using ..Systems,..Tools

"""
    reserved for future use
"""

function zhang2002model!(du,u,p,t)


    (Xp,dXdt0,M)=vectoXM(u)


    numofliquidslug =  Integer( (length(u) - 1)/5 )
    sys0 = p

    γ = sys0.liquidslug.γ
    ω0 = sys0.liquidslug.ω0
    ℘ = sys0.liquidslug.℘
    Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L)
#     Lliquidslug = XptoLliquidslug(Xp)
    height = getheight(Xp,sys0.tube.L2D,sys0.tube.alpha)
    Xpvapor = getXpvapor(Xp,sys0.tube.L)



    # get P from M and γ
    P = (M./Lvaporplug).^(γ)

    # get θ from P and γ
    θ = P.^((γ-1)/γ)


    for i = 1:numofliquidslug
        du[2*i-1] = u[2*numofliquidslug+2*i-1]
        du[2*i] = du[2*i-1]

        du[2*numofliquidslug + 2*i-1] = -32*u[2*numofliquidslug + 2*i-1] - (ω0[i]^2)*(0.5*(height[i][end]-height[i][1])) + ℘[i]*(P[i]-P[i+1])
        du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

    end



        du[4*numofliquidslug+1:length(u)] .= dMdtzhang2002model(Xpvapor,θ,sys0)

    return du

end

function dMdtzhang2002model(Xpvapor,θ,sys0)

    dMdt=zeros(length(Xpvapor))


    Xe = sys0.evaporator.Xe
    He = sys0.evaporator.He
    θe = sys0.evaporator.θe

    Xc = sys0.condenser.Xc
    Hc = sys0.condenser.Hc
    θc = sys0.condenser.θc

    Levapoverlap=XptoLoverlap(Xpvapor,Xe)
    Lcondoverlap=XptoLoverlap(Xpvapor,Xc)

    for i = 1:length(Xpvapor)
        if Lcondoverlap[i] < 1e-8
            dMdt[i] = He * (Xpvapor[i][end]-Xpvapor[i][1]) * (θe-θ[i])
        else
            dMdt[i] = -Hc * Lcondoverlap[i] * (θ[i]-θc)
        end
    end

    return dMdt

end
end
