module Boiling

using ..Systems,..Tools

export nucleateboiling,getnewδ,crossAtoδ,getcrossAδ,getnewθarrays,getnewXarrays,getinsertindex,getarrayindex,getnewXp,getnewδ,getnewM


function nucleateboiling(sys,Xvapornew,Pinsert)
    ρ = deepcopy(sys.liquid.ρ)
    d = deepcopy(sys.tube.d)
    Xp = deepcopy(sys.liquid.Xp)
    dXdt = deepcopy(sys.liquid.dXdt)
    δ = deepcopy(sys.vapor.δ)
    P = deepcopy(sys.vapor.P)
    L = sys.tube.L
    γ = sys.vapor.γ
    Xarrays = sys.liquid.Xarrays
    θarrays = sys.liquid.θarrays

    Lvaporplug =    XptoLvaporplug(Xp,sys.tube.L)
    M = P.^(1/γ).* Lvaporplug







    # Xvapornew = (1.50,1.60)
    # Pinsert = 1.5
    ρinsert = Pinsert.^(1/γ)
    Linsert = Xvapornew[end] - Xvapornew[1]
    Minsert = Pinsert.^(1/γ).* Linsert

    index = getinsertindex(Xp,Xvapornew)
    δnew = getnewδ(δ,index,Xvapornew,ρ,d,Minsert) # mass conservation
    Xpnew = getnewXp(Xp,index,Xvapornew)
    Mnew = getnewM(M,index,Minsert)

    Lvaporplugnew = XptoLvaporplug(Xpnew,L)
    Pnew = (Mnew./Lvaporplugnew).^γ

    Xarraysnew = getnewXarrays(index,Xp,Xpnew,Xarrays,L)
    θarraysnew = getnewθarrays(index,Xp,Xpnew,Xarrays,θarrays,L)

    dXdtnew = deepcopy(dXdt) # momentum conservation
    insert!(dXdtnew,index+1,dXdtnew[index])






    sysnew = deepcopy(sys)

    sysnew.liquid.Xp = Xpnew
    sysnew.liquid.dXdt = dXdtnew
    sysnew.liquid.Xarrays = Xarraysnew
    sysnew.liquid.θarrays = θarraysnew
    sysnew.vapor.P = Pnew
    sysnew.vapor.δ = δnew

    walltoliquid, liquidtowall = constructmapping(sysnew.liquid.Xarrays ,sysnew.wall.Xarray)
    sysnew.mapping = Mapping(walltoliquid,liquidtowall)

return sysnew
end


function getnewδ(δ,index,Xvapornew,ρ,d,Minsert)
    Linsert = Xvapornew[end] - Xvapornew[1]

    crossAfilms = getcrossAδ.([d],δ)
    insertcrossA = 0.5(crossAfilms[index] + crossAfilms[index+1]) - Minsert/ρ/Linsert
    insert!(crossAfilms,index+1,insertcrossA)

    δnew = crossAtoδ.([d],crossAfilms)
end


function crossAtoδ(d,crossA)
    C = crossA # nondimensional area
    δ = 1 - sqrt(1-C)
end


function getcrossAδ(d,δ)
    crossAouter = 1
    crossAinner = (1-δ)^2
    crossAδ     = crossAouter - crossAinner
end



 function getnewθarrays(index,Xp,Xpnew,Xarrays,θarrays,L)
    θarraysnew = deepcopy(θarrays)
    Xarraysnew = deepcopy(Xarrays)
    arrayindex = getarrayindex(Xpnew[index][2],Xarrays[index])

    θarraysnewleft = θarrays[index][1:arrayindex]
    θarraysnewright= θarrays[index][arrayindex+1:end]

    splice!(θarraysnew, index)
    insert!(θarraysnew, index,θarraysnewleft)
    insert!(θarraysnew, index+1,θarraysnewright)
end


function getnewXarrays(index,Xp,Xpnew,Xarrays,L)
    Xarraysnew = deepcopy(Xarrays)
    arrayindex = getarrayindex(Xpnew[index][2],Xarrays[index])

    Xarraysnewleft = LinRange(Xpnew[index][1],Xpnew[index][2],arrayindex)
    Xarraysnewright= LinRange(Xpnew[index+1][1],Xpnew[index+1][2],length(Xarrays[index])-arrayindex)

    splice!(Xarraysnew, index)
    insert!(Xarraysnew, index,Xarraysnewleft)
    insert!(Xarraysnew, index+1,Xarraysnewright)
end


function getinsertindex(Xp,Xvapornew)

for index = 1:length(Xp)
    if Xp[index][1] <= Xvapornew[1] && Xp[index][2] >= Xvapornew[2]
        return index
end
    end
        return NaN
end

function getarrayindex(X,Xarray)

for arrayindex = 1:length(Xarray)
    if X >= Xarray[arrayindex] && X <= Xarray[arrayindex+1]
        return arrayindex
end
    end
        return NaN
end

function getnewXp(Xp,index,Xvapornew)

    Xpnew = deepcopy(Xp)

    Linsert = Xvapornew[end] - Xvapornew[1]

    insertXp1=(Xp[index][1]-Linsert/2,Xvapornew[1])
    insertXp2=(Xvapornew[2],Xp[index][2]+Linsert/2)

    splice!(Xpnew, index)
    insert!(Xpnew, index,insertXp1)
    insert!(Xpnew, index+1,insertXp2)

    return Xpnew
end

function getnewδ(δ,Pinsert,index,Xvapornew)

    Linsert = Xvapornew[end] - Xvapornew[1]
    δnew = deepcopy(δ)

    δavg = (δ[index]+δ[index+1])/2
    insert!(δnew, index+1,δavg)

    return δnew
end



function getnewM(M,index,Minsert)

    Mnew = deepcopy(M)

    insert!(Mnew,index+1,Minsert)

    return Mnew
end


end
