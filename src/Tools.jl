module Tools

export getheight,XMtovec,XMδtovec,vectoXM,vectoXMδ,XptoLvaporplug,XptoLliquidslug,getXpvapor,XpvaportoLoverlap,ifamongone,ifamong,settemperature!,laplacian,constructXarrays,constructXarrays,walltoliquidmapping,liquidtowallmapping,truncate,constructmapping,duliquidθtovec,duwallθtovec,liquidθtovec,wallθtovec,updateXarrays,getcurrentsys,wallmodel,liquidmodel,nucleateboiling,getnewδ,crossAtoδ,getcrossAδ,getnewθarrays,getnewXarrays,getinsertindex,getarrayindex,getnewXp,getnewδ,getnewM,newvaporbubble

using ..Systems
using LinearAlgebra

"""
    This function is a sub-function of getheight. This function is to get the actural physical height for one interface
        X     ::   the location of one interface
        L2D   ::   the length of one bend to another bend (the length in 2D)
        angle ::   the inclination angle
"""

function getoneheight(X::Float64,L2D::Float64,angle::Float64)

    Integer(mod(div(X,L2D),2.0)) == 0 ? L2D - mod(X,L2D) : mod(X,L2D)

end

"""
    This function is to get the actural physical heights for all interfaces
        Xp    ::   the locations of all interfaces
        L2D   ::   the length of one bend to another bend (the length in 2D)
        angle ::   the inclination angle
"""

function getheight(Xp::Array{Tuple{Float64,Float64},1},L2D::Float64,angle::Float64)

    height=deepcopy(Xp)

    for i =1:length(Xp)
        height[i]=(getoneheight(Xp[i][1],L2D,angle), getoneheight(Xp[i][end],L2D,angle))
    end

    return height
end

"""
    This function is to transform Xp, dXdt of the interface, and M of the vapor to form our state vector u
        Xp    ::   the locations of all interfaces
        dXdt  ::   the 1D velocity of all interfaces
        M     ::   the mass of all vapors
"""

function XMtovec(Xp::Array{Tuple{Float64,Float64},1},dXdt::Array{Tuple{Float64,Float64},1},M::Array{Float64,1})
    if (length(Xp) == length(dXdt)) && (length(Xp) + 1 == length(M))

        u=zeros(5*length(Xp)+1)

        for i = 1:length(Xp)

            # input Xp
            u[2*i-1] = Xp[i][1]
            u[2*i] = Xp[i][end]

            # input dXdt
            u[2*length(Xp) + 2*i-1] = dXdt[i][1]
            u[2*length(Xp) + 2*i] = dXdt[i][end]
        end

        for i = 1:length(M)
            # input M
            u[4*length(Xp) + i] = M[i]
        end
    else
        println("the lengthes of X and dXdt and M do not match!")
    end

    return u

end

function XMδtovec(Xp,dXdt,M,δ)

    return ([XMtovec(Xp,dXdt,M);δ])
end

"""
    This function is to transform Xp, dXdt of the interface, and M of the vapor to form our state vector u
        u    ::   the state vector
"""

function vectoXM(u::Array{Float64,1})

    maxindex = Integer( (length(u) - 1)/5 )

    Xp = map(tuple, zeros(maxindex), zeros(maxindex))
    dXdt = map(tuple, zeros(maxindex), zeros(maxindex))
    M = zeros(maxindex+1)

    for i = 1:maxindex

        # input Xp
        Xp[i] = (u[2*i-1],u[2*i])

        # input dXdt
        dXdt[i] = (u[2*maxindex + 2*i-1],u[2*maxindex + 2*i])
    end

    for i = 1:(maxindex+1)

        # input M
        M[i] = u[4*maxindex + i]

    end

    return Xp,dXdt,M

end

"""
    This function is to transform Xp, dXdt of the interface, and M of the vapor to form our state vector u
        u    ::   the state vector
"""

function vectoXMδ(u::Array{Float64,1})

    maxindex = Integer( (length(u) - 2)/6 )

    Xp = map(tuple, zeros(maxindex), zeros(maxindex))
    dXdt = map(tuple, zeros(maxindex), zeros(maxindex))
    M = zeros(maxindex+1)
    δ = zeros(maxindex+1)

    for i = 1:maxindex

        # input Xp
        Xp[i] = (u[2*i-1],u[2*i])

        # input dXdt
        dXdt[i] = (u[2*maxindex + 2*i-1],u[2*maxindex + 2*i])
    end

    for i = 1:(maxindex+1)

        # input M
        M[i] = u[4*maxindex + i]
        δ[i] = u[5*maxindex + 1 + i]
    end

    return Xp,dXdt,M,δ

end

"""
    This function is to transform Xp of every interface, and L of the tube to form an array of vapor length
        Xp    ::   the locations of all interfaces
        L     ::   the length of the 1D tube
"""

function XptoLvaporplug(Xp::Array{Tuple{Float64,Float64},1},L::Float64)

    maxindex = length(Xp) + 1
    Lvaporplug = zeros(maxindex)

    Lvaporplug[1] = Xp[1][1]-0.0
    Lvaporplug[end] = L-Xp[end][end]

    if maxindex > 2
        for i = 2:maxindex-1

            Lvaporplug[i] = Xp[i][1] - Xp[i-1][end]

        end
    end

    return Lvaporplug

end

"""
    This function is to transform Xp of every interface to form an array of vapor length
        Xp    ::   the locations of all interfaces
"""

function XptoLliquidslug(Xp::Array{Tuple{Float64,Float64},1})

    Lliquidslug = zeros(length(Xp))


        for i = 1: length(Xp)

        Lliquidslug[i] = Xp[i][end] - Xp[i][1]

        end

    return Lliquidslug

end

"""
    The Xp was coupled by every liquid slug. For instance, if there is one liquid slug. Xp is a one-element tuple (Xp[1][1], Xp[1][2]).
    But sometimes we need Xp to be coupled by every vapor plug. For one liquid slug, we have two vapor plugs.
    So by adding 0 and L at the beginning and the end,
    we construct a two-element tuple ((0.0,Xp[1][1]) and ((Xp[1][2],L). Generally, for every N-element Xp, we construct an N+1 element Xpvapor
        Xp    ::   the locations of all interfaces, each element means a liquid slug.
        L     ::   the length of the 1D tube
"""

function getXpvapor(Xp,L)
    Xpvapor=deepcopy(Xp)

    Xpvapor[1]=(0.0,Xp[1][1])

    for i = 2:(length(Xp))
        Xpvapor[i]=(Xp[i-1][end],Xp[i][1])
    end

    push!(Xpvapor,(Xp[end][end],L))

    return Xpvapor
end

"""
    This function aims to get the overlapping length between the vapor plug and each evaporator/condenser section and sum them up.
    It can sum up all evaporator sections or condensor sections respectively, but it cannot sum both of them at the same time.
    For example:

    XpvaportoLoverlap(Xpvapor,Xe) sums all the overlapping regions for evaporator
    XpvaportoLoverlap(Xpvapor,Xc) sums all the overlapping regions for condensor

        Xpvapor    ::   the locations of all interfaces, each element means a vapor plug.
        Xce        ::   the locations of all evaporators/condensors, each element is a evaporators/condensors section
"""

function XpvaportoLoverlap(Xpvapor,Xce)
    Loverlap = zeros(length(Xpvapor))
    for i = 1:length(Xpvapor)
        for j = 1:length(Xce)
            if ifoverlap(Xpvapor[i],Xce[j])
            Loverlap[i] += oneoverlap(Xpvapor[i],Xce[j])
            end
        end
    end

    return Loverlap
end

"""
    This is a sub-function to solve the overlapping length between a single vapor section and a single evaporator/condensor section

    oneXpvapor ::  the locations of two ends of a vapor plug
    oneXce     ::  the locations of two ends of a evaporator/condensor
"""

function oneoverlap(oneXpvapor,oneXce)
    return min(oneXpvapor[end],oneXce[end]) - max(oneXpvapor[1],oneXce[1])
end

"""
    This is a sub-function to determine if there is a overlapping region between a single vapor section and a single evaporator/condensor section

    oneXpvapor ::  the locations of two ends of a vapor plug
    oneXce     ::  the locations of two ends of a evaporator/condensor
"""

function ifoverlap(oneXpvapor,oneXce)
    return ( (oneXpvapor[end] >= oneXce[1]) || (oneXpvapor[1] <= oneXce[end]) ) && (oneXpvapor[end] > oneXce[1]) && (oneXpvapor[1] < oneXce[end])
end

function ifamongone(value::Float64, range::Tuple{Float64,Float64})
    return (value >= range[1]) && (value <= range[2]) ? true : false
end

function ifamong(value::Float64, X::Array{Tuple{Float64,Float64},1})

    return Bool(sum(ifamongone.(value,X)))
end

function settemperature!(θᵣ,xvalue,sys0)

    if ifamong(xvalue, sys0.evaporator.Xe)
        θᵣ = sys0.evaporator.θe

    elseif ifamong(xvalue, sys0.condenser.Xc)
        θᵣ = sys0.condenser.θc
    end

    return θᵣ

end

function laplacian(u)
    unew = deepcopy(u)

    dl = ones(length(u)-1)
    dr = dl
    d  = -2*ones(length(u))

    A = Tridiagonal(dl, d, dr)

    unew = A*u

    unew[1]   = unew[2]
    unew[end] = unew[end-1]

    return (unew)
end

function constructXarrays(X0::Array{Tuple{Float64,Float64},1},N,θinitial,L)
    Xarrays=Array{Array{Float64, 1}, 1}(undef, length(X0))

    Lliquid = XptoLliquidslug(X0)

    Nliquid =  floor.(Int, N.*Lliquid./L)

    for i = 1:length(Xarrays)
        Xarrays[i] = range(X0[i][1], X0[i][2], length=Nliquid[i])
    end

    θarrays = deepcopy(Xarrays)
    for i = 1:length(θarrays)
        θarrays[i][:] .= θinitial
    end

    return(Xarrays,θarrays)
end

function constructXarrays(L::Float64,N,θinitial)
    Xwallarray = Array{Float64, 1}(undef, N)
    Xwallarray = range(0, L, length=N)

    θwallarray = deepcopy(Xwallarray)
    θwallarray = range(θinitial, θinitial, length=N)

    return(Xwallarray,θwallarray)
end

function walltoliquidmapping(Xwall,Xarrays)
for i = 1:length(Xarrays)
    if Xarrays[i][end] < Xwall

    else
        for j = 1:length(Xarrays[i])
            if j == 1 && Xarrays[i][j] >= Xwall
                    return (i,-1)
                    elseif Xarrays[i][j] >= Xwall && Xarrays[i][j-1] <= Xwall
                        return (i,j)
            end
        end
    end
end
    return (length(Xarrays)+1,-1) # for closed end tube
end



function liquidtowallmapping(Xliquidone,Xwallarray)

for i = 2:length(Xwallarray)
    if Xwallarray[i] >= Xliquidone && Xwallarray[i-1] <= Xliquidone
        return (i)
    end
end
    return (-1) # for closed end tube
end


function truncate(Xarrays::Array{Array{Float64,1},1})

    integerXarrays = Array{Array{Int64,1},1}(undef, length(Xarrays))

    for i =1:length(Xarrays)
        integerXarrays[i] = trunc.(Int, Xarrays[i])
    end
    return integerXarrays
end

function constructmapping(Xarrays,Xwallarray)
    walltoliquid = Array{Tuple{Int64,Int64},1}(undef, length(Xwallarray))

    for i = 1:length(Xwallarray)
        walltoliquid[i] = walltoliquidmapping(Xwallarray[i],Xarrays)
    end

    liquidtowall = truncate(Xarrays)

    for i = 1:length(Xarrays)
        for j = 1:length(Xarrays[i])
            liquidtowall[i][j] = liquidtowallmapping(Xarrays[i][j],Xwallarray)
        end
    end

    return walltoliquid,liquidtowall
end

function duliquidθtovec(duθarrays)
    return vcat(map(duwallθtovec, duθarrays)...)
end

function duwallθtovec(duθwall)
    return [0.0; duθwall]
end

function liquidθtovec(θarrays)
    return vcat(map(wallθtovec, θarrays)...)
end

function wallθtovec(θwall)
    return [-1e10; θwall]
end

function updateXarrays(Xp,θarrays)

    Xarrays = deepcopy(θarrays)

    for i = 1:length(Xp)
        Xarrays[i] = LinRange(Xp[i][1],Xp[i][2],length(θarrays[i]))
    end

    return Xarrays
end

function getcurrentsys(u,sys0)

        indexes = Int64[]
        θliquidrec = Array[]

        for i = 1:length(u)
            if abs(u[i]+1e10) <= 10^(-1)
                push!(indexes,i)
            end
        end


    Xp,dXdt,M,δ = vectoXMδ(u[1:indexes[1]-1])
    θwallrec = u[indexes[1]+1:indexes[2]-1]

    for i = 1:length(indexes)-2
    push!(θliquidrec, u[indexes[i+1]+1:indexes[i+2]-1])
    end
    push!(θliquidrec, u[indexes[end]+1:end])

    sysnew = deepcopy(sys0)

    sysnew.liquid.Xp = Xp
    sysnew.liquid.dXdt = dXdt
    sysnew.liquid.θarrays = θliquidrec
    sysnew.liquid.Xarrays = updateXarrays(Xp,sysnew.liquid.θarrays)


    Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L)
    γ = sys0.vapor.γ
    P = real.((M./Lvaporplug .+ 0im).^(γ))
    sysnew.vapor.P = P
    sysnew.vapor.δ = δ

    sysnew.wall.θarray = θwallrec

    walltoliquid, liquidtowall = constructmapping(sysnew.liquid.Xarrays ,sysnew.wall.Xarray)
    sysnew.mapping = Mapping(walltoliquid,liquidtowall)

    return sysnew
end

function wallmodel(θarray::Array{Float64,1},p::PHPSystem)
    sys = deepcopy(p)

    du = zero(deepcopy(θarray))

    γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ
    He = sys.evaporator.He
    dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]


    H = zero(deepcopy(θarray))
    θarray_temp_flow = zero(deepcopy(θarray))
    for i = 1:length(θarray)

        index = sys.mapping.walltoliquid[i]

        if index[2] == -1
            P = sys.vapor.P[index[1]]
            θarray_temp_flow[i] = real.((P .+ 0im).^((γ-1)/γ))

            H = He
        else
            θliquidarrays = sys.liquid.θarrays
            θarray_temp_flow[i] = θliquidarrays[index[1]][index[2]]

            H = Hₗ
        end
    end

#     print("θ=",θarray_temp_flow[1:20],"\n")


    du = sys.wall.α .* laplacian(θarray) ./ dx ./ dx + H .* (θarray_temp_flow - θarray) .* dx

#     du = sys.wall.α .* laplacian(θarray) ./ dx ./ dx

    return du
end

function liquidmodel(θarrays,p::PHPSystem)
    sys = deepcopy(p)

    du = zero.(deepcopy(θarrays))

    γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ



    θarray_temp_wall = zero.(deepcopy(θarrays))
    for i = 1:length(θarrays)

        dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]

        indexes = sys.mapping.liquidtowall[i]

        for j = 1:length(indexes)
            θarray_temp_wall[i][j] = sys.wall.θarray[indexes[j]]
        end

        du[i] = sys.wall.α .* laplacian(θarrays[i]) ./ dx ./ dx + Hₗ .* (θarray_temp_wall[i] - θarrays[i]) .* dx
    end


    return du
end

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

function newvaporbubble(sys)
    Xp = sys.liquid.Xp
    dXdt = sys.liquid.dXdt
    P = sys.vapor.P

    Lvaporplug = XptoLvaporplug(X0,sys.tube.L)
    M = P.^(1/γ).* Lvaporplug

    δ = sys.vapor.δ


    index = getinsertindex(Xp,Xvapornew)

    modifyδ!(δ,index,Xvapornew)

    insertXvapor!(Xp,index,Xvapornew)


end



end
