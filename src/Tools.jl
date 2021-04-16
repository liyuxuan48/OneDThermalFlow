module Tools

export getheight,XMtovec,vectoXM,XptoLvaporplug,XptoLliquidslug,getXpvapor,XpvaportoLoverlap,ifamongone,ifamong,settemperature!,laplacian,constructXarrays,constructXarrays,walltoliquidmapping,liquidtowallmapping,truncate,constructmapping,duliquidθtovec,duwallθtovec,liquidθtovec,wallθtovec,updateXarrays,getcurrentsys,wallmodel,liquidmodel

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

function updateXarrays(Xp,Xarrays)

    for i = 1:length(Xp)
        Xarrays[i] .= Xarrays[i] .- [Xarrays[i][1]] .+ [Xp[i][1]]
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


    Xp,dXdt,M = vectoXM(u[1:indexes[1]-1])
    θwallrec = u[indexes[1]+1:indexes[2]-1]

    for i = 1:length(indexes)-2
    push!(θliquidrec, u[indexes[i+1]+1:indexes[i+2]-1])
    end
    push!(θliquidrec, u[indexes[end]+1:end])

    sysnew = deepcopy(sys0)

    sysnew.liquid.Xp = Xp
    sysnew.liquid.dXdt = dXdt
    sysnew.liquid.Xarrays = updateXarrays(Xp,sys0.liquid.Xarrays)
    sysnew.liquid.θarrays = θliquidrec

    Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L)
    γ = sys0.vapor.γ
    P = real.((M./Lvaporplug .+ 0im).^(γ))
    sysnew.vapor.P = P

    sysnew.wall.θarray = θwallrec

    walltoliquid, liquidtowall = constructmapping(sysnew.liquid.Xarrays ,sysnew.wall.Xarray)
    sysnew.mapping = Mapping(walltoliquid,liquidtowall)

    return sysnew
end

function wallmodel(θarray::Array{Float64,1},p::PHPSystem)
    sys = p

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
    sys = p

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






end
