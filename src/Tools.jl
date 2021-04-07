module Tools

export getheight,XMtovec,vectoXM,XptoLvaporplug,XptoLliquidslug,getXpvapor,XpvaportoLoverlap

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





end
