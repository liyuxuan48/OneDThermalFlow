module Tools

export getheight,XMtovec,vectoXM,XptoLvaporplug,XptoLliquidslug,getXpvapor,XptoLoverlap,oneoverlap,ifoverlap

"""
    this function's inputs are uu and gamma

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is ̂E (total energy per volume)

        gamma is the heat capacity ratio

    this function's outputs is tourhoc

        tourhoc has 3 rows
        the 1st row is u (velocity)
        the 2nd row is ρ (density)
        the 3rd row is c (speed of sound)

    this function uses perfect gas state function to get tourhoc
"""
function getoneheight(X,L2D,alpha)

    Integer(mod(div(X,L2D),2.0)) == 0 ? L2D - mod(X,L2D) : mod(X,L2D)

end

"""
    this function's inputs are uu and gamma

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is ̂E (total energy per volume)

        gamma is the heat capacity ratio

    this function's outputs is tourhop

        tourhoc has 3 rows
        the 1st row is u (velocity)
        the 2nd row is ρ (density)
        the 3rd row is p (pressure)

    this function uses perfect gas state function to get tourhop
"""

function getheight(Xp,L2D,alpha)

    height=deepcopy(Xp)

    for i =1:length(Xp)
        height[i]=(getoneheight(Xp[i][1],L2D,alpha), getoneheight(Xp[i][end],L2D,alpha))
    end


    return height
end

function XMtovec(Xp,dXdt,M)
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

function vectoXM(u)

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

function XptoLvaporplug(Xp,L)

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

function XptoLliquidslug(Xp)

    Lliquidslug = zeros(length(Xp))


        for i = 1: length(Xp)

        Lliquidslug[i] = Xp[i][end] - Xp[i][1]

        end

    return Lliquidslug

end

function getXpvapor(Xp,L)
    Xpvapor=deepcopy(Xp)

    Xpvapor[1]=(0.0,Xp[1][1])

    for i = 2:(length(Xp))
        Xpvapor[i]=(Xp[i-1][end],Xp[i][1])
    end

    push!(Xpvapor,(Xp[end][end],L))

    return Xpvapor
end

function XptoLoverlap(Xpvapor,Xe)
    Loverlap = zeros(length(Xpvapor))
    for i = 1:length(Xpvapor)
        for j = 1:length(Xe)
            if ifoverlap(Xpvapor[i],Xe[j])
            Loverlap[i] += oneoverlap(Xpvapor[i],Xe[j])
            end
        end
    end

    return Loverlap
end

function oneoverlap(oneXpvapor,oneXe)
    return min(oneXpvapor[end],oneXe[end]) - max(oneXpvapor[1],oneXe[1])
end

function ifoverlap(oneXpvapor,oneXe)
    return ( (oneXpvapor[end] >= oneXe[1]) || (oneXpvapor[1] <= oneXe[end]) ) && (oneXpvapor[end] > oneXe[1]) && (oneXpvapor[1] < oneXe[end])
end





end
