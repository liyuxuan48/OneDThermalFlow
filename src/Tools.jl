module Tools

export getheight

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

function getheight(sys)
    Xp=sys.liquidslug.Xp
    L2D=sys.tube.L2D
    alpha=sys.tube.alpha

    height=deepcopy(Xp)

    for i =1:length(Xp)
        height[i]=(getoneheight(Xp[i][1],L2D,alpha), getoneheight(Xp[i][end],L2D,alpha))
    end


    return height
end

end
