module Postprocessing

export soltoResult

using ..Systems,..Tools

"""
    reserved for future use
"""

function soltoResult(sol,sys0) #only good for one calculation per time point, not good for onec calculation for all time

    γ = sys0.liquidslug.γ
    numofliquidslug =  Integer( (size(sol)[1]-1)/5  )

    MatrxXp=sol[1:2*numofliquidslug,:]
    MatrxdXdt=sol[2*numofliquidslug+1:4*numofliquidslug,:]
    M=sol[4*numofliquidslug+1:end,:]

    Xp = Array{Tuple{Float64, Float64}}((undef), Integer(size(MatrxXp)[1]/2), size(MatrxXp)[2])
    dXdt = Array{Tuple{Float64, Float64}}((undef), size(Xp))
    Lvaporplug = zeros(size(M))

    # transfer from 2D array to array of tuple for Xp and dXdt
    for i in 1:size(Xp)[1]
        for j in 1:size(Xp)[2]
            Xp[i,j]=((MatrxXp[2*i-1,j]),(MatrxXp[2*i,j]))
            dXdt[i,j]=((MatrxdXdt[2*i-1,j]),(MatrxdXdt[2*i,j]))
        end
    end

    for j in 1:size(M)[2]
        Lvaporplug[:,j]=XptoLvaporplug(Xp[:,j], sys0.tube.L)
    end

    # get P from M and γ
    P = (M./Lvaporplug).^(γ)

    # get θ from P and γ
    θ = P.^((γ-1)/γ)

    # convert a matrix to array of array
    Xp=mapslices(x->[x], MatrxXp, dims=2)[:]
    dXdt=mapslices(x->[x], MatrxdXdt, dims=2)[:]
    M=mapslices(x->[x], M, dims=2)[:]
    P=mapslices(x->[x], P, dims=2)[:]
    θ=mapslices(x->[x], θ, dims=2)[:]

    result=PHPResult(sol.t,Xp,dXdt,P,θ,M)
end
end
