module TimeMarching

export

using ..Systems

"""
    this function's inputs are t and uu

        t is the current time

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is ̂E (total energy per volume)

    this function's outputs are new t and uunew

        t is the current time

        uunew has the same structure as uu

    this function uses flux splitting scheme (Steger and Warming, 1981)
    and 1st order Euler time marching

    The current function is only good for u<=c (speed of sound)
"""


end
