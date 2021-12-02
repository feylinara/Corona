module Corona

using DynamicalSystems, LabelledArrays
# The compartmental model implemented here is

# \\begin{align*}
# \\frac{\\mathrm{d}S}{\\mathrm{d}} &= -\\frac{β S I}{N} - v S + \\overline r R + \\overline v V\\\\
# \\frac{\\mathrm{d}I}{\\mathrm{d}} &= \frac{\beta S I}{N} - (r + d) I\\\\
# \\frac{\\mathrm{d}R}{\\mathrm{d}} &= rI - \\overline r R\\\\
# \\frac{\\mathrm{d}V}{\\mathrm{d}} &= v S - \\overline v V\\\\
# \\frac{\\mathrm{d}D}{\\mathrm{d}} &= d I\\\\
# \\end{align*}

# with ``S`` being the sucseptible population, ``I`` infected, ``R`` recovered, ``V``
# vaccinated, and ``D`` the deaths. To better model current decision making it
# might be useful to add ``V_{boost}`` for people who have recieved their booster
# shot, and ``V_{waned}`` for those who are due to recieve it, rather than
# releasing them to ``S\. It might also be a good idea to revise the vaccination
# function to a logistic one. It might also be useful to implement a delay, but
# I do not know whether this is possible with `DynamicalSystems.jl` and do not
# have access to powersim
function SIRVD(x, p, t)
    β = p.β(t)
    recovery = p.recovery(t)
    deaths = p.deaths(t)
    vax = p.vax(t)
    resucseptible = p.resucseptible(t)
    vaxwane = p.wane(t)
    immunity = p.immunity(t)
    population = p.population
    S = x[1]
    I = x[2]
    R = x[3]
    V = x[4]
    Vwaned = x[5]
    Vboost = x[6]
    D = x[7]
    SVector{7}(
        (-β * S * I) / population - vax * S / population + resucseptible * R,
        (β * S * I) / population + β * Vwaned * I / population + β * immunity * V * I / population - (recovery + deaths) * I,
        recovery * I - resucseptible * R,
        vax * S / population - β * (1.0 - immunity) * V * I / population - vaxwane * V,
        vaxwane * V - β * Vwaned * I / population,
        0.0,
        deaths * I,
    )
end

"""
"""
function model(;
    IFR = 0.0041,
    R₀ = 2,
    vax = 1_000_000,
    infectionPeriod = 12,
    recoveryImmunity = 6 * 30,
    vaxStart = 365 * 1.5,
    vaxImmunity = 6 * 30,
    effectiveness = 0.7,
    boosterEffectiveness = 0.7)

    pop = 80000000
    s0 = [pop - 1, 1, 0.00, 0.00, 0.00, 0.00, 0.00]
    β = R₀ * (IFR + 1 / infectionPeriod)
    p = SLVector(
        β = t -> β,
        recovery = t -> 1 / infectionPeriod,
        deaths = t -> IFR,
        vax = t -> if t > vaxStart
            vax
        else
            0
        end,
        resucseptible = t -> 1 / (recoveryImmunity * 30),
        wane = t -> 1 / (vaxImmunity * 30),
        immunity = t -> effectiveness,
        population = pop
    )
    ContinuousDynamicalSystem(Corona.SIRVD, s0, p)
end
end