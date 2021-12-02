import Pkg; Pkg.activate(".");
import Corona
using WebIO, Interact, Blink, DynamicalSystems, Plots, LabelledArrays

Plots.set_theme(:juno)

app = @manipulate for IFR = widget(0:0.0005:0.01; label = "Deaths per Infection", value = 0.004),
    R₀ = widget(0.5:0.1:7; label = "Base R₀"),
    vax = widget(10_000:10_000:2_000_000; label = "Maximum Vaccinations per Day"),
    vaxstart = widget(0:30:600; label = "Start of Vaccinations (days after first infection)"),
    period = widget(3:1:14, label = "Length of Infection"),
    recoveryImmunity = widget(3:1:12, label = "Duration of Immunity after Infection"),
    vaxImmunity = widget(6:1:24, label = "Duration of Immunity after Vaccination")


    t = 365 * 3
    model = Corona.model(
        IFR = IFR[],
        R₀ = R₀[],
        vax = vax[],
        infectionPeriod = period[],
        recoveryImmunity = recoveryImmunity[],
        vaxStart = vaxstart[],
        vaxImmunity = vaxImmunity[] * 30,
        effectiveness = 0.7,
        boosterEffectiveness = 0.7
    );
    traj = trajectory(model, t; Δt = 1) |> Matrix;
    vbox(
        "After $t days, the Simulated Virus Disease has killed $(Int(floor(traj[size(traj)[1], 5]))) people",
        plot(traj;
            labels = ["Sucsptible" "Infected" "Recovered" "Vaccinated" "Booster needed" "Boostered" "Dead"],
            linewidth = 2
        )
    )
end

w = Window()
title(w, "Corona Simulation")
body!(w, app);