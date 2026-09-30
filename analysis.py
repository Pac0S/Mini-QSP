import numpy as np
import matplotlib.pyplot as plt
from model import simulate, set_parameter, set_parameters, residuals_ic50
from scipy.optimize import least_squares




params = {
    "ka": 1.0,          #Central compartment drug absorbtion rate
    "CL": 2.0,          #Drug clearance
    "V": 20.0,          #Drug gistribution volume
    "Imax": 0.8,        #Maximal Biomarker production inhibition
    "IC50": 0.2,        #Half maximal inhibitory concentration (Drug concentration inducing half inhibition effect)
    "kprod": 10.0,      #Biomarker production rate
    "kdeg": 0.1,        #Biomarker degradation rate
    "dose" : 50.0         #Drug dose
}


#48h time, with 12 mn interval between each discretized observation
times = np.linspace(0, 48, 241)

result = simulate(params=params, times=times)

t = result.t

Ag, Ac, B = result.y
C = Ac / params["V"]
B0 = params["kprod"]/params["kdeg"]

fig, axes = plt.subplots(3,1,figsize=(8,9), sharex=True)

axes[0].plot(t, Ag, label="Ag : Absorption compartment")
axes[0].plot(t, Ac, label="Ac : Central compartment")
axes[0].set_ylabel("Quantité (mg)")
axes[0].legend()

axes[1].plot(t, C, color="tab:orange")
axes[1].set_ylabel("Concentration C (mg/L)")

axes[2].plot(t, B/B0, color="tab:green")
axes[2].axhline(0.7, color="gray", linestyle="--", label="30% Reduction")
axes[2].set_ylabel("Biomarker B / B0")
axes[2].set_xlabel("Time (h)")
axes[2].legend()

fig.tight_layout()
plt.show()

#Test system behavior with 0, 20 and 50 mg as drug dose.

params_0 = set_parameter(params, "dose", 0.0)
params_20 = set_parameter(params, "dose", 20.0)
params_50 = set_parameter(params, "dose", 50.0)

result_0 = simulate(params=params_0, times=times)
result_20 = simulate(params=params_20, times=times)
result_50 = simulate(params=params_50, times=times)

Ag_0, Ac_0, B_0 = result_0.y
Ag_20, Ac_20, B_20 = result_20.y
Ag_50, Ac_50, B_50 = result_50.y

C_0 = Ac_0/params_0["V"]
C_20 = Ac_20/params_20["V"]
C_50 = Ac_50/params_50["V"]

fig, axes = plt.subplots(2,1,figsize=(8,9), sharex=True)
# Plot concentration curves
axes[0].plot(t, C_0, label="C with D = 0 mg")
axes[0].plot(t, C_20, label="C with D = 20 mg")
axes[0].plot(t, C_50, label="C with D = 50 mg")


########################
# Comparision of 3 doses
########################


# Show Cmax and Tmax
for dose, concentration, color, offset in [
    (20, C_20, "tab:orange", (20, -10)),
    (50, C_50, "tab:green", (20, -10)),
]:
    i_max = np.argmax(concentration)
    t_max = t[i_max]
    c_max = concentration[i_max]

    axes[0].axvline(
        x=t_max,
        color=color,
        linestyle=":",
        alpha=0.7,
    )
    axes[0].scatter(t_max, c_max, color=color, zorder=3)
    axes[0].annotate(
        f"{dose} mg : Cmax = {c_max:.2f} mg/L\nTmax = {t_max:.1f} h",
        xy=(t_max, c_max),
        xytext=offset,
        textcoords="offset points",
        color=color,
    )
    
axes[0].set_ylabel("Concentration C (mg/L)")
axes[0].legend()

# Plot B/B0 curves
axes[1].plot(t, B_0/B0, label="B / B0 with D = 0 mg")
axes[1].plot(t, B_20/B0, label="B / B0 with D = 20 mg")
axes[1].plot(t, B_50/B0, label="B / B0 with D = 50 mg")

# Plot 24h Marker
axes[1].axvline(
    x=24.0,
    color="tab:red",
    linestyle=":",
    alpha=0.7,
)

# Show B/B0 min
for dose, biomarker, color, offset in [
    (20, B_20, "tab:orange", (-50, 20)),
    (50, B_50, "tab:green", (-50, 20)),
]:
    i_min = np.argmin(biomarker)
    t_min = t[i_min]
    b_min = biomarker[i_min]
    
    # Show B/B0 min and T_min
    axes[1].axvline(
        x=t_min,
        color=color,
        linestyle=":",
        alpha=0.7,
    )
    axes[1].scatter(t_min, b_min/B0, color=color, zorder=3)
    axes[1].annotate(
        f"B/B0 min = {b_min/B0:.2f}\nTmin = {t_min:.1f} h",
        xy=(t_min, b_min/B0),
        xytext=offset,
        textcoords="offset points",
        color=color,
    )
    
# Show B/B0 at 24h
for dose, biomarker, color, offset in [
    (20, B_20, "tab:orange", (3, -10)),
    (50, B_50, "tab:green", (3, -10)),
]:
    i_24h = np.argmin(np.abs(t-24))
    b_24h = biomarker[i_24h]
    
    
    axes[1].scatter(24.0, b_24h/B0, color=color, zorder=3)
    axes[1].annotate(
        f"B/B0 at 24h = {b_24h/B0:.2f}",
        xy=(24.0, b_24h/B0),
        xytext=offset,
        textcoords="offset points",
        color=color,
    )
    
axes[1].axhline(0.7, color="gray", linestyle="--", label="30% Reduction")
axes[1].set_ylabel("Biomarker B / B0")
axes[1].set_xlabel("Time (h)")
axes[1].legend()

axes[1].set_xticks([0, 10, 20, 24, 30, 40, 50])
fig.tight_layout()
plt.show()

##########################################
# Exercice : Calibration of parameter IC50
##########################################

#Build a observation times array
t_obs = np.array([0, 2, 4, 8, 12, 18, 24, 36, 48])

#We use the original parameters : ka = 1 h-1, CL = 2 mg/h, V = 20L, Imax = 0.8, IC50 = 0. mg/ml, kprod = 10 UB/h, kdeg = 0.1 h-1, D = 50
indices_obs = np.searchsorted(t, t_obs)
np.testing.assert_allclose(t[indices_obs], t_obs, atol=1e-10)
B_theory = B_50[indices_obs]

#Add noise to B_theory to make an observed B array
rng = np.random.default_rng(69)
noise = rng.normal(
    loc=0.0,                # 0 mean error
    scale=3.0,              # STD = 3 UB
    size=len(B_theory),     # Error on each observation
)
B_obs = B_theory + noise

for time, theory, observe in zip(t_obs, B_theory, B_obs):
    print(f"{time:>4.0f} h : theory={theory:.2f} UB, observation={observe:.2f} UB")
    
# NL regression with least squares - optimization of IC50
fit = least_squares(
    fun=residuals_ic50,
    x0=[0.5],                  # estimation de départ, en mg/L
    bounds=([0.01], [2.0]),   # IC50 reste positif
    args=(params, t_obs, B_obs),
)
if not fit.success:
    raise RuntimeError(fit.message)
    
ic50_estimated = fit.x[0]
print(f"IC50 used to create data : {params['IC50']:.3f} mg/L")
print(f"IC50 estimated : {ic50_estimated:.3f} mg/L")
print(f"Residuals at observed times : {fit.fun}")

fitted_params = params.copy()
fitted_params["IC50"] = ic50_estimated

fitted_result = simulate(params=fitted_params, times=times)
B_fitted = fitted_result.y[2]

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(times, B_fitted, label=f"Adjusted model : IC50 = {ic50_estimated:.3f} mg/L")
ax.scatter(t_obs, B_obs, color="black", label="Synthetic observations", zorder=3)
ax.set_xlabel("Time (h)")
ax.set_ylabel("Biomarker B (UB)")
ax.legend()
fig.tight_layout()
plt.show()


####################################
# Exercise : Sensitivity variability
####################################

nb_individuals = 100

indiv_IC50 = rng.uniform(
    low=0.6 * params["IC50"],
    high=1.4 * params["IC50"],
    size=nb_individuals,
)

indiv_CL = rng.uniform(
    low=0.7 * params["CL"],
    high=1.3 * params["CL"],
    size=nb_individuals,
)

#Simulate individuals with variable IC50 and CL with dose = 20 mg
B_B0_indivs_20 = []
for i in range(nb_individuals):
    params_indiv = set_parameters(params, ["IC50", "CL", "dose"], [indiv_IC50[i], indiv_CL[i], 20.0])
    result = simulate(params=params_indiv, times=times)
    if not result.success:
        raise RuntimeError(f"Simulation of individual {i} : {result.message}")
    B = result.y[2]
    B_B0_indivs_20.append(B/B0)
    
B_B0_indivs_20 = np.stack(B_B0_indivs_20)

#Simulate individuals with variable IC50 and CL with dose = 50 mg
B_B0_indivs_50 = []
for i in range(nb_individuals):
    params_indiv = set_parameters(params, ["IC50", "CL", "dose"], [indiv_IC50[i], indiv_CL[i], 50.0])
    result = simulate(params=params_indiv, times=times)
    if not result.success:
        raise RuntimeError(f"Simulation of individual {i} : {result.message}")
    B = result.y[2]
    B_B0_indivs_50.append(B/B0)

B_B0_indivs_50 = np.stack(B_B0_indivs_50)

B_B0_24h_20 = []

#Get all values of B/B0 at 24h
i_24h = np.argmin(np.abs(times-24))
print(f"Selected time : {times[i_24h]} h")
B_B0_24h_20 = B_B0_indivs_20[:, i_24h]
B_B0_24h_50 = B_B0_indivs_50[:, i_24h]

nb_responding_20 = np.count_nonzero(B_B0_24h_20 <= 0.7)
nb_responding_50 = np.count_nonzero(B_B0_24h_50 <= 0.7)

#Show repartitions in a boxplot

fig, ax = plt.subplots(figsize=(7, 6))

ax.boxplot(
    [B_B0_24h_20, B_B0_24h_50],
    positions=[1, 2],
    showmeans=False,
)

ax.set_xticks(
    [1, 2],
    [
        f"20 mg\n{nb_responding_20}/{len(B_B0_24h_20)} responding",
        f"50 mg\n{nb_responding_50}/{len(B_B0_24h_50)} responding",
    ],
)

ax.axhline(0.7, color="gray", linestyle="--", label="Response threshold : $B(24\,h)/B_0 = 0{,}70$")
ax.set_ylabel("Relative response B/B0 at 24 h")
ax.set_title("Distribution of simulated response at 24h relative to drug dose")
ax.legend()

fig.tight_layout()
plt.show()
