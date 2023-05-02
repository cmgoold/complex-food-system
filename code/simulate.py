#!/usr/bin/env python3 

import numpy as np
import cmdstanpy as csp
import matplotlib.pyplot as plt
from scipy.stats import lognorm
from scipy.stats import gaussian_kde
import os
import dill

from ode.system import FoodSystem, DimensionlessFoodSystem
from ode.solve import Rk4

seed = 1234
rng = {"random_state": np.random.default_rng(seed)}
root = os.getcwd() + "/"
model_file = root + "stan/stan-sim-model.stan"
model = csp.CmdStanModel(stan_file=model_file)

parameters = {
    "a": 0.05,
    "b": 130,
    "e": 1/(2.5*12),
    "f": 24/12,
    "g": 110*0.75,
    "k": 0.3,
    "h": 130e6,
    "w": 0.3,
    "m": 1/6,
    "q": 160,
    "r": 1/6,
    "s": 1,
}

states = {
    "C": 400e3,
    "I": 130e6,
    "D": 130e6,
    "P": 160,
}

dt = 0.1
N = 120
times = np.arange(0, 120, dt)
system = FoodSystem(parameters=parameters, state_names=list(states.keys()))
dsystem = DimensionlessFoodSystem.from_dimensional(system, states["C"])
rk4 = Rk4(system, states, times)
rk4.solve(cache=True)

# Sensitivities
low, high = (0.1, 1.9)
n = 20
key_pars = ("a", "b", "e", "k", "q", "s", "w")
crs = {}
for par in key_pars:
    for i in np.linspace(parameters[par]*low, parameters[par]*high, n):
        sys = FoodSystem({**parameters, **{par: i}}, state_names=list(states.keys()))
        dsys = DimensionlessFoodSystem.from_dimensional(sys, states["C"])
        if par in crs:
            crs[par].append(dsys.critical_ratio)
        else:
            crs[par] = [dsys.critical_ratio]

def _add_noise(x, sigma):
    return lognorm(scale=x, s=sigma).rvs(**rng)

sigma = 0.1
fit_times = [t for t in times if t % 1 == 0]
noisy_states = [[_add_noise(i, sigma) for _, i in v.items()] for k, v in rk4.cache.items() if k % 1 == 0]
production = [v[0] * system.g * system.f for v in noisy_states]
imports = [_add_noise(system.k * system.h, sigma) for _ in production]
exports = [v[0] * system.k * system.g * system.f for v in noisy_states]

stan_data = {
    "n_t": len(fit_times),
    "ts": list(range(1, len(noisy_states) + 1)),
    "n_states": len(states),
    "n_parameters": len(parameters),
    "y": noisy_states,
    "production": production,
    "imports": imports,
    "exports": exports,
    "prior_only": 0
}

fit_file = root.replace("code/", "results/sim/fit.pkl")

if not os.path.exists(fit_file):
    fit = model.sample(data=stan_data, parallel_chains=4, inits=0, seed=seed)
    print(fit.diagnose())

    with open(fit_file, "wb") as f:
        dill.dump(fit, f)

else:
    fit = dill.load(open(fit_file, "rb"))

draws = fit.stan_variables()

## Make figure
fig_path = root.replace("code/", "results/plots/")

# posterior predictive
n = 50
indices = rng["random_state"].choice(range(len(fit_times)), n, replace=False)
state_sigmas = draws["sigma"][indices,:]
state_paths = draws["states"][indices, :, :]
mean_paths = state_paths.mean(axis=0)
pp = np.array([
    [
        _add_noise(path, state_sigmas[path_idx, i]) 
        for path_idx, path
        in enumerate(s.T)
    ]
    for i, s 
    in 
    enumerate(state_paths.T)
]).T

###### Figure 2
figure2 = plt.figure(figsize=(20, 12), layout="constrained")
plt.style.use('classic')
gs = figure2.add_gridspec(4, 16)

# Sensitivities
sens = figure2.add_subplot(gs[:2, :8])
sens.scatter(
    [i for i, _ in enumerate(crs)],
    [np.mean(np.diff(v)) for v in crs.values()],
    c="#1da1f2", s=100
)
sens.axhline(y=0, ls="dashed", c="k", lw=2)
sens.set_ylabel("Critical ratio gradient")
sens.grid()
plt.xticks(
    [i for i, _ in enumerate(crs)],
    labels=[k for k in crs.keys()]
)

# Stability regions
stability_axes = [
    figure2.add_subplot(gs[0, 8:12]),
    figure2.add_subplot(gs[0, 12:]),
    figure2.add_subplot(gs[1, 8:12]),
    figure2.add_subplot(gs[1, 12:]),
]

alphas = np.linspace(0.01, 2, 100)
kappas = np.linspace(0.01, 0.99, 100)
betas = [0.1, 0.5, 1, 1.5]
acs = {}
for beta in betas:
    for kappa in kappas:
        dsys = DimensionlessFoodSystem.from_dimensional(system, states["C"])
        dsys.parameters["beta"] = beta
        dsys.parameters["kappa"] = kappa
        acs[(beta, kappa)] = dsys

zeros = np.zeros(len(kappas))
threes = [3.0]*len(kappas)
for i, (ax, beta) in enumerate(zip(stability_axes, betas)):
    ac = {k: v for k, v in acs.items() if k[0] == beta}
    ax.fill_between([kappa for _, kappa in ac.keys()], zeros, [v.critical_alpha for v in ac.values()],
        color="red",
        label="Unsustainable",
    )
    ax.fill_between(
        [kappa for _, kappa in ac.keys()],
        [v.critical_alpha for v in ac.values()],
        [2 * v.gamma * (1 + v.beta)/(v.gamma + 2 * v.omega) for v in ac.values()],
        color="#1da1f2",
        label="Sustainable (import reliant)",
    )
    ax.fill_between(
        [kappa for _, kappa in ac.keys()],
        [2 * v.gamma * (1 + v.beta)/(v.gamma + 2 * v.omega) for v in ac.values()],
        [3.0]*len(ac),
        color="green",
        label="Sustainable (over-producing)",
    )
    ax.set_ylim((0, 3))
    ax.set_xlim((0.01, 0.99))
    if i == 0 or i == 2:
        ax.set_ylabel(r"$\alpha$")

    if i == 0:
        ax.text(0.6, 0.5, "Unsustainable", c="white")
        ax.text(0.1, 1.0, "Sustainable (importing)")
        ax.text(0.1, 2.0, "Sustainable (over-producing)")
    ax.title.set_text(rf"$\beta$ = {beta}")
    ax.set_xlabel(r"$\kappa$")

# Path simulations
sim_axes = [
    figure2.add_subplot(gs[2, :4]),
    figure2.add_subplot(gs[2, 4:8]),
    figure2.add_subplot(gs[2, 8:12]),
    figure2.add_subplot(gs[2, 12:]),
]
titles = {"C": rf"Captial ($v$)", "I": "Inventory ($x$)", "D": "Demand ($y$)", "P": "Price ($z$)"}
for i, ax in enumerate(sim_axes):
    ax.grid()
    k = list(states.keys())[i]
    for idx, path in enumerate(pp[:, :, i].T):
        ax.scatter(
            fit_times, 
            path/states[k], 
            edgecolors="#1da1f2", 
            facecolors="none",
            alpha=0.5, 
            label=None if idx else "Posterior predictive"
        )

    ax.scatter(
        fit_times, 
        [noisy_states[ii][i]/states[k] for ii, _ in enumerate(fit_times)], 
        facecolor="indianred", 
        edgecolor="indianred",
        alpha=0.6, 
        label="Simulated + noise"
    )
    for idx, path in enumerate(state_paths[:, :, i]):
        ax.plot(fit_times, path/states[k], c="yellow", lw=0.1, label=None if idx else "Mean posterior samples")
    ax.plot(fit_times, [rk4.cache[t][k]/states[k] for t in fit_times], c="k", lw=2, label="Simulated trajectory")
    ax.title.set_text(titles[k])
    ax.set_xlabel("Simulation time")
    ax.set_xlim((0, 120))
    if i == 0:
        ax.set_ylabel("State variable")
    if i == 1:
        legend = ax.legend(frameon=True, prop={"size": 10})
        for line in legend.get_lines():
            line.set_linewidth(1.0)

# Critical ratios
ratio_axes = [figure2.add_subplot(gs[3, :8]), figure2.add_subplot(gs[3, 8:])]
critical_ratio = draws["critical_ratio"]
surplus_ratio = draws["surplus_ratio"]

ratio_axes[0].grid()
ratio_axes[0].fill_between(
    sorted(critical_ratio), 
    [0]*len(critical_ratio),
    gaussian_kde(critical_ratio)(sorted(critical_ratio)), 
    color="#1da1f2", alpha=0.8,
)
ratio_axes[0].axvline(critical_ratio.mean(), ls="dashed", c="indianred", lw=3, label="Simulated value")
ratio_axes[0].title.set_text("Critical ratio")
ratio_axes[0].set_ylabel("Density")
ratio_axes[0].set_xlabel("Posterior estimate")
ratio_axes[0].legend(frameon=False)
ratio_axes[1].grid()
ratio_axes[1].fill_between(
    sorted(surplus_ratio), 
    [0]*len(surplus_ratio),
    gaussian_kde(surplus_ratio)(sorted(surplus_ratio)), 
    color="#1da1f2", alpha=0.8,
)
ratio_axes[1].axvline(surplus_ratio.mean(), ls="dashed", c="indianred", lw=3)
ratio_axes[1].title.set_text("Surplus ratio")
ratio_axes[1].set_xlabel("Posterior estimate")

# labels
plt.gcf().text(0, 1, "A", fontsize=20)

plt.savefig(fig_path + "figure-2.png", dpi=200)
