from biomass import create_model, optimize, run_simulation
from biomass.estimation import InitialPopulation

# create model objects from the original and reconstructed models
original = create_model("biomass_models/original")
reconst = create_model("biomass_models/reconstructed")

# optimize parameters for each model; this may take several hours~couple of days
print("Optimizing parameters for the original model...")
for x_id in range(1, 11):
    initpop = InitialPopulation(original, popsize=15).generate(n_proc=8, progress=False)
    optimize(
        original,
        x_id=x_id,
        disp_here=True if x_id == 1 else False,
        optimizer_options={
            "maxiter": 100,
            "popsize": 15,
            "init": initpop,
            "workers": 8,
        },
    )

print("Optimizing parameters for the reconstructed model...")
for x_id in range(1, 11):
    initpop = InitialPopulation(reconst, popsize=15).generate(n_proc=8, progress=False)
    optimize(
        reconst,
        x_id=x_id,
        disp_here=True if x_id == 1 else False,
        optimizer_options={
            "maxiter": 100,
            "popsize": 15,
            "init": initpop,
            "workers": 8,
        },
    )

# Visualize simulation results
run_simulation(original, viz_type="average", show_all=False, stdev=True)
run_simulation(reconst, viz_type="average", show_all=False, stdev=True)
