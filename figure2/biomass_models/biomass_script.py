from biomass import create_model, optimize, run_simulation
from biomass.estimation import InitialPopulation

# create model objects from the mcf7 and mdamb231 models
mcf7 = create_model("biomass_models/KEGG_erbb_MCF7")
mda = create_model("biomass_models/KEGG_erbb_MDAMB231")

# optimize parameters for each model
print("Optimizing parameters for MCF7 model")
for x_id in range(1, 11):
    initpop = InitialPopulation(mcf7, popsize=15).generate(n_proc=8, progress=False)
    optimize(
        mcf7,
        x_id=x_id,
        disp_here=True if x_id == 1 else False,
        optimizer_options={
            "maxiter": 100,
            "popsize": 15,
            "init": initpop,
            "workers": 8,
        },
    )

print("Optimizing parameters for MDAMB231 model")
for x_id in range(1, 11):
    initpop = InitialPopulation(mda, popsize=15).generate(n_proc=8, progress=False)
    optimize(
        mda,
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
run_simulation(mcf7, viz_type="average", show_all=False, stdev=True)
run_simulation(mda, viz_type="average", show_all=False, stdev=True)
