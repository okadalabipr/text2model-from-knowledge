from biomass import create_model, optimize
from modified_biomass import run_simulation
from biomass.estimation import InitialPopulation


# show simulation results for original model;
# parameter optimization was not performed for the original model
model = create_model("Raia2011_original")
run_simulation(model, viz_type="original", show_all=False, stdev=True)


for model_name in ["KEGG_JAK-STAT", "Raia2011_reconst"]:
    print(f"Creating {model_name}...")
    model = create_model(model_name)

    # optimize parameters for each model; this may take several hours~couple of days
    print(f"Optimizing parameters for {model_name}...")
    for x_id in range(1, 11):
        initpop = InitialPopulation(model, popsize=15).generate(n_proc=8, progress=False)
        optimize(
            model,
            x_id=x_id,
            disp_here=True if x_id == 1 else False,
            optimizer_options={
                "maxiter": 100,
                "popsize": 15,
                "init": initpop,
                "workers": 8,
            },
        )

    # visualize simulation results
    print(f"Visualizing simulation results for {model_name}...")
    run_simulation(model, viz_type="average", show_all=False, stdev=True)
