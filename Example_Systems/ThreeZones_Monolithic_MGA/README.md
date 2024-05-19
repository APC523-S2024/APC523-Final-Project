# Small New England: Three Zones, Slack Variables Example

**SmallNewEngland** is set of a simplified versions of the more detailed example system RealSystemExample. It is condensed for easy comprehension and quick testing of different components of the GenX. **SmallNewEngland/ThreeZones_Slack_Variables_Example**, a one-year example with hourly resolution, contains zones representing Massachusetts, Connecticut, and Maine. The ten represented resources include only natural gas, solar PV, wind, and lithium-ion battery storage. It additionally contains example input files establishing slack variables for policy constraints (e.g. the Capacity Reserve Margin, CO2 Cap, etc.). These slack variables allow the relevant constraints to be violated at the cost of a specified objective function penalty, which can be used to either identify problematic constraints without causing infeasibilities in GenX, or to set price caps beyond which policies are no longer enforced. These slack variables will only be created if the relevant input data (Capacity_reserve_margin_slack.csv, CO2_cap_slack.csv, Energy_share_requirement_slack.csv, or the 'PriceCap' column in 'Minimum_capacity_requirement.csv') are present. If any of these inputs are not present, GenX-Benders will isntantiate the relevant policy as a hard constraint, which will throw an infeasibility if violated. We are interested in running the simple baseline (monolithic MGA without Benders).

To run the model, first navigate to the example directory at `GenX-Benders/Example_Systems/ThreeZones_Monolithic_MGA`:

`cd("Example_Systems/ThreeZones_Monolithic_MGA")` in Julia, or

`cd GenX-Benders/Example_Systems/ThreeZones_Monolithic_MGA`

If not running the later .sh on the cluster, you must also export the SLURM_CPUS_PER_TASK variable

`export SLURM_CPUS_PER_TASK=5`
   
Next, ensure that your settings in `GenX_settings.yml` are correct. The default settings use the solver Gurobi (`Solver: Gurobi`), time domain reduced input data (`TimeDomainReduction: 1, MonolithicMGA: 1`). Other optional policies include minimum capacity requirements, a capacity reserve margin, and more. A rate-based carbon cap of 50 gCO<sub>2</sub> per kWh is specified in the `CO2_cap.csv` input file.

Once the settings are confirmed, run the model with the `jobscript_distributed.sh` script in the example directory:

`sh jobscript_distributed.sh`

Once the model has completed, results will write to the `Results` directory. Compare results using

`sh Opts_test.sh`