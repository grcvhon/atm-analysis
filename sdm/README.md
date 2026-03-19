# Species distribution modelling
This directory contains code and input data for developing an updated species distribution model (SDM) for the Short-nosed (<i>Aipysurus apraefrontalis</i>) and Leaf-scaled sea snakes (<i>A. foliosquama</i>). In developing these species-specific SDMs, we will incorporate the following layers as predictor variables: (1) environmental and habitat layers, (2) [genetic layer](https://github.com/grcvhon/atm-analysis/tree/master/genetic_layer), (3) [passage probability layer](https://github.com/grcvhon/atm-analysis/tree/master/passage_layer).<br>
<br>
The scripts were divided into three sequential files (so far). Scripts will have steps that will write the output into file to avoid need of running it from the start every time.<br>
<br>
The scripts are as follows:<br>
1) [Generating input for modelling step](https://github.com/grcvhon/atm-analysis/blob/b86915b87133145658a11a1411aab8e990bc08bc/sdm/scripts/01_generate_modelling_input.R)
2) [Modelling step](https://github.com/grcvhon/atm-analysis/blob/b86915b87133145658a11a1411aab8e990bc08bc/sdm/scripts/02_optimise_model.R) (where the optimised model is chosen, variable importance calculated, and response curves plotted)
3) [Prediction step](https://github.com/grcvhon/atm-analysis/blob/b86915b87133145658a11a1411aab8e990bc08bc/sdm/scripts/03_predict.R) (where the optimised model is used to predict species distribution)<br>
<br>

## 1) Generating input for modelling step
This section processes and generates the following input required for the modelling step:

<b><i>Model extent</b></i>: Australia's northwest shelf<br>
<br>
<b><i>Occurrence dataset</b></i><br>

The occurrence data were compiled into a master dataset ([`ATM_master-occurrence-dataset.csv`](https://github.com/grcvhon/atm-analysis/tree/master/sdm/occurrence-data/ATM_master-occurrence-dataset.csv)). These records of Short-nosed and Leaf-scaled sea snakes were sourced from the following:<br>
* `trawled_seasnakes.xlsx`(2009, 2014 - 2021)<br>
* Atlas of Living Australia (downloaded on 21 January 2025)
* DPIRD Fisheries surveys in Exmouth Gulf and Shark Bay, Western Australia (June 2024 - December 2025)
* One spreadsheet (KLS Lab catalogue records)

Occurrence data are point coordinates. When captured via trawling (with start and end coordinates), the mid latitude and mid longitude values were used as proxy for point coordinates; and were calculated as the average latitude and longitude values of the start and end coordinates.<br>
<br>
<b><i>Bias layer</b></i><br>

This layer was generated using sea snake occurrence data compiled in `2024-12-01_SnakeOcc(02).csv`. This layer was written to file.<br>
<br>
<b><i>Background points</b></i><br>

This layer was generated for each species using a combination of:
* 1,000 random background points within the bias layer raster that excludes areas where the species is known to be present
* 1,000 random points within the northwest shelf<br>

The layer was written to file. The seed number was not noted.<br>
<br>
<b><i>Predictor variables</b></i><br>

An initial stack of predictor variables was generated using the following:<br>
* Mean salinity
* Bathymetry
* Mean sea surface temperature
* Distance to land
* Distance to reef
* Distance to freshwater
* Non-passage probability (scaled 0-1)
* Ancestry coefficient (scaled to represent two clusters in same layer, -1 to 0 to 1)
<br>

More on the ancestry coefficient layer: the values were scaled between -1 to 1 to represent the two clusters in the same layer. That is, anything between 0 and 1 belongs to cluster 1 (coastal) and anything between -1 and 0 belongs to cluster 2. By doing this, all of K1 values (0-1, blue in right plot) and all of K2 values (0-1, red in right plot) are represented in the genetic layer for the SDM. It also avoids constraining what’s coastal or offshore based on 0.5 value.

<p align = center>
<img src="https://raw.githubusercontent.com/grcvhon/atm-analysis/master/sdm/plots/clust.png">
</p>

These variables were tested for multicollinearity. The final rasterstack of predictor variables was written to file.<br>
<br>

## 2) Modelling step
Using all the input from step 1, here the optimised model is chosen and output such as variable importance plots and response curves are generated. The optimised model for each species was written into an RDS file so that the optimisation step is only run once.

Variable importance plots and response curves for each species are shown below:<br>
<br>

<div align = "center">
<b>Leaf-scaled sea snake (<i>Aipysurus foliosquama</i>)</b>
</div>
<p align = center>
<img src="https://raw.githubusercontent.com/grcvhon/atm-analysis/master/sdm/plots/leaf_scaled_opt_mod_variable_importance_plot.png", width = 49%, height = 49%>
<img src="https://raw.githubusercontent.com/grcvhon/atm-analysis/master/sdm/plots/leaf_scaled_opt_mod_response_curve_plot.png", width = 49%, height = 49%>
</p>
<br>

<div align = "center">
<b>Short-nosed sea snake (<i>Aipysurus apraefrontalis</i>)</b>
</div>
<p align = center>
<img src="https://raw.githubusercontent.com/grcvhon/atm-analysis/master/sdm/plots/short_nosed_opt_mod_variable_importance_plot.png", width = 49%, height = 49%>
<img src="https://raw.githubusercontent.com/grcvhon/atm-analysis/master/sdm/plots/short_nosed_opt_mod_response_curve_plot.png", width = 49%, height = 49%>
</p>
<br>

## 3) Prediction step
This step is where the optimised model is used to predict species distribution. The ancestry coefficient layer is also modified here in order to predict distributions in coastal and offshore locations.

Recall that the ancestry coefficient layer has a 
