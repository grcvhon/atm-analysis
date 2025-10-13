## Monitoring of the short-nosed and Leaf-scaled sea Snakes to track trajectory

This repository contains input data and code used for developing species distribution models (SDM) for the project: <i><b>Monitoring of the Short-nosed and Leaf-scaled sea snakes to track trajectory</i></b>.

Overall, this project aims to 1) determine baseline populations and develop monitoring methods through a mark-recapture program; and 2) generate more robust species distribution models for two Critically Endangered sea snake species (EPCB Act 1999): the Short-nosed sea snake (<i>Aipysurus apraefrontalis</i>) and the Leaf-scaled sea snake (<i>A. foliosquama</i>)

We look to incorporate genetic information in generating species distribution models. Doing so may progress our understanding of current species-specific distributional ranges as well as identify potential key locations for future surveys.

Our SDM approach is summarised in the diagram below:


```mermaid
flowchart TD;
    A["Occurrence data"]
    B["Habitat/environment layers"]
    C["Genetic layer"]
    D["Passage probability layer"]
    E["SDM analysis"]
    F["MaxEnt"]
    G["Random Forest"]
    H["General Additive 
    Mixed Models
    (GAMM)"]
    I["Ensemble"]
    
    style A fill:white,stroke:grey,color:black,font-size:15px
    style B fill:white,stroke:grey,color:black,font-size:15px
    style C fill:white,stroke:grey,color:black,font-size:15px
    style D fill:white,stroke:grey,color:black,font-size:15px
    style E fill:white,stroke:grey,color:black,font-size:15px
    style F fill:white,stroke:grey,color:black,font-size:15px
    style G fill:white,stroke:grey,color:black,font-size:15px
    style H fill:white,stroke:grey,color:black,font-size:15px
    style I fill:white,stroke:grey,color:black,font-size:15px

    A--
    <p style="font-size:14px">
    <b>Species-specific data</b><br> 
    (<i>Aipysurus laevis</i>,<br>
    <i>A. apraefrontalis</i>,<br>
    <i>A. foliosquama</i>)
    </p>--->E;

    B--
    <p style="font-size:14px">
    Sources:<br>
    UNEP, IMOS,<br> 
    CSIRO, BioOracle
    </p>--->E;
    
    C--
    <p style="font-size:14px">
    Interpolated ancestry<br>
    coefficient values<br>
    for <i>A. laevis</i><br>
    (<i>algatr</i>, R)
    </p>--->E;
    
    D--
    <p style="font-size:14px">
    Mean passage probability<br>
    values based on ocean<br>
    current direction<br>
    (<i>gdistance</i>, R)
    </p>--->E;
    
    E-->F;
    E-->G;
    E-->H;
    E-->I;
```
##

#### List of contents
* [Generate spatial genetic layer for species distribution modelling](https://github.com/grcvhon/atm-analysis/tree/master/genetic_layer)
* [Generate spatial passage layer for species distribution modelling](https://github.com/grcvhon/atm-analysis/tree/master/passage_layer)

##
<sub>
MERIT Project ID: SNS-TST-P2 | Programme: Saving Native Species | Sub-programme: Tracking Species Trajectories | Commonwealth Contract Services Reference: ATM_2023_0715<br>
This project is funded by the Australian Government’s Saving Native Species Program and delivered by the University of Adelaide.
</sub>
