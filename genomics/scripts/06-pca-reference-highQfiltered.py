import ipyrad.analysis as ipa
import pandas as pd
import toyplot
import toyplot.pdf

### scripting and comments based on:
### https://github.com/a-lud/sea-snake-dart/blob/main/scripts/06-pca.ipynb

## drop
problems = pd.read_csv(
    "/hpcfs/users/a1235304/atm/data/popmaps/dropped_samples-highQfiltered.tsv",
    sep = "\t",
    comment="#",
    names = ["sample", "population"]
)
problems = dict(problems.groupby("population")["sample"].apply(list))

# Accessory function for dropping samples
def drop_samples(probs: dict, pop_dict: dict) -> dict:
    for k,v in probs.items():
        if k in pop_dict.keys():
            for sample in v:
                if sample in pop_dict[k]:
                    print(f"Dropping: {sample}")
                    pop_dict[k].remove(sample)
    return pop_dict

# Aipysurus foliosquama

## Filtered HDF5 file
afo_50kb_LD = "/hpcfs/users/a1235304/atm/results/ipyrad/AFO-reference_outfiles/AFO-reference.highQ.filtered.LD50k.snps.hdf5"

## Read the A. foliosquama population file as a Pandas dataframe
populations = pd.read_csv(
    "/hpcfs/users/a1235304/atm/data/popmaps/AFO-popmap.tsv",
    comment = "#",
    sep=" ",
    names=["sample", "grouping"]
)

## Convert the Pandas data frame to a dictionary
imap = dict(populations.groupby("grouping")["sample"].apply(list))

## 50% of samples with data in each group
minmap = {i: 0.5 for i in imap}

## Generate the `pca` object
## This object takes all the data we've specified above and runs the PCA analysis.
## `minmap` variable is set at 50% threshold so that the imputation software can do its thing.
pca = ipa.pca(
    data = afo_50kb_LD,
    imap = imap,
    minmap = minmap,
    mincov = 0.75,
    impute_method = "sample"
)
pca.run()

## Save the principal component values as a dataframe so they can be plotted in R, and 
## also save the variance components as a single column table.

## store the PC axes as a dataframe
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df_variance = pd.DataFrame(pca.variances[0])

## write the PC axes to a CSV file
df.to_csv("/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AFO-reference-highQfiltered_pca.csv")
df_variance.to_csv("/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AFO-reference-highQfiltered_pca-variance.csv")

## Generate the figure using `pca.draw()`
figure, _ = pca.draw(
    0, 1,
    width=600, height=500,
    label = "Aipysurus foliosquama (R,Q,filt): PC1 vs PC2"
)

## set bg colour to white
figure.style.update({"background-color": "white"})

## Increase the figure width to prevent legend from being cut off
figure.width = 650

## save as PDF
toyplot.pdf.render(figure, "/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AFO-reference-highQfiltered_PCA-1_2.pdf")

# Drop AFO

## Filtered HDF5 file
afo_50kb_LD = "/hpcfs/users/a1235304/atm/results/ipyrad/AFO-reference_outfiles/AFO-reference.highQ.filtered.LD50k.snps.hdf5"

## Read the A. foliosquama population file as a Pandas dataframe
populations = pd.read_csv(
    "/hpcfs/users/a1235304/atm/data/popmaps/AFO-popmap.tsv",
    comment = "#",
    sep=" ",
    names=["sample", "grouping"]
)

## Convert the Pandas data frame to a dictionary
imap = dict(populations.groupby("grouping")["sample"].apply(list))

# Remove samples based on problems dict
imap = drop_samples(probs=problems, pop_dict=imap)

## 50% of samples with data in each group
minmap = {i: 0.5 for i in imap}

## Generate the `pca` object
## This object takes all the data we've specified above and runs the PCA analysis.
## `minmap` variable is set at 50% threshold so that the imputation software can do its thing.
pca = ipa.pca(
    data = afo_50kb_LD,
    imap = imap,
    minmap = minmap,
    mincov = 0.75,
    impute_method = "sample"
)
pca.run()

## Save the principal component values as a dataframe so they can be plotted in R, and 
## also save the variance components as a single column table.

## store the PC axes as a dataframe
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df_variance = pd.DataFrame(pca.variances[0])

## write the PC axes to a CSV file
df.to_csv("/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AFO-reference-highQfiltered-drop_pca.csv")
df_variance.to_csv("/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AFO-reference-highQfiltered-drop_pca-variance.csv")

## Generate the figure using `pca.draw()`
figure, _ = pca.draw(
    0, 1,
    width=600, height=500,
    label = "Aipysurus foliosquama (R,Q,filt,drop): PC1 vs PC2"
)

## set bg colour to white
figure.style.update({"background-color": "white"})

## Increase the figure width to prevent legend from being cut off
figure.width = 650

## save as PDF
toyplot.pdf.render(figure, "/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AFO-reference-highQfiltered-drop_PCA-1_2.pdf")

# Aipysurus apraefrontalis

## Filtered HDF5 file
aap_50kb_LD = "/hpcfs/users/a1235304/atm/results/ipyrad/AAP-reference_outfiles/AAP-reference.highQ.filtered.LD50k.snps.hdf5"

## Read the A. foliosquama population file as a Pandas dataframe
populations = pd.read_csv(
    "/hpcfs/users/a1235304/atm/data/popmaps/AAP-popmap.tsv",
    comment = "#",
    sep=" ",
    names=["sample", "grouping"]
)

## Convert the Pandas data frame to a dictionary
imap = dict(populations.groupby("grouping")["sample"].apply(list))

## 50% of samples with data in each group
minmap = {i: 0.5 for i in imap}

## Generate the `pca` object
## This object takes all the data we've specified above and runs the PCA analysis.
## `minmap` variable is set at 50% threshold so that the imputation software can do its thing.
pca = ipa.pca(
    data = aap_50kb_LD,
    imap = imap,
    minmap = minmap,
    mincov = 0.75,
    impute_method = "sample"
)
pca.run()

## Save the principal component values as a dataframe so they can be plotted in R, and 
## also save the variance components as a single column table.

## store the PC axes as a dataframe
df = pd.DataFrame(pca.pcaxes[0], index=pca.names)
df_variance = pd.DataFrame(pca.variances[0])

## write the PC axes to a CSV file
df.to_csv("/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AAP-reference-highQfiltered_pca.csv")
df_variance.to_csv("/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AAP-reference-highQfiltered_pca-variance.csv")

## Generate the figure using `pca.draw()`
figure, _ = pca.draw(
    0, 1,
    width=600, height=500,
    label = "Aipysurus apraefrontalis (R,Q,filt): PC1 vs PC2"
)

## set bg colour to white
figure.style.update({"background-color": "white"})

## Increase the figure width to prevent legend from being cut off
figure.width = 650

## save as PDF
toyplot.pdf.render(figure, "/hpcfs/users/a1235304/atm/results/ipyrad/population-structure/pca/AAP-reference-highQfiltered_PCA-1_2.pdf")