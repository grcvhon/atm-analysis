import ipyrad.analysis as ipa

# Aipysurus foliosquama

## Unfiltered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AFO-reference.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AFO-reference_outfiles/AFO-reference.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AFO-reference_outfiles/',
    ld_block_size=50000
)

converter.run()