import ipyrad.analysis as ipa

# Aipysurus foliosquama
print("Starting conversion of vcf to hdf5 for Aipysurus foliosquama...")
print(" ")

print("### CONVERTING: AFO - reference - Unfiltered ipyrad output to hdf5 ###")
## AFO - reference - Unfiltered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AFO-reference.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AFO-reference_outfiles/AFO-reference.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AFO-reference_outfiles/',
    ld_block_size=50000
)
converter.run()
print("### COMPLETE: AFO - reference - Unfiltered ###")

print(" ")

print("### CONVERTING: AFO - reference - Filtered ipyrad output to hdf5 ###")
## AFO - reference - Filtered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AFO-reference.highQ.filtered.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AFO-reference_outfiles/AFO-reference.highQ.filtered.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AFO-reference_outfiles/',
    ld_block_size=50000
)
converter.run()
print("### COMPLETE: AFO - reference - Filtered ###")

print(" ")

print("### CONVERTING: AFO - de novo - Unfiltered ipyrad output to hdf5 ###")
## AFO - de novo - Unfiltered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AFO-denovo.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AFO-denovo_outfiles/AFO-denovo.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AFO-denovo_outfiles/',
    ld_block_size=50000
)
converter.run()
print("### COMPLETE: AFO - de novo - Unfiltered ###")

print(" ")

print("### CONVERTING: AFO - de novo - Filtered ipyrad output to hdf5 ###")
## AFO - de novo - Filtered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AFO-denovo.highQ.filtered.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AFO-denovo_outfiles/AFO-denovo.highQ.filtered.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AFO-denovo_outfiles/',
    ld_block_size=50000
)
converter.run()
print("### COMPLETE: AFO - de novo - Filtered ###")

print(" ")
print("Aipysurus foliosquama: All conversion complete.")

print(" ")
# Aipysurus apraefrontalis
print("Starting conversion of vcf to hdf5 for Aipysurus apraefrontalis...")
print(" ")

print("### CONVERTING: AAP - reference - Unfiltered ipyrad output to hdf5 ###")
## AAP - reference - Unfiltered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AAP-reference.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AAP-reference_outfiles/AAP-reference.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AAP-reference_outfiles/',
    ld_block_size=50000
)
converter.run()
print("### COMPLETE: AAP - reference - Unfiltered ###")

print(" ")

print("### CONVERTING: AAP - reference - Filtered ipyrad output to hdf5 ###")
## AAP - reference - Filtered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AAP-reference.highQ.filtered.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AAP-reference_outfiles/AAP-reference.highQ.filtered.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AAP-reference_outfiles/',
    ld_block_size=50000
)
converter.run()
print("### COMPLETE: AAP - reference - Filtered ###")

print(" ")

print("### CONVERTING: AAP - de novo - Unfiltered ipyrad output to hdf5 ###")
## AAP - de novo - Unfiltered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AAP-denovo.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AAP-denovo_outfiles/AAP-denovo.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AAP-denovo_outfiles/',
    ld_block_size=50000
)
converter.run()
print("### COMPLETE: AAP - de novo - Unfiltered ###")

print(" ")

print("### CONVERTING: AAP - de novo - Filtered ipyrad output to hdf5 ###")
## AAP - de novo - Filtered ipyrad output
converter = ipa.vcf_to_hdf5(
    name="AAP-denovo.highQ.filtered.LD50k",
    data="/hpcfs/users/a1235304/atm/results/ipyrad/AAP-denovo_outfiles/AAP-denovo.highQ.filtered.vcf.gz",
    workdir='/hpcfs/users/a1235304/atm/results/ipyrad/AAP-denovo_outfiles/',
    ld_block_size=50000
)
converter.run()
print("### COMPLETE: AAP - de novo - Filtered ###")

print(" ")
print("Aipysurus apraefrontalis: All conversion complete.")