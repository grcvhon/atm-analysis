# Omniscape X Circuitscape
# Created 12 June 2025
# Vhon Garcia

# This is also my first time using Julia lang so I have notes below
# about basic commands
# - pwd() print present working dir
# - cd("path/to/dir") change working dir
# - plot() ... capable of plotting, will open new window to print plot

#################
# Omniscape
#################
# Following https://docs.circuitscape.org/Omniscape.jl/latest/usage/
# Circuitscape.jl installation
using Pkg
Pkg.add("Omniscape")

#################
# Circuitscape
#################
# Following https://github.com/Circuitscape/Circuitscape.jl

# Circuitscape.jl installation
using Pkg
Pkg.add("Circuitscape")

# Check if all the tests are passing by doing the following:
Pkg.test("Circuitscape")
