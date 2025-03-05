## Setting up `pixi` package manager
This document describes steps for setting up `pixi` in your Phoenix HPC account (instructions from Dr Alastair Ludington)<br>
<br>

`pixi` is a cross-platform, multi-language package manager and workflow tool built on the foundation of the conda ecosystem. 
It provides developers with an exceptional experience similar to popular package managers like `cargo` or `yarn`, but for any language. (https://github.com/prefix-dev/pixi)<br>
<br>

The main benefit of using `pixi` is the freedom to call commands from conda packages (e.g. `bcftools`,`plink2`,`fastp`) without having to run `conda activate <environment/software>` every time you need it. The conda packages will be installed globally (see below) so that software can be invoked directly (i.e., as if they are already loaded).<br>
<br>

---

#### Set-up

##### 1) Install `pixi` as per website
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

##### 2) In your home directory, make the config directory
```bash
mkdir .config/pixi/
```

##### 3) Create the `config.toml` file
```bash
nano .config/pixi/config.toml
```

##### 4) Copy the lines in config.toml 
```nginx
# Channels
default_channels = ["conda-forge", "bioconda"]
detached-environments = "/home/a#######/hpcfs/.pixi/envs" # change a-number to user
```

##### 5 Append lines below to `.bashrc` file `(cd to ~/.bashrc)`
```nginx
# pixi
export PATH="/home/a1235304/.pixi/bin:/hpcfs/users/$USER/.pixi/bin:$PATH"
export PIXI_CACHE_DIR="/hpcfs/users/$USER/.cache/rattler/cache"
export PIXI_HOME="/hpcfs/users/$USER/.pixi"
```
To apply changes, run `source .bashrc` or restart terminal.

---

#### Install a conda package globally
Here is an example on how to install a conda package globally with pixi.
```bash
pixi global install fastp
```
This command should install `fastp` software. Check that `fastp` was installed successfully and can be called directly by running `fastp`.
