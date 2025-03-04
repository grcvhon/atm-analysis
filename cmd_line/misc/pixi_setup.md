## Setting up `pixi` package manager
This document describes steps for setting up `pixi` in your Phoenix HPC account (instructions from Alastair Ludington)<br>
<br>

From https://github.com/prefix-dev/pixi:<br>
`pixi` is a cross-platform, multi-language package manager and workflow tool built on the foundation of the conda ecosystem. 
It provides developers with an exceptional experience similar to popular package managers like `cargo` or `yarn`, but for any language.<br>
<br>

#### 1) Install `pixi` as per website
```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

#### 2) In your home directory, make the config directory
```bash
mkdir .config/pixi/
```

#### 3) Create the `config.toml` file
```bash
nano .config/pixi/config.toml
```

#### 4) Copy the lines in config.toml 
```nginx
# Channels
default_channels = ["conda-forge", "bioconda"]
detached-environments = "/home/a#######/hpcfs/.pixi/envs" # change a-number to user
```

#### 5 Append lines below to `.bashrc` file `(cd to ~/.bashrc)`
```nginx
# pixi
export PATH="/home/a1235304/.pixi/bin:/hpcfs/users/$USER/.pixi/bin:$PATH"
export PIXI_CACHE_DIR="/hpcfs/users/$USER/.cache/rattler/cache"
export PIXI_HOME="/hpcfs/users/$USER/.pixi"
```
To apply changes, run `source .bashrc` or restart terminal.
