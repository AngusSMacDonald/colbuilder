# Install UCSF Chimera

Install UCSF Chimera (classic, NOT ChimeraX). Make sure the executable is executable and added to your `$PATH` (e.g. via `chmod +x chimera...`, `./chimera...`).
Add the following to the end of `~/.bashrc`: `export PATH="/home/user/Colbuilder/chimera/bin:$PATH"`


# Install Modeller

Download Modeller 10.5:

`wget https://salilab.org/modeller/10.5/modeller-10.5.tar.gz`

Alternatively, download from:
https://salilab.org/modeller/old_releases.html
(select "Other Linux/Unix..." download)

Note: you must obtain an academic license from the Modeller website (free).

Untar and install:

`tar -zxvf modeller-10.5.tar.gz`
`cd modeller-10.5`
`./Install`

During installation:
- Choose option #2 for x86_64
- Note the install path (e.g. /home/user/bin/modeller10.5/)

Append the following lines to the end of `~/.bashrc` (adjust paths and user):
```
export PYTHONPATH="/home/user/bin/modeller10.5/lib/x86_64-intel8/python3.3:$PYTHONPATH"
export PYTHONPATH="/home/user/bin/modeller10.5/modlib:$PYTHONPATH"
export LD_LIBRARY_PATH="/home/user/bin/modeller10.5/lib/x86_64-intel8:$LD_LIBRARY_PATH"
export PATH="$HOME/.local/bin:$PATH"
```

# Create conda environment and install ColBuilder
*In this case an editable environment for development. If you want a fixed local install, do `pip install .` instead, and change environment name to reflect.*
```
conda create -n colbuilder-dev python=3.9
conda activate colbuilder-dev

git clone https://github.com/AngusSMacDonald/colbuilder.git
cd colbuilder
pip install -e .
```

# Install dependencies
```
conda install conda-forge::pymol-open-source
conda install bioconda::muscle
```

# Verify installation

`colbuilder --help`

# Running ColBuilder

To run colbuilder commands, you must currently be in the directory where the
repository was cloned (i.e. the `colbuilder` directory) as this is an editable (live) install. You could instead do a permanent install (`pip install .`) at the start, and add this to PATH.

# Test run

Find a sample test configuration from the ColBuilder GitHub instructions,
save it as `test.yaml`, and run:

`colbuilder --config_file test.yaml`

