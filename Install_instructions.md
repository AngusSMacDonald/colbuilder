
#Install UCSF Chimera (NOT ChimeraX) and make sure to add to $PATH (chmod +x chimera..., ./chimera...)
# e.g. in ~/.bashrc, add at end: export PATH="/home/user/Colbuilder/chimera/bin:$PATH"


wget https://salilab.org/modeller/10.5/modeller-10.5.tar.gz
#or go to https://salilab.org/modeller/old_releases.html and grab "Other Linux/Unix..." download
#Need to go to Modeller website and get academic license (is free)

#untar (tar -zxvf) and cd into directory, ./Install to run installer
#Pick option #2 to get x86_64 install
#Note install path e.g.: /home/user/bin/modeller10.5/...
#Open ~/.bashrc file and append to end following lines:

#export PYTHONPATH="/home/user/bin/modeller10.5/lib/x86_64-intel8/python3.3:$PYTHONPATH"
#export PYTHONPATH="/home/user/bin/modeller10.5/modlib:$PYTHONPATH"
#export LD_LIBRARY_PATH="/home/user/bin/modeller10.5/lib/x86_64-intel8:$LD_LIBRARY_PATH"
#export PATH="$HOME/.local/bin:$PATH"

conda create -n colbuilder python=3.9
conda activate colbuilder
git clone https://github.com/graeter-group/colbuilder.git
cd colbuilder
pip install .

conda install conda-forge::pymol-open-source
conda install bioconda::muscle


colbuilder --help

#However, to run colbuilder commands etc you need to be in this colbuilder directory where the git clone installed too
#Maybe you can add to path?

#Also for me pymol didn't want to install until after modeller was installed... so just keep trying


#Find sample test script from github instructions, save as test.yaml, run
colbuilder --config_file test.yaml
