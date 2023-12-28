Create a conda environment within which to run the script. This assumes you have a conda install capable of installing Py3.8.0 (what was used here)
The following will create an environment and allow it to be displayed in the jupyterlab kernel selection window
conda create --name myenv_3_8FT python=3.8.0
conda activate myenv_3_8_FT or source activate myenv_3_8_FT if on Palmetto
conda install jupyter
# create a new kernel so that  you can switch between different Jupyter kernels representing different Conda environments without 
# needing to restart the Jupyter Notebook server
python -m ipykernel install --user --name python_custom --display-name "<name of your env>"
packages needed (* need to be installed, as they aren't defaults in conda installs)
numpy*
pandas*
matplotlib*
json       
mdtraj*
fretraj*
copy


installed in this order:
conda install numpy
conda install pandas
conda install matplotlib
conda install fretraj -c conda-forge        (this installs mdtraj as well; if not, can pre-install it)
if needed:
python -m pip install copy, json, 
this should allow running the notebook
