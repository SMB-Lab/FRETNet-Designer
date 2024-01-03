Setup instructions:

1. Create a conda environment with python 3.8 within which to run the notebook.
    $ conda create --name myenv_3_8FT python=3.8.0
    
2. Activate the environment.
    conda activate myenv_3_8_FT or source activate myenv_3_8_FT if on Palmetto

3. Install jupyter in the new environment.
    $ conda install jupyter

4. Creating a new kernel associated with the new environment which will allow you run the notebook in the env.
    $ python -m ipykernel install --user --name python_custom --display-name "<name of your env>"

5. Install the required python packages in the environment.
    $ python -m pip install numpy==1.24.4
    $ python -m pip install pandas==2.0.3
    $ python -m pip install matplotlib==3.7.4
    $ python -m pip install mdtraj==1.9.9
    $ python -m pip install fretraj==0.2.10
    $ python -m pip install json==2.0.9




# when Frank tried the following instructions, they did not work. matplotlib was installed with an outdated version which used np.float.
installed in this order:
conda install numpy
conda install pandas
conda install matplotlib
conda install fretraj -c conda-forge        (this installs mdtraj as well; if not, can pre-install it)
if needed:
python -m pip install copy, json, 
this should allow running the notebook
