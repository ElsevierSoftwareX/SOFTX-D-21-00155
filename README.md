**sas_temper**

SAS data analysis package that uses simulated annealing to fit data using the sasmodels package.
The program automatically performs multiple runs of the fitting to check the reproducibility of the result.

The program is written in python3.

The code began as sa_fitter, which was developed as part of ORNL LDRD project 8235.

**System Requirements**
If your system does not presently have python 3 and the required libraries,
you can create a specific environment using conda and the sas_temper_env.yaml file here.

```
conda env create -f sas_temper_env.yaml
```

To activate the environment so you can run sas_temper, type the following.

```
conda activate sas_temper_env
```

**Installation**
pip install --upgrade git+https://code.ornl.gov/wt3/sas_temper.git

or, if you want a specific branch

pip install --upgrade git+https://code.ornl.gov/wt3/sas_temper.git@<branch_name>

**Usage**
An example is available in the `examples` directory:

```
cd examples
sas_temper test_sphere_nz.yaml
```
