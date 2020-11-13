from setuptools import setup

REQUIREMENTS_DOT_TXT="requirements.txt"

def required_packages():
    return [format_dependency(dep) for dep in requirements_txt()]

def format_dependency(dependency):
    package = dependency.split('#egg=')[-1]
    return package.replace('-','==')

def requirements_txt():
    with open(REQUIREMENTS_DOT_TXT) as f:
        return f.read().splitlines()

def external_dependency_links():
    return filter(is_external_dependency, requirements_txt())

def is_external_dependency(dependency):
    return dependency.startswith('git')

setup(
    name='sas_temper',
    version='0.2.2',
    description='SAS data analysis using simulated annealing and reproducibility characterization.  Uses sasmodels package',
    packages=['sas_temper'],
    scripts=['scripts/sas_temper'],
    instal_requires=required_packages(),
    dependency_link=external_dependency_links(),
    package_data={'': [REQUIREMENTS_DOT_TXT]}
    )
