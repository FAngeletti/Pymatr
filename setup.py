#! /usr/bin/env python
from distutils.core import setup

setup(
	name="Pymatr",
	version="1.0",
	description='Python matrix-correlated random variable library',
	author='Florian Angeletti',
	author_email='angeletti@achronie.fr',
	url='https://github.com/FAngeletti/Pymatr.git',
	packages=['Pymatr'],
	package_dir = {'Pymatr': 'src'},
	package_data={'Pymatr' : ["src/tests/*.py"]},
	requires=['numpy','sympy','matplotlib'],
	provides="pymatr"

)
