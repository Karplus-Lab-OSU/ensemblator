from setuptools import setup
from setuptools import find_packages

setup(name='ensemblator',
      description='Ensemblator v3',
      long_description='A tool for analyzing protein structure ensembles.',
      version='3.0.1',
      url='https://github.com/Karplus-Lab-OSU/ensemblator',
      license='CCA 4.0 IPL',
      packages=find_packages('src'),
      package_dir={'': 'src'},
      zip_safe=False
      )
