from setuptools import setup,find_packages
setup(name='BF2PALM',
      version='0.0.37',
      description='Transfer building footprint to 2D DEM for LES simulation',
      author='Jiachen Lu',
      author_email='jiachensc@gmail.com',
      requires= ['numpy','matplotlib','shapely','pandas','scipy','osmnx'], 
      install_requires= ['numpy','matplotlib','shapely','pandas','scipy','osmnx','netCDF4'], 
      packages=["BF2PALM"],
      license="MIT",
      python_requires=">=3.6",
      )
