from setuptools import setup,find_packages
setup(name='OSM2PALM',
      version='0.0.21',
      description='Transfer building footprint to 2D DEM for LES simulation',
      author='Jiachen Lu',
      author_email='jiachensc@gmail.com',
      requires= ['numpy','matplotlib','shapely','pandas','scipy','osmnx'], 
      install_requires= ['numpy','matplotlib','shapely','pandas','scipy','osmnx','netCDF4'], 
      packages=["OSM2PALM"],
      license="MIT",
      python_requires=">=3.6",
      )
