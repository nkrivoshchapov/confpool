from setuptools import setup

setup(
    name='pyxyz',
    version='0.3{ver}.12',
    author='Nikolai Krivoshchapov',
    python_requires='==3.{ver}.*',
    packages=['pyxyz'],
    package_data={'pyxyz': ['__init__.py', 'confpool.so', 'libgslcblas.so.0', 'test/*']}
)
