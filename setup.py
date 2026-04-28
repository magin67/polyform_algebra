from setuptools import setup, find_packages

setup(
    name='polyform_algebra',
    version='0.1',
    packages=find_packages(),
    install_requires=['sympy', 'numpy', 'networkx'],
)