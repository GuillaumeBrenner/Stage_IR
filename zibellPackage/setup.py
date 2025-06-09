from setuptools import setup, find_packages

setup(
    name='zibellpackage',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy'
    ],
    author='Stage-IR',
    description='Package de simulation et estimation pour le modèle ZI-Bell avec MDPDE',
)
