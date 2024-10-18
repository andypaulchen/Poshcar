from setuptools import setup, find_packages

setup(
    name="poshcar",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "ase",
        "pyvis",
        "rdkit",
        "scipy",
        "pymatgen",
        "chgnet"
    ],
)
