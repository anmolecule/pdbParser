from setuptools import setup, find_packages, install_requires

setup(
    name='pdbParser',
    packages=find_packages(),
    version=1.0
    author="Ozge Yoluk"
    email="ozgeyoluk@proteinart.net"
    install_requires=[
        "numpy",
        "biopython",
        "requests",
        "argparse"
        ]
)
