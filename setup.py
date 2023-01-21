from setuptools import setup
"""
Description of how to make a python package

https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html

"""


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="pyDocking",
    version="0.1",
    long_description=readme(),
    description="A collection of drug discovery tools",
    url="https://github.com/tuhinmallick/pyDocking",
    author="Tuhin Mallick",
    author_email="tuhin.mllk@gmail.com",
    license="GPL-3.0",
    packages=["pyDocking"],
    install_requires=[
        "numpy",
        "pandas",
        "pubchempy",
        "mdtraj",
        "rdkit",
    ],
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.5",
)
