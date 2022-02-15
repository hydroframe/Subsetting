import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="parflow_subsetter",
    version="1.0.0",
    author="HydroFrame Team",
    author_email="parflow@parflow.org",
    description="A set of tools for clipping ParFlow model inputs and outputs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hydroframe/Subsetting",
    packages=setuptools.find_packages(exclude=("tests",)),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['pyyaml>=5.3.0', 'pandas>=1.0', 'parflowio>=0.0.4', 'pftools>=1.0.0'],
    package_data={
        # Include any *.yaml, *.tcl files found in the "data" subdirectory
        # of the "pfsubset.subset" package, also:
        "pfsubset.subset": ["data/*.yaml", "data/*.tcl"],
    }
)
