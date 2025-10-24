import setuptools
from pathlib import Path

with Path("README.md").open(encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="saltshaker",
    version="0.1.0",
    author="OUS AMG",
    description="Pattern classification and visualization for MitoSAlt",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aksenia/mitosalt_custom",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "matplotlib>=3.3.0",
        "biopython>=1.78",
        "logistro>=2.0.0"
    ],
    packages=setuptools.find_packages(),
    package_data={
        'saltshaker': ['data/*.bed'],
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "saltshaker=saltshaker.__main__:main",
        ],
    },
    python_requires=">=3.8",
)