import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="smmuty",
    version="0.0.3",
    author="Songpeng Zu",
    author_email="szu@health.ucsd.edu",
    description="Utilities for single-cell multi-omics data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/beyondpie/smmtools/tree/main/inst/smmuty",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        "leidenalg",
        "numpy",
        "scrublet",
        "sklearn",
        "igraph",
        "pysam",
        "scipy",
        "matplotlib",
        "fastcluster"
    ],
    python_requires=">=3.6",
)
