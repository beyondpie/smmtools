import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="smmuty",
    version="0.0.1",
    author="Songpeng Zu",
    author_email="szu@health.ucsd.edu",
    description="utilities for single-cell multi-omics data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/beyondpie/smmtools/inst",
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
        "pysam"
    ],
    python_requires=">=3.6",
)
