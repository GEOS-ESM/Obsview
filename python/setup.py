from setuptools import setup, find_packages
setup(
    name="obsview",
    version="0.1.0",
    description="Python port of GEOS-ESM Obsview observation visualization tool",
    author="GEOS-ESM Team",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=["numpy>=1.20.0", "xarray>=0.18.0", "netcdf4>=1.5.0", "matplotlib>=3.3.0", "scipy>=1.6.0"],
    extras_require={"dev": ["pytest>=6.0", "pytest-cov>=2.12.0"]},
)
