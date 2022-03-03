import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pympc",
    author="Joe Lyman",
    description="minor planet checking",
    packages=setuptools.find_packages(),
    install_requires=["ephem>=3.7.7.1", "astropy>=4.2", "pandas>=1.0", "requests>=2.0"],
    license="GNU General Public License v3 (GPLv3)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lyalpha/pympc",
    entry_points={"console_scripts": ["minor_planet_check=pympc.pympc:_minor_planet_check"]},
    classifiers=["Programming Language :: Python :: 3"],
    python_requires=">=3.6",
)
