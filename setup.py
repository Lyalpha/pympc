import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pympc",
    author="Joe Lyman",
    description="minor planet checking",
    packages=setuptools.find_packages(),
    install_requires=[
        "ephem>=4.1",
        "astropy>=4.2",
        "numpy>=1.18",
        "pandas>=1.0",
        "pyerfa>=2.0.0.0",
        "requests",
        "rich>=14.0.0",
        "loguru>=0.7.0"
    ],
    license="GNU General Public License v3 (GPLv3)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lyalpha/pympc",
    entry_points={
        "console_scripts": ["minor_planet_check=pympc.pympc:_console_script"]
    },
    classifiers=["Programming Language :: Python :: 3"],
    python_requires=">=3.8",
)
