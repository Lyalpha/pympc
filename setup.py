import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='pympc',
    version='0.6.0',
    author='Joe Lyman',
    description='minor planet checking',
    packages=setuptools.find_packages(),
    install_requires=[
        "ephem",
        "astropy>=3.2",
        "pandas>=1.0",
        "requests>=2.0",
    ],
    license='GNU General Public License v3 (GPLv3)',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lyalpha/pympc",
    download_url="https://github.com/Lyalpha/pympc/archive/v0.6.0.tar.gz",
    classifiers=[
         "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
)
