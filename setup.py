from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    readme = f.read()


setup(
    name='BGWAS',

    # Version:
    version='0.0.1',

    description='Helper functions for GWAS analysis',
    long_description=readme,

    # The project's main homepage.
    url='https://github.com/b-schubert/BGWAS',

    # Author details
    author='Benjamin Schubert',
    author_email='benni.schubert@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: GWAS, Genome analysis, Developer',
        'Topic :: Genomics :: Analysis functions',

        # The license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # The supported Python versions (other than development version were
        # not yet tested. Especially we should check for Python 3 support
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],

    # What FRED2 relates to:
    keywords='GWAS',

    # Specify  packages via find_packages() and exclude the tests and
    # documentation:
    packages=find_packages(),

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #include_package_data=True,

    #package_data is a lie: http://stackoverflow.com/questions/7522250/how-to-include-package-data-with-setuptools-distribute

    # 'package_data' is used to also install non package data files
    # see http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    # example:
    #data_files=data_files,

    # Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    # IMPORTANT: script names need to be in lower case ! ! ! (otherwise
    # deinstallation does not work)
    #entry_points={
    #    'console_scripts': [
    #        'epitopeprediction=Fred2.Apps.EpitopePrediction:main',
    #    ],
    #},



    # Run-time dependencies. (will be installed by pip when FRED2 is installed)
    install_requires=['setuptools>=18.2', 'pandas', "intervaltree", "scikit-learn",
                      "statsmodels", "matplotlib", "seaborn", "plotly", "numpy", "scipy"],

)