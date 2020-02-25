# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, '', 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='generatr',
    version='0.3',
    python_requires='>3.7',
    description='Dynamic multi-loci/multi-repeat tract microsatellite reference sequence generator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/helloabunai/RefGeneratr',
    author='Alastair Maxwell/University of Glasgow',
    author_email='alastair.maxwell@glasgow.ac.uk',
    license='GPLv3',
    classifiers=[],
    keywords='XML FASTA Genetic-references Bioinformatics Data-analysis',
    packages=find_packages(exclude=['input',
									'lib',
									'generatr.egg-info',
									'build',
									'dist',
									'logs'
									]),
    install_requires=['lxml'],
    package_data={'': ['dtdvalidate/example_input.xml',
					   'dtdvalidate/xml_rules.dtd']},
	include_package_data=True,
    entry_points={
        'console_scripts': ['generatr=generatr.generatr:main',],
    },
)