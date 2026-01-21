from setuptools import setup, find_packages


setup(
    name="ProcessCuration",
    version="1.1",
    author="Delphine Lariviere",
    description="Process manually curated genome assembly",
    url="https://github.com/Delphine-L/vgl-curation",
    scripts=['src/split_agp.py','src/chromosome_assignment.py','src/sak_generation.py'],
    license='MIT License',    
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "biopython >= 1.85",
        "pandas >= 2.3", 
        "natsort >= 8.4.0",
    ],
    long_description=open('README.rst').read(),
    entry_points={
        'console_scripts': [
            'split_agp=split_agp:main', 
            'chromosome_assignment=chromosome_assignment:main',
            'sak_generation=sak_generation:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
