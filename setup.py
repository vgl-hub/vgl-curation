from setuptools import setup


setup(
    name="ProcessCuration",
    version="0.1.0",
    author="Delphine Lariviere",
    description="Process manually curated genome assembly",
    url="https://github.com/Delphine-L/vgl-curation",
    scripts=['ProcessCuration/AGPcorrect.py','ProcessCuration/hap_split.py','ProcessCuration/unloc.py','ProcessCuration/chromosome_assignment.py','ProcessCuration/sak_generation.py'],
    license='MIT License',
    packages=['ProcessCuration'],
    install_requires=[
        "biopython >= 1.85",
        "pandas >= 2.3", 
        "natsort >= 8.4.0",
    ],
    long_description=open('README.rst').read(),
    entry_points={
        'console_scripts': [
            'AGPcorrect=AGPcorrect:main', 
            'hap_split=hap_split:main',
            'unloc=unloc:main',
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
