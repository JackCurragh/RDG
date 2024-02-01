
from setuptools import setup, find_packages
from pathlib import Path

with open("README.md") as readme_file:
    readme = readme_file.read()

# Get the path to the requirements.txt file
requirements_path = Path(__file__).parent / 'requirements.txt'

# Read the contents of the requirements file
with open(requirements_path) as f:
    requirements = f.read().splitlines()

test_requirements = [
    "pytest>=3",
]

setup(
    author="Jack Tierney",
    author_email="jackcurragh@gmail.com",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    description="""
    A python command-line utility for Ribosome Decision Graph (RDG) generation.
    """,
    entry_points={
        "console_scripts": [
            "RDG=RDG.cli:rdg_cli",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords="RDG, Ribosome Decision Graph",
    name="RDG",
    packages=find_packages(
    ),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/jackcurragh/RDG",
    version="0.1.6",
    zip_safe=False,
)
