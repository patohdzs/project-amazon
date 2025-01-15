from setuptools import setup, find_packages

setup(
    name="project-amazon",
    version="0.1",
    packages=find_packages(include=["pysrc", "pysrc.*"]),  # Includes 'pysrc' and its submodules
    install_requires=[],
)
