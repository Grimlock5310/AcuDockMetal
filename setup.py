"""Minimal setup for AcuDockMetal package."""

from pathlib import Path
from setuptools import setup, find_packages

requirements = Path("requirements.txt").read_text().strip().splitlines()

setup(
    name="acudockmetal",
    version="0.1.0",
    description="Metal-Aware epsilon-Certified Molecular Docking System",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=requirements,
)
