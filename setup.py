from setuptools import setup, find_packages

setup(
    name="shakelab",
    version="0.0.2",
    packages=find_packages(include=["shakelab", "shakelab.*"]),
    scripts=['scripts/fdsnclient.py'],
)

