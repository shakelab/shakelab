import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="quakelab",
    version="0.0.1",
    author="Valerio Poggi",
    author_email="vpoggi@inogs.it",
    description="Python tools for engineering seismology",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/klunk386/QuakeLab",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL License",
        "Operating System :: OS Independent",
    ],
)