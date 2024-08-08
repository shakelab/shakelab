from setuptools import setup
from setuptools import Extension
from Cython.Build import cythonize

libmseed_path = "EarthScope/libmseed-3.1.3"

extensions = [
    Extension(
        name="libmseed",
        sources=["libmseed.pyx"],
        libraries=["mseed"],
        library_dirs=[libmseed_path],
        include_dirs=[libmseed_path]
    )
]

setup(
    name="miniseed",
    ext_modules=cythonize(extensions)
)
