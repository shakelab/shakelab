import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='shakelab',
    version='0.0.1',
    author='Valerio Poggi',
    author_email='vpoggi@inogs.it',
    description='Python tools for engineering seismology',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/shakelab/shakelab',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research'
    ],
    python_requires='>=3.6',
    install_requires=[
        'setuptools',
        'numpy',
        'scipy',
        'shapely',
        'matplotlib'
    ]
)
