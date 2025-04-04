<img alt="ShakeLab - Tools for Engineering Seismology" class="right" style="width: 100%" src="https://raw.githubusercontent.com/klunk386/shakelab/master/logo/shakelab-header.png#1" />

[![AGPLv3](https://www.gnu.org/graphics/agplv3-88x31.png)](https://www.gnu.org/licenses/agpl.html)
[![Build Status](https://travis-ci.org/shakelab/shakelab.svg?branch=master)](https://travis-ci.org/shakelab/shakelab)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/shakelab.svg)](https://pypi.python.org/pypi/shakelab)
[![PyPI Version](https://img.shields.io/pypi/v/shakelab.svg)](https://pypi.python.org/pypi/shakelab)

# ShakeLab

**ShakeLab** is an open-source Python library that provides tools for engineering seismology.  
It is designed to assist researchers and practitioners in the analysis and interpretation of ground motion data, generation of shakemaps, site response calculations, and rapid seismic impact assessment.

⚠️ **Note:** _The project is under active development and some functionalities may change._

##  Installation

You can install ShakeLab in different ways depending on your needs.

### 1️⃣ Standard Installation from PyPI

Install the latest stable version from [PyPI](https://pypi.org/project/shakelab/):

```bash
pip install shakelab
```

or explicitly with Python 3:

```bash
pip3 install shakelab
```

**Note:** The version available on PyPI may not include the latest features and fixes as the code is under active development. If you want the most up-to-date version, consider using one of the installation methods below.

### 2️⃣ Installation from GitHub (latest code)

To install the most recent version directly from the GitHub repository:

```bash
pip install git+https://github.com/shakelab/shakelab.git
```

To upgrade an existing installation:

```bash
pip install --upgrade git+https://github.com/shakelab/shakelab.git
```

### 3️⃣ Developer Installation from a Local Repository

If you want to modify the code and test it locally:

```bash
git clone https://github.com/shakelab/shakelab.git
cd shakelab
pip install -r requirements.txt -e .
```

This will install the package in "editable mode", so any changes to the source code will be reflected immediately.

> ✅ Optionally, create and activate a virtual environment before installation:
> ```bash
> python3 -m venv venv
> source venv/bin/activate  # On Windows: venv\Scripts\activate
> ```

##  Dependencies

ShakeLab relies on the following Python packages:

- [NumPy](https://numpy.org/) >= 1.24.0
- [SciPy](https://scipy.org/) >= 1.10.0
- [Shapely](https://pypi.org/project/Shapely/) >= 2.0.0
- [Requests](https://pypi.org/project/requests/) >= 2.32.0
- [cymseed3](https://pypi.org/project/cymseed3/) >= 0.1.4


## Documentation

Official documentation is available at: [shakelab.org](http://shakelab.org)

##  License

ShakeLab is free software, released under the **GNU Affero General Public License v3.0**.  
You can redistribute it and/or modify it under the terms of the license.

See <https://www.gnu.org/licenses/agpl-3.0.html> for full details.

## ⚠️ Disclaimer

ShakeLab is distributed in the hope that it will be useful, but **without any warranty**, including the implied warranty of merchantability or fitness for a particular purpose.

The authors assume no liability for any use of the software.
