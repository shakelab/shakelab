# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# ****************************************************************************
"""
"""
import numpy as _np


def cast_value(value, dtype, default=None):
    """
    Cast a value to the specified data type, handling both scalars and
    sequences.

    Parameters
    ----------
    value : any
        The input value to cast. It can be a scalar, list, or NumPy array.
    dtype : type
        The target data type (e.g., int, float, complex). Applied directly
        to scalars or element-wise to sequences.
    default : any, optional
        The value to return if `value` is considered empty (None, '', etc.).

    Returns
    -------
    casted : dtype or np.ndarray
        A scalar or NumPy array with all elements cast to the target dtype.

    Examples
    --------
    >>> cast_value("3.14", float)
    3.14
    >>> cast_value([1, 2, 3], float)
    array([1., 2., 3.])
    >>> cast_value([1+2j, 3-4j], complex)
    array([1.+2.j, 3.-4.j])
    >>> cast_value(None, float, default=0.0)
    0.0
    """
    # Handle null-like values
    if value in [None, 'None', 'none', 'NaN', 'nan', '', []]:
        return default

    # If already of target type, return as-is
    if isinstance(value, dtype):
        return value

    # Scalar casting
    if _np.isscalar(value):
        return dtype(value)

    # Sequence casting (element-wise)
    if isinstance(value, (list, tuple, _np.ndarray)):
        return _np.array([dtype(v) for v in value])

    raise TypeError(f"Cannot cast value of type {type(value)}")


def a_round(number, decimals=3):
    """
    Round value or an array to a given decimal place
    (it works also for single scalar and list)

    :param [float, list, tuple, numpy.ndarray] number:
        Input; can be a number or a sequence

    :param int decimals:
        Rounding decimals

    :return [float, list, tuple, numpy.ndarray] number:
        Rounded output
    """
    if isinstance(number, (list, tuple, _np.ndarray)):
        for i, n in enumerate(number):
            number[i] = round(n, decimals)
    else:
        number = round(number, decimals)

    return number


def lin_stat(data):
    """
    Compute mean and standard deviation of data array
    assuming normal (linear) statistic.

    :param list or numpy.ndarray data:
        The input dataset

    :return float (mn, sd):
        Mean and standard deviation
    """
    mn = _np.mean(data, axis=0)
    sd = _np.std(data, axis=0)

    return (mn, sd)


def log_stat(data):
    """
    Compute mean and standard deviation of data array
    assuming log-normal statistic.

    :param list or numpy.ndarray data:
        The input dataset

    :return float (mn, sd):
        Mean and standard deviation

    """
    mn = _np.exp(_np.mean(_np.log(data), axis=0))
    sd = _np.exp(_np.std(_np.log(data), axis=0))

    return (mn, sd)


def slice(data, index=[]):
    """
    Extract items from list using an index list.

    :param list data:
        Input data to be sliced

    :param int or list index:
        Indexes of the list items

    :return list data_slice:
        List with the sliced data
    """
    if is_empty(index):
        return data
    else:
        if not isinstance(index, list):
            index = [index]
        data_slice = [data[i] for i in index]

    return data_slice


def is_empty(value):
    """
    Boolean check if variable of arbitrary format is empty or None.

    :param [arbitrary] value:
        Arbitrary input value for the check

    :return boolean empty:
        True if no or empty value is given
    """
    C0 = (value == [])
    C1 = (value == '')
    C2 = (value == None)
    C3 = (value == 'None')
    C4 = (value != value)

    empty = (C0 or C1 or C2 or C3 or C4)

    return empty


def none_check(number):
    """
    Replace with none if empty.
    """
    number = None if is_empty(number) else number

    return number

def serialize_ndarray(obj):
    """
    Recursively convert NumPy arrays and complex numbers
    for JSON serialization.
    """
    if isinstance(obj, dict):
        return {k: serialize_ndarray(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [serialize_ndarray(v) for v in obj]
    elif isinstance(obj, _np.ndarray):
        return serialize_ndarray(obj.tolist())
    elif isinstance(obj, complex):
        return {'real': obj.real, 'imag': obj.imag}
    return obj

def deserialize_complex(obj):
    """
    Recursively convert {'real': x, 'imag': y} back to complex numbers.
    """
    if isinstance(obj, dict):
        if 'real' in obj and 'imag' in obj and len(obj) == 2:
            return complex(obj['real'], obj['imag'])
        return {k: deserialize_complex(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [deserialize_complex(v) for v in obj]
    return obj