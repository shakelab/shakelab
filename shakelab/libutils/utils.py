# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
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
# ****************************************************************************
"""
General-purpose helpers shared across ShakeLab.

The module groups small utilities for value casting, empty-value
checks, rounding, basic statistics, item selection and JSON-friendly
(de)serialization of NumPy arrays and complex numbers.

Conventions
-----------
- A value is *empty* when it is ``None``, NaN, an empty string, one of
  the null markers ``'None'``, ``'none'``, ``'NaN'``, ``'nan'``, or an
  empty list, tuple or array (see :func:`is_empty`).
- Functions never modify their inputs in place.
"""

from __future__ import annotations

import warnings
from typing import Any, Callable, Sequence

import numpy as np

__all__ = [
    "cast_value",
    "a_round",
    "lin_stat",
    "log_stat",
    "pick",
    "slice",
    "is_empty",
    "none_check",
    "serialize_ndarray",
    "deserialize_complex",
]

_NULL_STRINGS = frozenset({"", "None", "none", "NaN", "nan"})


# ----------------------------------------------------------------------------
# Empty values and casting
# ----------------------------------------------------------------------------


def is_empty(value: Any) -> bool:
    """
    Check whether a value of arbitrary type is empty.

    Parameters
    ----------
    value
        Value to check.

    Returns
    -------
    bool
        True for ``None``, NaN scalars, null-marker strings and empty
        lists, tuples or arrays; False otherwise.
    """
    if value is None:
        return True

    if isinstance(value, str):
        return value in _NULL_STRINGS

    if isinstance(value, np.ndarray):
        return value.size == 0

    if isinstance(value, (list, tuple)):
        return len(value) == 0

    try:
        return bool(value != value)
    except (TypeError, ValueError):
        return False


def none_check(value: Any) -> Any:
    """
    Replace an empty value with ``None``.

    Parameters
    ----------
    value
        Value to check.

    Returns
    -------
    Any
        ``None`` if ``value`` is empty, otherwise ``value`` unchanged.
    """
    return None if is_empty(value) else value


def cast_value(value: Any, dtype: Callable, default: Any = None) -> Any:
    """
    Cast a scalar or a sequence to a target type.

    Parameters
    ----------
    value
        Scalar, list, tuple or NumPy array to cast.
    dtype
        Target type (e.g. ``int``, ``float``, ``complex``, ``str``),
        applied to scalars or element-wise to sequences.
    default
        Value returned when ``value`` is empty (see :func:`is_empty`).

    Returns
    -------
    Any
        The cast scalar, or a NumPy array for sequence inputs.

    Raises
    ------
    TypeError
        If ``value`` is neither a scalar nor a sequence.

    Examples
    --------
    >>> cast_value("3.14", float)
    3.14
    >>> cast_value([1, 2, 3], float)
    array([1., 2., 3.])
    >>> cast_value(None, float, default=0.0)
    0.0
    """
    if is_empty(value):
        return default

    if isinstance(dtype, type) and isinstance(value, dtype):
        return value

    if np.isscalar(value):
        return dtype(value)

    if isinstance(value, (list, tuple, np.ndarray)):
        return np.array([dtype(v) for v in value])

    raise TypeError(f"Cannot cast value of type {type(value).__name__}")


# ----------------------------------------------------------------------------
# Rounding and statistics
# ----------------------------------------------------------------------------


def a_round(number: Any, decimals: int = 3) -> Any:
    """
    Round a scalar or a sequence to a given number of decimals.

    Parameters
    ----------
    number
        Scalar, list, tuple or NumPy array. Nested lists and tuples are
        rounded recursively.
    decimals
        Number of decimal places.

    Returns
    -------
    Any
        Rounded copy of ``number`` with the same container type.
    """
    if isinstance(number, np.ndarray):
        return np.round(number, decimals)

    if isinstance(number, (list, tuple)):
        return type(number)(a_round(n, decimals) for n in number)

    return round(number, decimals)


def lin_stat(data: Any) -> tuple[Any, Any]:
    """
    Compute mean and standard deviation assuming normal statistics.

    Parameters
    ----------
    data
        Input samples; statistics are computed along the first axis.

    Returns
    -------
    tuple
        Arithmetic mean and (population) standard deviation.
    """
    data = np.asarray(data, dtype=float)
    return np.mean(data, axis=0), np.std(data, axis=0)


def log_stat(data: Any) -> tuple[Any, Any]:
    """
    Compute mean and standard deviation assuming log-normal statistics.

    Parameters
    ----------
    data
        Strictly positive input samples; statistics are computed along
        the first axis.

    Returns
    -------
    tuple
        Geometric mean and geometric standard deviation (a
        multiplicative factor, ``exp`` of the standard deviation of the
        natural logarithm).
    """
    logs = np.log(np.asarray(data, dtype=float))
    return np.exp(np.mean(logs, axis=0)), np.exp(np.std(logs, axis=0))


# ----------------------------------------------------------------------------
# Item selection
# ----------------------------------------------------------------------------


def pick(data: Sequence[Any], index: Any = None) -> list[Any]:
    """
    Select items from a sequence by index.

    Parameters
    ----------
    data
        Input sequence.
    index
        Integer index or iterable of indices. If empty or ``None``,
        all items are returned.

    Returns
    -------
    list
        Selected items, in the order given by ``index``.
    """
    if is_empty(index):
        return list(data)

    if isinstance(index, (int, np.integer)):
        index = [index]

    return [data[i] for i in index]


def slice(data: Sequence[Any], index: Any = None) -> list[Any]:
    """
    Deprecated alias of :func:`pick`.

    .. deprecated::
        Use :func:`pick`; this name shadows the built-in ``slice``.
    """
    warnings.warn(
        "shakelab.libutils.utils.slice is deprecated; use pick",
        DeprecationWarning,
        stacklevel=2,
    )
    return pick(data, index)


# ----------------------------------------------------------------------------
# JSON-friendly serialization
# ----------------------------------------------------------------------------


def serialize_ndarray(obj: Any) -> Any:
    """
    Recursively convert NumPy objects and complex numbers for JSON.

    Arrays and tuples become lists, NumPy scalars become Python scalars
    and complex numbers become ``{'real': x, 'imag': y}`` mappings.

    Parameters
    ----------
    obj
        Object to convert (dicts, lists and tuples are traversed).

    Returns
    -------
    Any
        JSON-serializable copy of ``obj``.
    """
    if isinstance(obj, dict):
        return {k: serialize_ndarray(v) for k, v in obj.items()}

    if isinstance(obj, (list, tuple)):
        return [serialize_ndarray(v) for v in obj]

    if isinstance(obj, np.ndarray):
        return serialize_ndarray(obj.tolist())

    if isinstance(obj, np.generic):
        return serialize_ndarray(obj.item())

    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}

    return obj


def deserialize_complex(obj: Any) -> Any:
    """
    Recursively convert ``{'real': x, 'imag': y}`` mappings to complex.

    Parameters
    ----------
    obj
        Object to convert (dicts and lists are traversed).

    Returns
    -------
    Any
        Copy of ``obj`` with complex numbers restored.
    """
    if isinstance(obj, dict):
        if set(obj) == {"real", "imag"}:
            return complex(obj["real"], obj["imag"])
        return {k: deserialize_complex(v) for k, v in obj.items()}

    if isinstance(obj, list):
        return [deserialize_complex(v) for v in obj]

    return obj
