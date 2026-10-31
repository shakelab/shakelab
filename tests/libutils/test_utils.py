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
Unit tests for :mod:`shakelab.libutils.utils`.

The tests cover empty-value detection, casting of scalars and
sequences, rounding without side effects, normal and log-normal
statistics, item selection (including the deprecated ``slice`` alias)
and the JSON-friendly (de)serialization round trip.
"""

from __future__ import annotations

import json
import math

import numpy as np
import pytest

from shakelab.libutils import utils
from shakelab.libutils.utils import (
    a_round,
    cast_value,
    deserialize_complex,
    is_empty,
    lin_stat,
    log_stat,
    none_check,
    pick,
    serialize_ndarray,
)


# ----------------------------------------------------------------------------
# Empty values
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value",
    [None, "", "None", "none", "NaN", "nan", [], (), np.array([]),
     float("nan"), np.float64("nan")],
)
def test_is_empty_true(value):
    """Null-like values are detected as empty."""
    assert is_empty(value)


@pytest.mark.parametrize(
    "value",
    [0, 0.0, False, "0", "text", [0], (None,), np.array([1.0, 2.0]),
     np.array([np.nan]), {}],
)
def test_is_empty_false(value):
    """Meaningful values, including arrays, are not empty."""
    assert not is_empty(value)


def test_none_check():
    """Empty values become None, others are returned unchanged."""
    assert none_check("") is None
    assert none_check(np.array([])) is None
    data = [1, 2]
    assert none_check(data) is data


# ----------------------------------------------------------------------------
# Casting
# ----------------------------------------------------------------------------


def test_cast_value_scalars():
    """Scalars and numeric strings are cast to the target type."""
    assert cast_value("3.14", float) == 3.14
    assert cast_value("7", int) == 7
    assert cast_value(2.5, str) == "2.5"
    assert cast_value("1+2j", complex) == 1 + 2j


def test_cast_value_returns_same_object_when_already_typed():
    """A value already of the target type is returned as is."""
    text = "abc"
    assert cast_value(text, str) is text


def test_cast_value_sequences():
    """Lists, tuples and arrays are cast element-wise to arrays."""
    out = cast_value(["1", "2", "3"], float)
    np.testing.assert_array_equal(out, [1.0, 2.0, 3.0])

    out = cast_value((1 + 2j, 3 - 4j), complex)
    np.testing.assert_array_equal(out, [1 + 2j, 3 - 4j])

    out = cast_value(np.array([1, 2]), float)
    assert out.dtype == float


def test_cast_value_empty_returns_default():
    """Empty inputs, arrays included, return the default."""
    assert cast_value(None, float) is None
    assert cast_value("nan", float, default=0.0) == 0.0
    assert cast_value(np.array([]), float, default=-1) == -1


def test_cast_value_invalid_type():
    """Unsupported containers raise TypeError."""
    with pytest.raises(TypeError):
        cast_value({"a": 1}, float)


# ----------------------------------------------------------------------------
# Rounding and statistics
# ----------------------------------------------------------------------------


def test_a_round_scalar():
    """Scalars are rounded to the requested decimals."""
    assert a_round(3.14159, 2) == 3.14


def test_a_round_preserves_container_and_input():
    """Sequences are rounded into a new container of the same type."""
    data = [1.23456, 2.34567]
    out = a_round(data, 2)
    assert out == [1.23, 2.35]
    assert data == [1.23456, 2.34567]

    assert a_round((1.23456, [2.34567]), 1) == (1.2, [2.3])

    arr = np.array([1.23456, 2.34567])
    np.testing.assert_allclose(a_round(arr, 3), [1.235, 2.346])
    assert arr[0] == 1.23456


def test_lin_stat():
    """Mean and standard deviation are computed along the first axis."""
    mn, sd = lin_stat([[1.0, 10.0], [3.0, 30.0]])
    np.testing.assert_allclose(mn, [2.0, 20.0])
    np.testing.assert_allclose(sd, [1.0, 10.0])


def test_log_stat():
    """Geometric mean and multiplicative standard deviation."""
    mn, sd = log_stat([1.0, 100.0])
    assert math.isclose(mn, 10.0)
    assert math.isclose(sd, 10.0)


# ----------------------------------------------------------------------------
# Item selection
# ----------------------------------------------------------------------------


def test_pick():
    """Items are selected by integer or iterable index."""
    data = ["a", "b", "c"]
    assert pick(data, 1) == ["b"]
    assert pick(data, [2, 0]) == ["c", "a"]
    assert pick(data, (0, 1)) == ["a", "b"]
    assert pick(data, np.int64(2)) == ["c"]


def test_pick_without_index_returns_all():
    """An empty or missing index selects all items."""
    data = ["a", "b"]
    assert pick(data) == data
    assert pick(data, []) == data


def test_slice_is_deprecated_alias():
    """The legacy slice name warns and behaves like pick."""
    with pytest.warns(DeprecationWarning):
        assert utils.slice(["a", "b"], 0) == ["a"]


# ----------------------------------------------------------------------------
# Serialization
# ----------------------------------------------------------------------------


def test_serialize_ndarray_is_json_compatible():
    """Arrays, NumPy scalars, tuples and complex become JSON types."""
    obj = {
        "a": np.array([[1.0, 2.0]]),
        "b": np.int64(3),
        "c": (1, 2),
        "d": np.array([1 + 2j]),
        "e": "text",
    }
    out = serialize_ndarray(obj)

    assert out == {
        "a": [[1.0, 2.0]],
        "b": 3,
        "c": [1, 2],
        "d": [{"real": 1.0, "imag": 2.0}],
        "e": "text",
    }
    json.dumps(out)


def test_serialization_round_trip():
    """Complex values survive serialization and deserialization."""
    obj = {"poles": [1 + 2j, -3j], "gain": 2.0}
    text = json.dumps(serialize_ndarray(obj))
    assert deserialize_complex(json.loads(text)) == obj


def test_deserialize_keeps_other_dicts():
    """Dicts that are not complex mappings are kept as dicts."""
    obj = {"real": 1.0, "imag": 2.0, "unit": "m"}
    assert deserialize_complex(obj) == obj
