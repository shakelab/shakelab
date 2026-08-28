# ****************************************************************************
#
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
#
# ****************************************************************************
"""
Base classes and context containers for ShakeLab ground-motion models.

The GMPE interface separates rupture properties from site-dependent
properties. Rupture parameters describe one earthquake rupture and are
normally scalar. Site parameters and distance measures describe the
evaluation sites and are stored as one-dimensional arrays.

Each GMPE declares the parameters it requires. This lets provider layers
prepare only the physical quantities needed by the selected model while
keeping the GMPE independent from event, geometry, database, and mapping
objects used elsewhere in ShakeLab.
"""

from __future__ import annotations

import abc as _abc
import json as _json
import os as _os

from dataclasses import dataclass as _dataclass
from types import MappingProxyType as _MappingProxyType
from typing import Any as _Any
from typing import Dict as _Dict
from typing import Iterator as _Iterator
from typing import Mapping as _Mapping
from typing import Optional as _Optional

import numpy as _np


class _Context(_Mapping[str, _Any]):
    """
    Immutable mapping with attribute access to context parameters.
    """

    __slots__ = ("_data",)

    def __init__(self, **parameters: _Any):
        object.__setattr__(
            self,
            "_data",
            _MappingProxyType(dict(parameters)),
        )

    def __getitem__(self, key: str) -> _Any:
        return self._data[key]

    def __iter__(self) -> _Iterator[str]:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __getattr__(self, name: str) -> _Any:
        try:
            return self._data[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __repr__(self) -> str:
        values = ", ".join(
            f"{key}={value!r}"
            for key, value in self._data.items()
        )
        return f"{self.__class__.__name__}({values})"

    @property
    def parameters(self) -> tuple[str, ...]:
        """
        Return the names of the parameters stored in the context.
        """
        return tuple(self._data.keys())

    def to_dict(self) -> _Dict[str, _Any]:
        """
        Return a shallow dictionary copy of the context parameters.
        """
        return dict(self._data)


class RuptureContext(_Context):
    """
    Physical parameters describing one earthquake rupture.

    Typical fields include ``mag``, ``rake``, ``dip``, ``ztor``, ``width``
    and ``hypo_depth``. The class deliberately does not prescribe a fixed
    parameter schema because different GMPEs require different quantities.

    Rupture parameters are normally scalar. Their scientific validity and
    required presence are checked by the selected GMPE.
    """


class SitesContext(_Context):
    """
    Site-dependent parameters and distance measures for a set of sites.

    Parameters
    ----------
    size
        Number of evaluation sites. If omitted, the size is inferred from
        the first non-scalar parameter. If all supplied parameters are
        scalar, the context represents one site.
    **parameters
        Site properties and distance measures. Scalar values are broadcast
        to all sites. Array-like values must be one-dimensional and have
        length ``size``.

    Notes
    -----
    Values are stored as read-only NumPy arrays. This guarantees a common
    vectorized representation for quantities such as ``vs30``, ``rjb``,
    ``rrup`` and ``rx`` and prevents accidental in-place modification of
    shared input data.
    """

    __slots__ = ("_size",)

    def __init__(
        self,
        size: _Optional[int] = None,
        **parameters: _Any,
    ):
        site_count = self._resolve_size(
            size,
            parameters,
        )

        normalized: _Dict[str, _np.ndarray] = {}

        for name, value in parameters.items():
            normalized[name] = self._as_site_array(
                name,
                value,
                site_count,
            )

        super().__init__(**normalized)
        object.__setattr__(self, "_size", site_count)

    @property
    def size(self) -> int:
        """
        Return the number of evaluation sites.
        """
        return self._size

    @staticmethod
    def _resolve_size(
        size: _Optional[int],
        parameters: _Mapping[str, _Any],
    ) -> int:
        """
        Resolve and validate the number of sites.
        """
        if size is not None:
            if isinstance(size, bool) or not isinstance(
                size,
                (int, _np.integer),
            ):
                raise TypeError("SitesContext size must be an integer.")

            site_count = int(size)

            if site_count < 0:
                raise ValueError(
                    "SitesContext size must be non-negative."
                )

            return site_count

        lengths = set()

        for name, value in parameters.items():
            array = _np.asarray(value)

            if array.ndim == 0:
                continue

            if array.ndim != 1:
                raise ValueError(
                    f"SitesContext parameter {name!r} must be scalar "
                    "or one-dimensional."
                )

            lengths.add(array.size)

        if len(lengths) > 1:
            values = ", ".join(
                str(length)
                for length in sorted(lengths)
            )
            raise ValueError(
                "SitesContext parameters have incompatible lengths: "
                f"{values}."
            )

        if lengths:
            return lengths.pop()

        return 1

    @staticmethod
    def _as_site_array(
        name: str,
        value: _Any,
        size: int,
    ) -> _np.ndarray:
        """
        Convert one parameter to a read-only site array.
        """
        array = _np.asarray(value)

        if array.ndim == 0:
            array = _np.broadcast_to(
                array,
                (size,),
            ).copy()

        elif array.ndim == 1:
            if array.size != size:
                raise ValueError(
                    f"SitesContext parameter {name!r} has length "
                    f"{array.size}, expected {size}."
                )

            array = array.copy()

        else:
            raise ValueError(
                f"SitesContext parameter {name!r} must be scalar "
                "or one-dimensional."
            )

        array.setflags(write=False)
        return array


@_dataclass(frozen=True)
class GroundMotionResult:
    """
    Result returned by a ground-motion model.

    Parameters
    ----------
    mean
        Natural logarithm of the median intensity measure.
    sigma
        Total standard deviation of ln(IM), when available.
    tau
        Inter-event standard deviation of ln(IM), when available.
    phi
        Intra-event standard deviation of ln(IM), when available.

    Notes
    -----
    Values may be scalars or arrays. Conversion from logarithmic median to
    linear intensity measure belongs to the caller/provider layer.
    """

    mean: _Any
    sigma: _Optional[_Any] = None
    tau: _Optional[_Any] = None
    phi: _Optional[_Any] = None


class GMPE(metaclass=_abc.ABCMeta):
    """
    Base class for ground-motion prediction equation models.

    GMPE subclasses declare required and optional inputs in three semantic
    categories:

    ``REQUIRES_RUPTURE_PARAMETERS``
        Properties of the earthquake rupture.

    ``REQUIRES_SITE_PARAMETERS``
        Site properties independent from rupture-to-site geometry.

    ``REQUIRES_DISTANCES``
        Rupture-to-site distance measures.

    Optional counterparts use the same categories. Rupture parameters are
    supplied through :class:`RuptureContext`; site parameters and distances
    are supplied through :class:`SitesContext`.
    """

    _COEFF_FILE = None
    _COEFF_SET = "default"

    REFERENCE_VELOCITY = None
    MAGNITUDE_TYPE = None

    REQUIRES_RUPTURE_PARAMETERS = frozenset()
    REQUIRES_SITE_PARAMETERS = frozenset()
    REQUIRES_DISTANCES = frozenset()

    OPTIONAL_RUPTURE_PARAMETERS = frozenset()
    OPTIONAL_SITE_PARAMETERS = frozenset()
    OPTIONAL_DISTANCES = frozenset()

    def __init__(self):
        self._validate_parameter_declarations()

        self.import_coeff_from_json(
            self._COEFF_FILE,
            self._COEFF_SET,
        )

    @property
    def required_rupture_parameters(self) -> frozenset[str]:
        """
        Return required rupture parameters.
        """
        return frozenset(self.REQUIRES_RUPTURE_PARAMETERS)

    @property
    def optional_rupture_parameters(self) -> frozenset[str]:
        """
        Return optional rupture parameters.
        """
        return frozenset(self.OPTIONAL_RUPTURE_PARAMETERS)

    @property
    def required_site_parameters(self) -> frozenset[str]:
        """
        Return required non-distance site parameters.
        """
        return frozenset(self.REQUIRES_SITE_PARAMETERS)

    @property
    def optional_site_parameters(self) -> frozenset[str]:
        """
        Return optional non-distance site parameters.
        """
        return frozenset(self.OPTIONAL_SITE_PARAMETERS)

    @property
    def required_distances(self) -> frozenset[str]:
        """
        Return required distance measures.
        """
        return frozenset(self.REQUIRES_DISTANCES)

    @property
    def optional_distances(self) -> frozenset[str]:
        """
        Return optional distance measures.
        """
        return frozenset(self.OPTIONAL_DISTANCES)

    @property
    def required_site_fields(self) -> frozenset[str]:
        """
        Return all required fields stored in SitesContext.
        """
        return (
            self.required_site_parameters
            | self.required_distances
        )

    @property
    def optional_site_fields(self) -> frozenset[str]:
        """
        Return all optional fields stored in SitesContext.
        """
        return (
            self.optional_site_parameters
            | self.optional_distances
        )

    def validate_contexts(
        self,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> None:
        """
        Validate contexts against the requirements of the GMPE.

        Raises
        ------
        TypeError
            If the supplied context objects have the wrong type.
        ValueError
            If one or more required parameters are missing.
        """
        if not isinstance(rupture, RuptureContext):
            raise TypeError(
                "rupture must be a RuptureContext instance."
            )

        if not isinstance(sites, SitesContext):
            raise TypeError(
                "sites must be a SitesContext instance."
            )

        missing_rupture = sorted(
            self.required_rupture_parameters
            - rupture.keys()
        )

        missing_sites = sorted(
            self.required_site_fields
            - sites.keys()
        )

        messages = []

        if missing_rupture:
            messages.append(
                "rupture: " + ", ".join(missing_rupture)
            )

        if missing_sites:
            messages.append(
                "sites: " + ", ".join(missing_sites)
            )

        if messages:
            raise ValueError(
                f"{self.__class__.__name__} is missing required "
                "context parameter(s): "
                + "; ".join(messages)
                + "."
            )

    @_abc.abstractmethod
    def ground_motion(
        self,
        imt: str,
        rupture: RuptureContext,
        sites: SitesContext,
    ) -> GroundMotionResult:
        """
        Evaluate one intensity-measure type for one rupture and many sites.

        Implementations should call :meth:`validate_contexts` before
        accessing required parameters.
        """
        raise NotImplementedError

    def import_coeff_from_json(
        self,
        json_file: str,
        coeff_set: str = "default",
    ) -> None:
        """
        Load model coefficients from a JSON file.

        The coefficient file must be stored in the ``data`` directory next
        to the GMPE package and follow the existing ShakeLab schema.
        """
        if not json_file:
            raise ValueError(
                f"{self.__class__.__name__} does not define _COEFF_FILE."
            )

        full_path = _os.path.dirname(__file__)
        path_file = _os.path.join(
            full_path,
            "data",
            json_file,
        )

        with open(path_file, encoding="utf-8") as jf:
            data = _json.load(jf)["coefficients"][coeff_set]

        self.keys = list(data["keys"])
        self.coeff: _Dict[str, list[float]] = {}

        for imt, values in data["type"].items():
            self.coeff[imt] = [
                float(value)
                for value in values
            ]

    def get_coefficients(
        self,
        imt: str,
    ) -> _Dict[str, float]:
        """
        Return coefficients corresponding to an intensity-measure type.
        """
        try:
            values = self.coeff[imt]
        except KeyError as exc:
            available = ", ".join(self.list_imts())
            raise KeyError(
                f"Invalid intensity measure type {imt!r}. "
                f"Available: {available}."
            ) from exc

        return dict(zip(self.keys, values))

    def list_imts(self) -> list[str]:
        """
        Return the intensity-measure types supported by the GMPE.
        """
        return list(self.coeff.keys())

    def _validate_parameter_declarations(self) -> None:
        """
        Validate class-level context parameter declarations.
        """
        groups = {
            "REQUIRES_RUPTURE_PARAMETERS":
                self.REQUIRES_RUPTURE_PARAMETERS,
            "REQUIRES_SITE_PARAMETERS":
                self.REQUIRES_SITE_PARAMETERS,
            "REQUIRES_DISTANCES":
                self.REQUIRES_DISTANCES,
            "OPTIONAL_RUPTURE_PARAMETERS":
                self.OPTIONAL_RUPTURE_PARAMETERS,
            "OPTIONAL_SITE_PARAMETERS":
                self.OPTIONAL_SITE_PARAMETERS,
            "OPTIONAL_DISTANCES":
                self.OPTIONAL_DISTANCES,
        }

        normalized = {}

        for name, values in groups.items():
            try:
                values = frozenset(values)
            except TypeError as exc:
                raise TypeError(
                    f"{name} must be an iterable of strings."
                ) from exc

            if any(
                not isinstance(value, str) or not value.strip()
                for value in values
            ):
                raise ValueError(
                    f"{name} must contain non-empty strings."
                )

            normalized[name] = values

        required_rupture = normalized[
            "REQUIRES_RUPTURE_PARAMETERS"
        ]
        optional_rupture = normalized[
            "OPTIONAL_RUPTURE_PARAMETERS"
        ]
        required_site = normalized[
            "REQUIRES_SITE_PARAMETERS"
        ]
        optional_site = normalized[
            "OPTIONAL_SITE_PARAMETERS"
        ]
        required_distances = normalized[
            "REQUIRES_DISTANCES"
        ]
        optional_distances = normalized[
            "OPTIONAL_DISTANCES"
        ]

        overlap = required_rupture & optional_rupture
        self._raise_overlap(
            overlap,
            "required and optional rupture parameters",
        )

        overlap = required_site & optional_site
        self._raise_overlap(
            overlap,
            "required and optional site parameters",
        )

        overlap = required_distances & optional_distances
        self._raise_overlap(
            overlap,
            "required and optional distances",
        )

        rupture_fields = required_rupture | optional_rupture
        site_fields = (
            required_site
            | optional_site
            | required_distances
            | optional_distances
        )

        self._raise_overlap(
            rupture_fields & site_fields,
            "rupture and site context fields",
        )

        self._raise_overlap(
            (required_site | optional_site)
            & (required_distances | optional_distances),
            "site parameters and distances",
        )

    def _raise_overlap(
        self,
        overlap: frozenset[str],
        description: str,
    ) -> None:
        """
        Raise for ambiguous parameter declarations.
        """
        if not overlap:
            return

        names = ", ".join(sorted(overlap))
        raise ValueError(
            f"{self.__class__.__name__} has overlapping "
            f"{description}: {names}."
        )


__all__ = [
    "GMPE",
    "RuptureContext",
    "SitesContext",
    "GroundMotionResult",
]
