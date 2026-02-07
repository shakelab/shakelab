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
Fragility models and collections for seismic impact scenarios.

This module provides a compact representation of fragility information used
to compute probabilities of exceedance (PoE) as a function of an intensity
measure (IM).

Design overview
---------------
The implementation is organized in three layers.

1) Single-damage-state components
   - BaseDamageModel: abstract PoE(IM) interface for one damage state
   - LognormalDamageModel: continuous lognormal model (theta, beta)
   - DiscreteDamageModel: tabulated model with interpolation

2) Multi-state fragility model
   - FragilityModel: model "by author/source" for a given taxonomy and IMT,
     containing multiple damage states (custom scale). Damage states are
     identified by strings (e.g., "D1", "D2", ...) with an explicit order.

3) Model collection
   - FragilityCollection: container of multiple FragilityModel instances,
     retrievable by model id and loadable from a JSON database.

JSON database
-------------
The expected JSON database uses a stable top-level structure:
- type: "ShakeLabFragility"
- schema_version: "1.0.0"
- metadata: object (at least name, date)
- models: array of model objects

No plotting utilities are included in this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import norm


__all__ = [
    "BaseDamageModel",
    "LognormalDamageModel",
    "DiscreteDamageModel",
    "DamageScale",
    "IMBounds",
    "FragilityModel",
    "FragilityCollection",
]


_SCHEMA_TYPE = "ShakeLabFragility"
_SCHEMA_VERSION = "1.0.0"
_DEFAULT_TAIL_KEY = "GT_LAST"


def _is_number(x: Any) -> bool:
    """Return True for int/float-like values (excluding bool)."""
    return isinstance(x, (int, float, np.number)) and not isinstance(x, bool)


def _as_1d_float(x: Any, name: str) -> np.ndarray:
    """Convert input to a 1D float array."""
    arr = np.asarray(x, dtype=float)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be scalar or 1D array-like.")
    return arr


def _clip01(p: np.ndarray) -> np.ndarray:
    """Clip probabilities to [0, 1]."""
    return np.clip(p, 0.0, 1.0)


def _require_str(dct: Mapping[str, Any], key: str, path: str) -> str:
    """Read a required non-empty string."""
    val = dct.get(key)
    if not isinstance(val, str) or not val.strip():
        raise ValueError(f"{path}.{key} must be a non-empty string.")
    return val.strip()


def _require_dict(dct: Mapping[str, Any], key: str, path: str) -> Dict[str, Any]:
    """Read a required object."""
    val = dct.get(key)
    if not isinstance(val, dict):
        raise ValueError(f"{path}.{key} must be an object.")
    return val


class BaseDamageModel:
    """
    Abstract base class for single-damage-state models.

    Subclasses implement poe(im) returning the probability of exceedance
    for the associated damage state.
    """

    def poe(self, im: Any) -> np.ndarray:
        """
        Probability of exceedance for deterministic IM values.

        Parameters
        ----------
        im
            Intensity measure values (scalar or 1D array-like).

        Returns
        -------
        numpy.ndarray
            PoE values in [0, 1], same length as input (or length 1).
        """
        raise NotImplementedError

    def poe_lognormal_im(self, im_median: float, im_beta: float) -> float:
        """
        Convolve PoE with a lognormal IM uncertainty (generic fallback).

        Parameters
        ----------
        im_median
            Median IM (> 0).
        im_beta
            Log-standard deviation of IM in natural log units (>= 0).

        Returns
        -------
        float
            Expected PoE in [0, 1].
        """
        if not _is_number(im_median) or im_median <= 0.0:
            raise ValueError("im_median must be a number > 0.")
        if not _is_number(im_beta) or im_beta < 0.0:
            raise ValueError("im_beta must be a number >= 0.")

        if im_beta == 0.0:
            return float(self.poe(im_median)[0])

        npts = 121
        z = np.linspace(-6.0, 6.0, npts)
        im = float(im_median) * np.exp(float(im_beta) * z)
        w = norm.pdf(z)
        poe_vals = self.poe(im)
        val = float(np.trapz(poe_vals * w, z))
        return float(np.clip(val, 0.0, 1.0))

    def poe_distribution(
        self,
        im_values: Sequence[float],
        probabilities: Sequence[float],
    ) -> float:
        """
        Convolve PoE with a discrete IM distribution.

        Parameters
        ----------
        im_values
            IM support values (N).
        probabilities
            Probabilities for each IM value (N). They will be normalized
            if they do not sum to 1.

        Returns
        -------
        float
            Expected PoE in [0, 1].
        """
        im = _as_1d_float(im_values, "im_values")
        pr = _as_1d_float(probabilities, "probabilities")

        if im.size != pr.size:
            raise ValueError("im_values and probabilities must have same size.")
        if im.size == 0:
            raise ValueError("im_values must not be empty.")
        if np.any(~np.isfinite(im)) or np.any(~np.isfinite(pr)):
            raise ValueError("im_values/probabilities must be finite.")
        if np.any(pr < 0.0):
            raise ValueError("probabilities must be >= 0.")

        ssum = float(pr.sum())
        if ssum <= 0.0:
            raise ValueError("probabilities sum must be > 0.")
        pr = pr / ssum

        poe_vals = self.poe(im)
        val = float(np.sum(pr * poe_vals))
        return float(np.clip(val, 0.0, 1.0))


@dataclass(frozen=True)
class LognormalDamageModel(BaseDamageModel):
    """
    Continuous lognormal single-damage-state model.

    Parameters
    ----------
    theta
        Median capacity (same units as IM), must be > 0.
    beta
        Log-standard deviation in natural log units, must be > 0.
    """

    theta: float
    beta: float

    def __post_init__(self) -> None:
        if not _is_number(self.theta) or self.theta <= 0.0:
            raise ValueError("theta must be a number > 0.")
        if not _is_number(self.beta) or self.beta <= 0.0:
            raise ValueError("beta must be a number > 0.")

    def poe(self, im: Any) -> np.ndarray:
        im_arr = _as_1d_float(im, "im")
        out = np.zeros_like(im_arr, dtype=float)

        msk = im_arr > 0.0
        if np.any(msk):
            z = (np.log(im_arr[msk] / float(self.theta))) / float(self.beta)
            out[msk] = norm.cdf(z)

        return _clip01(out)

    def poe_lognormal_im(self, im_median: float, im_beta: float) -> float:
        """
        Closed-form convolution for lognormal IM uncertainty.

        Parameters
        ----------
        im_median
            Median IM (> 0).
        im_beta
            Log-standard deviation of IM in natural log units (>= 0).

        Returns
        -------
        float
            Expected PoE in [0, 1].
        """
        if not _is_number(im_median) or im_median <= 0.0:
            raise ValueError("im_median must be a number > 0.")
        if not _is_number(im_beta) or im_beta < 0.0:
            raise ValueError("im_beta must be a number >= 0.")

        beta_tot = float(np.sqrt(float(self.beta) ** 2 + float(im_beta) ** 2))
        z = np.log(float(im_median) / float(self.theta)) / beta_tot
        return float(np.clip(norm.cdf(z), 0.0, 1.0))


@dataclass(frozen=True)
class DiscreteDamageModel(BaseDamageModel):
    """
    Discrete single-damage-state model defined by (IM, PoE) pairs.

    Parameters
    ----------
    im
        Strictly increasing IM values (N).
    poe
        PoE values in [0, 1] (N).
    log_im
        If True, interpolation is performed in ln(IM) space. In this case,
        IM must be strictly > 0.
    """

    im: Tuple[float, ...]
    poe: Tuple[float, ...]
    log_im: bool = False

    _interp: interp1d = field(init=False, repr=False)

    def __post_init__(self) -> None:
        im = _as_1d_float(self.im, "im")
        poe = _as_1d_float(self.poe, "poe")

        if im.size != poe.size:
            raise ValueError("im and poe must have the same length.")
        if im.size < 2:
            raise ValueError("Discrete models require at least 2 points.")
        if np.any(~np.isfinite(im)) or np.any(~np.isfinite(poe)):
            raise ValueError("im/poe must be finite.")
        if np.any(np.diff(im) <= 0.0):
            raise ValueError("im must be strictly increasing.")
        if np.any(poe < 0.0) or np.any(poe > 1.0):
            raise ValueError("poe values must be in [0, 1].")

        if self.log_im:
            if np.any(im <= 0.0):
                raise ValueError("im must be > 0 when log_im=True.")
            x = np.log(im)
        else:
            x = im

        y = _clip01(poe)
        fill = (float(y[0]), float(y[-1]))

        object.__setattr__(
            self,
            "_interp",
            interp1d(
                x,
                y,
                kind="linear",
                bounds_error=False,
                fill_value=fill,
                assume_sorted=True,
            ),
        )

    def poe(self, im: Any) -> np.ndarray:
        im_arr = _as_1d_float(im, "im")
        out = np.zeros_like(im_arr, dtype=float)

        msk = im_arr > 0.0
        if np.any(msk):
            x = np.log(im_arr[msk]) if self.log_im else im_arr[msk]
            out[msk] = self._interp(x)

        return _clip01(out)


@dataclass(frozen=True)
class DamageScale:
    """Damage scale definition: id + ordered levels."""

    id: str
    levels: Tuple[str, ...]

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], path: str) -> "DamageScale":
        ds_id = _require_str(data, "id", path)
        levels_raw = data.get("levels")

        if not isinstance(levels_raw, list) or not levels_raw:
            raise ValueError(f"{path}.levels must be a non-empty array.")

        levels: List[str] = []
        for i, lv in enumerate(levels_raw):
            if not isinstance(lv, str) or not lv.strip():
                raise ValueError(f"{path}.levels[{i}] must be a string.")
            levels.append(lv.strip())

        if len(set(levels)) != len(levels):
            raise ValueError(f"{path}.levels must not contain duplicates.")

        return cls(id=ds_id, levels=tuple(levels))

    def to_dict(self) -> Dict[str, Any]:
        return {"id": self.id, "levels": list(self.levels)}


@dataclass(frozen=True)
class IMBounds:
    """IM validity bounds: min >= 0, max > min."""

    min: float
    max: float

    @classmethod
    def from_dict(cls, data: Mapping[str, Any], path: str) -> "IMBounds":
        vmin = data.get("min")
        vmax = data.get("max")

        if not _is_number(vmin) or float(vmin) < 0.0:
            raise ValueError(f"{path}.min must be a number >= 0.")
        if not _is_number(vmax) or float(vmax) <= 0.0:
            raise ValueError(f"{path}.max must be a number > 0.")
        if float(vmax) <= float(vmin):
            raise ValueError(f"{path}.max must be > min.")

        return cls(min=float(vmin), max=float(vmax))

    def to_dict(self) -> Dict[str, Any]:
        return {"min": self.min, "max": self.max}


@dataclass
class FragilityModel:
    """
    Multi-damage-state fragility model for a given taxonomy and IMT.

    A FragilityModel represents one source/author model. It contains a
    damage scale (ordered levels) and a single-damage-state model for
    each level.

    Parameters
    ----------
    id
        Unique model identifier.
    taxonomy
        Taxonomy string used to associate this model to exposure assets.
    imt
        Intensity measure type (e.g., "PGA", "SA(0.3)").
    model_type
        Model family identifier (e.g., "lognormal_continuous", "discrete").
    damage_scale
        Damage scale definition (id + ordered levels).
    im_bounds
        IM validity bounds (min/max).
    damage_models
        Mapping: level_id -> BaseDamageModel.
    metadata
        Optional free-form metadata (author, year, reference, notes, ...).
    """

    id: str
    taxonomy: str
    imt: str
    model_type: str
    damage_scale: DamageScale
    im_bounds: IMBounds
    damage_models: Dict[str, BaseDamageModel]
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.id, str) or not self.id.strip():
            raise ValueError("id must be a non-empty string.")
        if not isinstance(self.taxonomy, str) or not self.taxonomy.strip():
            raise ValueError("taxonomy must be a non-empty string.")
        if not isinstance(self.imt, str) or not self.imt.strip():
            raise ValueError("imt must be a non-empty string.")
        if not isinstance(self.model_type, str) or not self.model_type.strip():
            raise ValueError("model_type must be a non-empty string.")
        if not isinstance(self.damage_models, dict) or not self.damage_models:
            raise ValueError("damage_models must be a non-empty dict.")
        if not isinstance(self.metadata, dict):
            raise ValueError("metadata must be a dict.")

        expected = set(self.damage_scale.levels)
        found = set(self.damage_models.keys())
        if found != expected:
            missing = sorted(expected - found)
            extra = sorted(found - expected)
            msg = "damage_models must match damage_scale.levels exactly."
            if missing:
                msg += f" Missing: {missing}."
            if extra:
                msg += f" Extra: {extra}."
            raise ValueError(msg)

        for key, val in self.damage_models.items():
            if not isinstance(key, str) or not key.strip():
                raise ValueError("Damage level ids must be non-empty strings.")
            if not isinstance(val, BaseDamageModel):
                raise ValueError("Each value must be a BaseDamageModel.")

    def poe(self, level: str, im: Any) -> np.ndarray:
        """Probability of exceedance for a given damage level."""
        if level not in self.damage_models:
            raise KeyError(f"Unknown damage level: {level}")
        return self.damage_models[level].poe(im)

    def poe_all(self, im: Any) -> Dict[str, np.ndarray]:
        """PoE for all damage levels (deterministic IM)."""
        out: Dict[str, np.ndarray] = {}
        for lv in self.damage_scale.levels:
            out[lv] = self.damage_models[lv].poe(im)
        return out

    def poe_lognormal_im(
        self,
        level: str,
        im_median: float,
        im_beta: float,
    ) -> float:
        """Expected PoE for a level with lognormal IM uncertainty."""
        if level not in self.damage_models:
            raise KeyError(f"Unknown damage level: {level}")
        return self.damage_models[level].poe_lognormal_im(im_median, im_beta)

    def poe_lognormal_im_all(
        self,
        im_median: float,
        im_beta: float,
    ) -> Dict[str, float]:
        """Expected PoE for all levels with lognormal IM uncertainty."""
        out: Dict[str, float] = {}
        for lv in self.damage_scale.levels:
            out[lv] = self.damage_models[lv].poe_lognormal_im(
                im_median,
                im_beta,
            )
        return out

    def poe_distribution(
        self,
        level: str,
        im_values: Sequence[float],
        probabilities: Sequence[float],
    ) -> float:
        """Expected PoE for a level with discrete IM uncertainty."""
        if level not in self.damage_models:
            raise KeyError(f"Unknown damage level: {level}")
        return self.damage_models[level].poe_distribution(
            im_values,
            probabilities,
        )

    def poe_distribution_all(
        self,
        im_values: Sequence[float],
        probabilities: Sequence[float],
    ) -> Dict[str, float]:
        """Expected PoE for all levels with discrete IM uncertainty."""
        out: Dict[str, float] = {}
        for lv in self.damage_scale.levels:
            out[lv] = self.damage_models[lv].poe_distribution(
                im_values,
                probabilities,
            )
        return out

    def state_probabilities(
        self,
        im: Any,
        tail_key: str = _DEFAULT_TAIL_KEY,
    ) -> Dict[str, np.ndarray]:
        """
        Convert exceedance PoE into mutually exclusive state probabilities.

        For ordered levels L1..Ln with exceedance PoE(Lk):
        - P(L1) = 1 - PoE(L1)
        - P(Lk) = PoE(L(k-1)) - PoE(Lk) for k=2..n
        - P(tail) = PoE(Ln)

        The output includes an additional tail key (default: "GT_LAST").
        """
        levels = self.damage_scale.levels
        if not levels:
            raise ValueError("damage_scale.levels is empty.")

        poe_stack = np.vstack(
            [self.damage_models[lv].poe(im) for lv in levels]
        )

        out: Dict[str, np.ndarray] = {}
        out[levels[0]] = _clip01(1.0 - poe_stack[0])

        for i in range(1, len(levels)):
            out[levels[i]] = _clip01(poe_stack[i - 1] - poe_stack[i])

        out[tail_key] = _clip01(poe_stack[-1])
        return out

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        path: str = "model",
    ) -> "FragilityModel":
        """Build a FragilityModel from a dictionary."""
        mid = _require_str(data, "id", path)
        taxonomy = _require_str(data, "taxonomy", path)
        imt = _require_str(data, "imt", path)
        model_type = _require_str(data, "model_type", path)

        ds = DamageScale.from_dict(
            _require_dict(data, "damage_scale", path),
            f"{path}.damage_scale",
        )
        bnd = IMBounds.from_dict(
            _require_dict(data, "im_bounds", path),
            f"{path}.im_bounds",
        )

        damage_models: Dict[str, BaseDamageModel] = {}

        if model_type == "lognormal_continuous":
            params = _require_dict(data, "parameters", path)
            for lv in ds.levels:
                p = params.get(lv)
                if not isinstance(p, dict):
                    raise ValueError(f"{path}.parameters.{lv} must be an object.")
                theta = p.get("theta")
                beta = p.get("beta")
                if not _is_number(theta) or float(theta) <= 0.0:
                    msg = f"{path}.parameters.{lv}.theta must be > 0."
                    raise ValueError(msg)
                if not _is_number(beta) or float(beta) <= 0.0:
                    msg = f"{path}.parameters.{lv}.beta must be > 0."
                    raise ValueError(msg)
                damage_models[lv] = LognormalDamageModel(
                    theta=float(theta),
                    beta=float(beta),
                )

        elif model_type == "discrete":
            tables = _require_dict(data, "tables", path)
            for lv in ds.levels:
                t = tables.get(lv)
                if not isinstance(t, dict):
                    raise ValueError(f"{path}.tables.{lv} must be an object.")
                im = t.get("im")
                poe = t.get("poe")
                log_im = t.get("log_im", False)

                if not isinstance(log_im, bool):
                    msg = f"{path}.tables.{lv}.log_im must be boolean."
                    raise ValueError(msg)
                if not isinstance(im, list) or len(im) < 2:
                    msg = f"{path}.tables.{lv}.im must be len >= 2."
                    raise ValueError(msg)
                if not isinstance(poe, list) or len(poe) < 2:
                    msg = f"{path}.tables.{lv}.poe must be len >= 2."
                    raise ValueError(msg)
                if len(im) != len(poe):
                    msg = f"{path}.tables.{lv}: im/poe size mismatch."
                    raise ValueError(msg)

                damage_models[lv] = DiscreteDamageModel(
                    im=tuple(float(v) for v in im),
                    poe=tuple(float(v) for v in poe),
                    log_im=log_im,
                )
        else:
            raise ValueError(f"{path}.model_type unsupported: {model_type}")

        reserved = {
            "id",
            "taxonomy",
            "imt",
            "model_type",
            "damage_scale",
            "im_bounds",
            "parameters",
            "tables",
        }
        meta = {k: v for k, v in dict(data).items() if k not in reserved}

        return cls(
            id=mid,
            taxonomy=taxonomy,
            imt=imt,
            model_type=model_type,
            damage_scale=ds,
            im_bounds=bnd,
            damage_models=damage_models,
            metadata=meta,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the model to a JSON-ready dictionary."""
        out: Dict[str, Any] = {
            "id": self.id,
            "taxonomy": self.taxonomy,
            "imt": self.imt,
            "model_type": self.model_type,
            "damage_scale": self.damage_scale.to_dict(),
            "im_bounds": self.im_bounds.to_dict(),
        }
        out.update(self.metadata)

        if self.model_type == "lognormal_continuous":
            params: Dict[str, Any] = {}
            for lv in self.damage_scale.levels:
                mdl = self.damage_models[lv]
                if not isinstance(mdl, LognormalDamageModel):
                    raise TypeError("Expected LognormalDamageModel.")
                params[lv] = {"theta": mdl.theta, "beta": mdl.beta}
            out["parameters"] = params

        elif self.model_type == "discrete":
            tables: Dict[str, Any] = {}
            for lv in self.damage_scale.levels:
                mdl = self.damage_models[lv]
                if not isinstance(mdl, DiscreteDamageModel):
                    raise TypeError("Expected DiscreteDamageModel.")
                tables[lv] = {
                    "im": list(mdl.im),
                    "poe": list(mdl.poe),
                    "log_im": bool(mdl.log_im),
                }
            out["tables"] = tables

        else:
            raise ValueError(f"Unsupported model_type: {self.model_type}")

        return out


@dataclass
class FragilityCollection:
    """
    Collection of multi-damage-state fragility models.

    A FragilityCollection is a container of FragilityModel instances,
    retrievable by their unique model id. It supports loading and saving
    the full database as JSON-ready dictionaries.

    Parameters
    ----------
    metadata
        Free-form collection metadata (name, date, description, ...).
    models
        Mapping model_id -> FragilityModel.
    """

    metadata: Dict[str, Any] = field(default_factory=dict)
    models: Dict[str, FragilityModel] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.metadata, dict):
            raise ValueError("metadata must be a dict.")
        if not isinstance(self.models, dict):
            raise ValueError("models must be a dict.")

        for mid, model in self.models.items():
            if not isinstance(mid, str) or not mid.strip():
                raise ValueError("Model ids must be non-empty strings.")
            if not isinstance(model, FragilityModel):
                raise ValueError("models values must be FragilityModel.")
            if model.id != mid:
                raise ValueError("models key must match model.id.")

    def add(self, model: FragilityModel) -> None:
        """Add or replace a model by id."""
        self.models[model.id] = model

    def get(self, model_id: str) -> FragilityModel:
        """Get a model by id."""
        return self.models[model_id]

    def list_ids(self) -> List[str]:
        """Return all model ids (sorted)."""
        return sorted(self.models.keys())

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "FragilityCollection":
        """Build a collection from a JSON-ready dictionary."""
        if not isinstance(data, dict):
            raise ValueError("Root JSON value must be an object.")

        tval = data.get("type")
        if tval != _SCHEMA_TYPE:
            raise ValueError(f"root.type must be '{_SCHEMA_TYPE}'.")

        sval = data.get("schema_version")
        if sval != _SCHEMA_VERSION:
            msg = f"root.schema_version must be '{_SCHEMA_VERSION}'."
            raise ValueError(msg)

        meta = data.get("metadata")
        if not isinstance(meta, dict):
            raise ValueError("root.metadata must be an object.")

        models_raw = data.get("models")
        if not isinstance(models_raw, list) or not models_raw:
            raise ValueError("root.models must be a non-empty array.")

        models: Dict[str, FragilityModel] = {}
        for i, mdata in enumerate(models_raw):
            if not isinstance(mdata, dict):
                raise ValueError(f"models[{i}] must be an object.")
            model = FragilityModel.from_dict(mdata, path=f"models[{i}]")
            if model.id in models:
                raise ValueError(f"Duplicate model id: {model.id}")
            models[model.id] = model

        return cls(metadata=dict(meta), models=models)

    @classmethod
    def from_json(cls, path: str) -> "FragilityCollection":
        """Load a collection from a JSON file."""
        with open(path, "r", encoding="utf-8") as fobj:
            data = json.load(fobj)
        return cls.from_dict(data)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize the collection to a JSON-ready dictionary."""
        models = [self.models[mid].to_dict() for mid in self.list_ids()]
        return {
            "type": _SCHEMA_TYPE,
            "schema_version": _SCHEMA_VERSION,
            "metadata": dict(self.metadata),
            "models": models,
        }

    def to_json(self, path: str, indent: int = 2) -> None:
        """Write the collection to a JSON file."""
        with open(path, "w", encoding="utf-8") as fobj:
            json.dump(self.to_dict(), fobj, indent=indent, ensure_ascii=False)
            fobj.write("\n")
