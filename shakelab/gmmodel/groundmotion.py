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
ShakeScenario: ground-motion interface and providers.

This module defines the ground-motion layer used by ShakeScenario. It provides
a stable API to compute intensity measures (IMs) at arbitrary sites for a given
scenario event, independently of any impact/fragility/exposure logic.

Main components
---------------
- ScenarioEvent:
  Minimal seismic source container (hypocentre, magnitude, optional metadata)
  providing convenience accessors (epicentre, depth_km).

- GroundMotionProvider:
  Factory/registry to instantiate a selected ground-motion backend:
  - "gmpe": analytical/statistical GMPEs implemented in ShakeLab
  - "shakemap": precomputed products (placeholder / skeleton)
  - "plugin": external or numerical simulators (placeholder / skeleton)

- GroundMotionContext:
  Lightweight context that binds a ScenarioEvent to a provider and exposes
  helpers to evaluate IM at multiple sites or at a single lon/lat location.

Output convention
-----------------
All providers return a pair (im, sigma_ln) for each site:
- im: median IM value in linear space (NOT logarithmic). Units are assumed to be
  consistent by design (no unit conversion is performed here).
- sigma_ln: standard deviation of ln(IM), dimensionless. If not available,
  providers should return 0.0.

Distance metric
---------------
For GMPE-based providers, the distance metric is taken from the GMPE itself
(via `DISTANCE_METRIC`) and is used consistently for both the GMPE input and any
diagnostic distance computations. Hypocentral distances rely on WgsPoint 3-D
distance methods (using point elevation).

Notes
-----
None
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite

from typing import (
    Any,
    Callable,
    Dict,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
)

from shakelab.libutils.geodeticN.primitives import WgsPoint


# ---------------------------------------------------------------------------
# Scenario primitives
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ScenarioEvent:
    """
    Seismic event container for impact scenarios.

    This class represents a single seismic source used in deterministic
    or scenario-based impact calculations.

    Conventions
    -----------
    - Coordinates use WGS84 longitude/latitude in degrees.
    - Elevation is expressed in meters.

    Attributes
    ----------
    magnitude
        Magnitude value, assumed compatible with the selected GMPE.
    hypocentre
        Hypocentral location as a WgsPoint instance.
    origin_time
        Optional origin time (ISO string or Date-like object).
    mechanism
        Optional focal mechanism mapping (e.g. strike/dip/rake).
    event_id
        Optional event identifier.
    magnitude_type
        Optional magnitude type (e.g. "Mw", "ML").
    """

    magnitude: float
    hypocentre: WgsPoint
    origin_time: Optional[Any] = None
    mechanism: Optional[Mapping[str, float]] = None
    event_id: Optional[str] = None
    magnitude_type: Optional[str] = None

    def __post_init__(self) -> None:
        self.validate()

    # ------------------------------------------------------------------
    # Derived geometric properties
    # ------------------------------------------------------------------

    @property
    def longitude(self) -> float:
        """Return hypocentre longitude (degrees)."""
        return float(self.hypocentre.longitude)

    @property
    def latitude(self) -> float:
        """Return hypocentre latitude (degrees)."""
        return float(self.hypocentre.latitude)

    @property
    def depth_km(self) -> float:
        """
        Return depth in kilometers (positive downward).

        Depth is derived from hypocentre elevation as:
            depth_km = -elevation / 1000.
        """
        return -float(self.hypocentre.elevation) / 1000.0

    @property
    def epicentre(self) -> WgsPoint:
        """
        Return epicentre (surface projection) as a WgsPoint.

        The returned point has the same longitude and latitude as the
        hypocentre, with elevation set to 0.0 m.
        """
        return WgsPoint(
            longitude=self.longitude,
            latitude=self.latitude,
            elevation=0.0,
        )

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_lonlat_depth(
        cls,
        magnitude: float,
        longitude: float,
        latitude: float,
        depth_km: float,
        **kwargs: Any,
    ) -> "ScenarioEvent":
        """
        Build a ScenarioEvent from longitude, latitude and depth.

        Parameters
        ----------
        magnitude
            Event magnitude.
        longitude, latitude
            Hypocentre coordinates in degrees (WGS84).
        depth_km
            Depth in kilometers (positive downward).

        Returns
        -------
        ScenarioEvent
            New instance with hypocentre elevation set to
            -depth_km * 1000 meters.
        """
        hypo = WgsPoint(
            longitude=float(longitude),
            latitude=float(latitude),
            elevation=-float(depth_km) * 1000.0,
        )
        return cls(
            magnitude=float(magnitude),
            hypocentre=hypo,
            **kwargs,
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate(self) -> None:
        """
        Validate event geometry and numerical consistency.

        The validation is intentionally lightweight and checks only:
        - numeric type and finiteness of magnitude,
        - presence and type of hypocentre,
        - finite longitude, latitude and elevation,
        - longitude/latitude ranges.
        """
        if not isinstance(self.magnitude, (int, float)) or isinstance(
            self.magnitude, bool
        ):
            raise TypeError("magnitude must be a number.")

        if not isfinite(float(self.magnitude)):
            raise ValueError("magnitude must be finite.")

        if not isinstance(self.hypocentre, WgsPoint):
            raise TypeError("hypocentre must be a WgsPoint.")

        lon = float(self.hypocentre.longitude)
        lat = float(self.hypocentre.latitude)
        elev = float(self.hypocentre.elevation)

        if not (isfinite(lon) and isfinite(lat) and isfinite(elev)):
            raise ValueError(
                "hypocentre longitude, latitude and elevation "
                "must be finite."
            )

        if not (-180.0 <= lon <= 180.0):
            raise ValueError("longitude out of range [-180, 180].")

        if not (-90.0 <= lat <= 90.0):
            raise ValueError("latitude out of range [-90, 90].")

        if self.mechanism is not None and not isinstance(
            self.mechanism, Mapping
        ):
            raise TypeError("mechanism must be a mapping or None.")


# ---------------------------------------------------------------------------
# Ground-motion providers (factory + interface)
# ---------------------------------------------------------------------------

class GroundMotionProvider:
    """
    Factory/registry for ground-motion evaluation backends.

    This class registers provider constructors (GMPE, ShakeMap, plugins, ...)
    and instantiates them from an identifier plus an optional configuration
    mapping.

    Output convention
    -----------------
    Concrete providers must implement::

        evaluate(
            imt: str,
            sites: Sequence[WgsPoint],
            event: ScenarioEvent,
            **kwargs: Any,
        ) -> Tuple[List[float], List[float]]

    and return:
    - im: median intensity-measure value (linear, NOT logarithmic)
    - sigma_ln: standard deviation of ln(IM) (dimensionless)

    No unit conversion is performed here.
    """

    _REGISTRY: Dict[str, Callable[..., "_BaseProvider"]] = {}

    @classmethod
    def register(cls, provider_id: str) -> Callable[..., Any]:
        """
        Decorator to register a provider constructor.

        Parameters
        ----------
        provider_id
            Registry identifier for the provider (e.g. "gmpe", "shakemap",
            "plugin"). It must be a non-empty string.

        Returns
        -------
        Callable
            A decorator that registers the decorated constructor.
        """
        if not isinstance(provider_id, str) or not provider_id.strip():
            raise ValueError("provider_id must be a non-empty string.")

        def _decorator(fn: Callable[..., "_BaseProvider"]) -> Callable[..., Any]:
            cls._REGISTRY[provider_id] = fn
            return fn

        return _decorator

    @classmethod
    def available_ids(cls) -> List[str]:
        """
        Return the list of registered provider IDs (sorted).
        """
        return sorted(cls._REGISTRY.keys())

    @classmethod
    def from_id(
        cls,
        provider_id: str,
        *,
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """
        Instantiate a provider by registry id.

        Parameters
        ----------
        provider_id
            The registered provider identifier.
        config
            Provider-specific configuration mapping. The mapping is copied
            into a plain dict and expanded as keyword arguments to the
            registered constructor.

        Returns
        -------
        _BaseProvider
            The instantiated provider.

        Raises
        ------
        KeyError
            If `provider_id` is not registered.
        TypeError
            If `config` is provided but is not a mapping.
        """
        if provider_id not in cls._REGISTRY:
            raise KeyError(
                f"Unknown provider_id: {provider_id!r}. Available: "
                f"{', '.join(cls.available_ids())}"
            )

        if config is None:
            cfg: Dict[str, Any] = {}
        else:
            if not isinstance(config, Mapping):
                raise TypeError("config must be a mapping (dict-like).")
            cfg = dict(config)

        return cls._REGISTRY[provider_id](**cfg)

    # ------------------------------------------------------------------
    # Convenience constructors
    # ------------------------------------------------------------------

    @classmethod
    def gmpe(
        cls,
        gmpe_name: str,
        *,
        distance_approx: str = "ellipsoid",
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """
        Convenience constructor for the GMPE provider.

        Parameters
        ----------
        gmpe_name
            Canonical GMPE name or alias, as defined by the ShakeLab GMPE
            registry (registry.json).
        distance_approx
            Approximation model forwarded to hypocentral distance calculations.
            Typical values: "ellipsoid", "sphere".
        config
            Optional extra configuration mapping forwarded to the GMPE provider
            constructor (and, ultimately, to the GMPE constructor).

        Returns
        -------
        _BaseProvider
            The instantiated GMPE provider.
        """
        cfg = dict(config or {})
        cfg["gmpe_name"] = gmpe_name
        cfg["distance_approx"] = distance_approx
        return cls.from_id("gmpe", config=cfg)

    @classmethod
    def shakemap(
        cls,
        *,
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """
        Convenience constructor for precomputed ShakeMap providers.

        Parameters
        ----------
        config
            Optional provider-specific configuration mapping.

        Returns
        -------
        _BaseProvider
            The instantiated ShakeMap provider.
        """
        return cls.from_id("shakemap", config=config)

    @classmethod
    def plugin(
        cls,
        plugin_id: str,
        *,
        config: Optional[Mapping[str, Any]] = None,
    ) -> "_BaseProvider":
        """
        Convenience constructor for plugin/simulation providers.

        Parameters
        ----------
        plugin_id
            Identifier used by the plugin provider to select the backend
            implementation.
        config
            Optional provider-specific configuration mapping.

        Returns
        -------
        _BaseProvider
            The instantiated plugin provider.
        """
        cfg = dict(config or {})
        cfg["plugin_id"] = plugin_id
        return cls.from_id("plugin", config=cfg)


class _BaseProvider:
    """
    Provider interface (internal base class).

    The concrete providers (GMPE / ShakeMap / Plugin) must implement a single
    method, :meth:`evaluate`, with a stable output convention.

    Output convention
    -----------------
    The provider must return:

        (im, sigma_ln)

    where, for each requested site:
    - im is the **median** intensity-measure value (linear, NOT logarithmic),
      expressed in **SI units** (no unit conversion is performed here).
    - sigma_ln is the standard deviation of ln(IM) (dimensionless). If the
      provider cannot supply uncertainty, it should return 0.0.

    Notes
    -----
    Many GMPEs return (mu_ln, sigma_ln); in that case return im = exp(mu_ln)
    while keeping sigma_ln unchanged.
    """

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """
        Evaluate ground motion at the given sites.

        Parameters
        ----------
        imt
            Intensity measure type (e.g. "PGA", "PGV", "SA(0.3)").
        sites
            Target sites as WgsPoint objects.
        event
            ScenarioEvent instance (updated API: `hypocentre`, `epicentre`,
            `depth_km`, etc.).
        **kwargs
            Provider-specific options.

        Returns
        -------
        (im, sigma_ln)
            im: list of median IM values (linear, SI units).
            sigma_ln: list of lognormal sigmas in ln-space (dimensionless).

        Raises
        ------
        NotImplementedError
            If the provider does not implement this method.
        """
        raise NotImplementedError


@GroundMotionProvider.register("gmpe")
class _GmpeProvider(_BaseProvider):
    """
    GMPE-backed ground motion provider.

    The GMPE implementation is selected by name (or alias) through the
    ShakeLab GMPE registry located in `shakelab.gmmodel.gmpe.registry`.

    Output convention
    -----------------
    Returns:
        (im, sigma_ln)

    where:
    - im is the median IM in linear units (NOT logarithmic)
    - sigma_ln is the standard deviation of ln(IM)

    No unit conversion is performed here. The GMPE is assumed to provide a
    consistent unit system (ideally SI) by design.

    Parameters
    ----------
    gmpe_name
        Canonical GMPE name or alias as defined in registry.json.
    gmpe_kwargs
        Extra keyword arguments forwarded to the GMPE constructor.

    Notes
    -----
    - Distance metric is taken from self._gmpe.DISTANCE_METRIC.
    - Many GMPEs return (mu_ln, sigma_ln). This provider converts the mean
      to linear median IM using:
          im = exp(mu_ln)
    """

    def __init__(
        self,
        gmpe_name: str,
        distance_approx: str = "ellipsoid",
        **gmpe_kwargs: Any,
    ) -> None:
        from shakelab.gmmodel.gmpe.registry import create_gmpe

        self._gmpe = create_gmpe(gmpe_name, **gmpe_kwargs)
        approx = str(distance_approx).strip().lower()
        self._distance_approx = approx or "ellipsoid"

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """
        Evaluate median IM and sigma_ln at the given sites.

        Parameters
        ----------
        imt
            Intensity measure type (e.g. "PGA", "PGV", "SA(0.3)").
        sites
            Target sites as WgsPoint objects.
        event
            ScenarioEvent instance (updated API).
        **kwargs
            Reserved for provider-specific runtime options (currently unused).

        Returns
        -------
        (im, sigma_ln)
            im: list of median IM values (linear).
            sigma_ln: list of lognormal sigmas in ln-space.
        """
        _ = kwargs  # reserved for future runtime options

        metric = self._get_distance_metric()
        im: List[float] = []
        sigma_ln: List[float] = []

        for site in sites:
            dist_km = self._compute_distance_km(
                event,
                site,
                metric,
                approx=self._distance_approx,
            )

            mean_ln, sig_ln = self._gmpe.ground_motion(
                imt,
                event.magnitude,
                dist_km,
            )

            mean_ln_f = float(mean_ln)
            sig_ln_f = float(sig_ln)

            if not isfinite(mean_ln_f):
                im.append(float("nan"))
            else:
                im.append(exp(mean_ln_f))

            if not isfinite(sig_ln_f) or sig_ln_f < 0.0:
                sig_ln_f = float("nan")
            sigma_ln.append(sig_ln_f)

        return im, sigma_ln

    def _get_distance_metric(self) -> str:
        """
        Return the distance metric required by the GMPE.
    
        The metric is taken from `self._gmpe.DISTANCE_METRIC` (expected to be
        always defined and standard in ShakeLab GMPEs).
        """
        metric = str(getattr(self._gmpe, "DISTANCE_METRIC")).strip().lower()
    
        # Normalize standard synonyms
        if metric in {"repi", "epicentral"}:
            return "epicentral"
    
        if metric in {"rhypo", "hypocentral"}:
            return "hypocentral"
    
        # If ShakeLab guarantees standard values, you may prefer raising here
        # to catch unexpected entries early.
        raise ValueError(
            f"Unsupported GMPE distance metric: {metric!r}. "
            "Expected 'epicentral/repi' or 'hypocentral/rhypo'."
        )

    @staticmethod
    def _compute_distance_km(
        event: ScenarioEvent,
        site: WgsPoint,
        metric: str,
        approx: str = "ellipsoid",
    ) -> float:
        """
        Compute distance in km according to the requested metric.

        Parameters
        ----------
        event
            ScenarioEvent instance.
        site
            Target site as WgsPoint.
        metric
            Distance metric identifier.
        approx
            Approximation model for hypocentral distance, forwarded to
            WgsPoint.hypocentral_distance_to(). Supported: "ellipsoid",
            "sphere".

        Returns
        -------
        float
            Distance in km.

        Notes
        -----
        geodeticN primitives operate in meters; here we convert to km.
        Hypocentral distance uses point elevations (3-D straight-line
        distance in ECEF/ellipsoid or sphere approximation).
        """
        if metric == "epicentral":
            return float(event.epicentre.epicentral_distance_to(site)) / 1000.0
        
        if metric == "hypocentral":
            return float(
                event.hypocentre.hypocentral_distance_to(
                    site, approx=approx
                )
            ) / 1000.0
        
        raise ValueError(f"Unsupported distance metric: {metric!r}")


@GroundMotionProvider.register("shakemap")
class _ShakeMapProvider(_BaseProvider):
    """
    Precomputed ShakeMap provider (skeleton).

    This provider is expected to:
    - load a ShakeMap grid/point dataset from disk or a service
    - interpolate values to requested sites
    - return **linear** median IM in SI units and sigma of ln(IM)
      (if available; else sigma=0)

    Parameters
    ----------
    path
        Path (or URL) to the ShakeMap source.
    interp
        Interpolation method identifier (provider-specific).
    **kwargs
        Provider-specific options.
    """

    def __init__(self, path: str, interp: str = "bilinear", **kwargs: Any) -> None:
        self.path = path
        self.interp = interp
        self.kwargs = dict(kwargs)

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """
        Evaluate median IM and sigma_ln at the given sites.

        Notes
        -----
        This skeleton returns zeros. A real implementation should:
        - read IM values for `imt` in SI units
        - interpolate to each requested site
        - return sigma_ln if provided by the ShakeMap product, else 0.0
        """
        im = [0.0 for _ in sites]
        sigma_ln = [0.0 for _ in sites]
        return im, sigma_ln


@GroundMotionProvider.register("plugin")
class _PluginProvider(_BaseProvider):
    """
    Plugin/simulation provider (skeleton).

    This provider is expected to dispatch to a plugin system, e.g.:
    - numerical wave propagation
    - analytical / Green's function models
    - external services

    The provider must return **linear** median IM in SI units and sigma_ln.

    Parameters
    ----------
    plugin_id
        Identifier used by the plugin system to select the backend.
    **kwargs
        Provider-specific options.
    """

    def __init__(self, plugin_id: str, **kwargs: Any) -> None:
        self.plugin_id = plugin_id
        self.kwargs = dict(kwargs)

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        event: ScenarioEvent,
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """
        Evaluate median IM and sigma_ln at the given sites.

        Notes
        -----
        This skeleton returns zeros. A real implementation should:
        - dispatch to the selected plugin backend
        - compute IM in SI units (linear)
        - provide sigma_ln if available, else 0.0
        """
        im = [0.0 for _ in sites]
        sigma_ln = [0.0 for _ in sites]
        return im, sigma_ln


# ---------------------------------------------------------------------------
# GroundMotionContext
# ---------------------------------------------------------------------------

@dataclass
class GroundMotionContext:
    """
    Bind an event with a provider and expose uniform evaluation methods.

    This is a thin convenience layer that stores a `ScenarioEvent` and a
    ground-motion provider, and offers helpers for evaluating IM at multiple
    sites or a single site.

    Output convention
    -----------------
    All methods return:
    - im: median intensity-measure value (linear, NOT logarithmic)
    - sigma_ln: standard deviation of ln(IM) (dimensionless)
    """

    event: ScenarioEvent
    provider: _BaseProvider

    def evaluate(
        self,
        imt: str,
        sites: Sequence[WgsPoint],
        **kwargs: Any,
    ) -> Tuple[List[float], List[float]]:
        """
        Evaluate IM and sigma_ln for an IMT at multiple sites.

        Parameters
        ----------
        imt
            Intensity measure type (e.g. "PGA", "PGV", "SA(0.3)").
        sites
            Sequence of target sites as WgsPoint objects.
        **kwargs
            Provider-specific options forwarded to the provider.

        Returns
        -------
        (im, sigma_ln)
            im: list of median IM values (linear).
            sigma_ln: list of lognormal sigmas in ln-space.
        """
        return self.provider.evaluate(imt, sites, self.event, **kwargs)

    def evaluate_at_site(
        self,
        imt: str,
        lon: float,
        lat: float,
        elevation_m: float = 0.0,
        **kwargs: Any,
    ) -> Tuple[float, float]:
        """
        Evaluate IM and sigma_ln for an IMT at a single WGS84 site.

        Parameters
        ----------
        imt
            Intensity measure type (e.g. "PGA", "PGV", "SA(0.3)").
        lon, lat
            Site longitude and latitude (WGS84).
        elevation_m
            Site elevation in meters.
        **kwargs
            Provider-specific options forwarded to the provider.

        Returns
        -------
        (im, sigma_ln)
            im: median IM (linear).
            sigma_ln: sigma of ln(IM).
        """
        site = WgsPoint(
            longitude=float(lon),
            latitude=float(lat),
            elevation=float(elevation_m),
        )
        im, sigma_ln = self.evaluate(imt, [site], **kwargs)
        return float(im[0]), float(sigma_ln[0])

