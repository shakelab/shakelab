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
Utilities for reading/writing OpenQuake NRML fragments used in ShakeLab.

This module provides:
- Conversion of in-memory fragility collections to NRML 0.5 XML.
- Generation of NRML 0.5 exposure models from project data.
- Streaming parse of NRML 0.5 exposure models into plain Python dicts.

Notes
-----
* The writer functions emit NRML 0.5 with a default namespace for NRML and
  the 'gml' prefix for GML.
* The reader handles default NRML namespaces explicitly (qualified tags).
"""

from __future__ import annotations

from typing import Any, Dict, Iterator, List, Optional
from xml.etree import ElementTree as ET
from typing import Dict, List, Mapping, Optional, Tuple, Union
from datetime import date

from shakelab.engineering.exposure.exposure import (
    Typology,
    Asset,
    Exposure,
)
#from shakelab.engineering.fragility_old import (
#    FragilityModelParametric,
#    FragilityModelDiscrete,
#)

Number = Union[int, float]

_NRML_NS = "http://openquake.org/xmlns/nrml/0.5"
_GML_NS = "http://www.opengis.net/gml"


# --------------------------------------------------------------------------- #
# XML helpers
# --------------------------------------------------------------------------- #

def indent(elem: ET.Element, level: int = 0) -> None:
    """
    Pretty-print an XML Element in-place with indentation.

    Parameters
    ----------
    elem : xml.etree.ElementTree.Element
        Root element to indent.
    level : int, optional
        Initial indentation level.
    """
    i = "\n" + level * "    "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "    "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for child in list(elem):
            indent(child, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def _safe_float(x: Optional[str], default: float = 0.0) -> float:
    """
    Convert a string to float, returning a default value on failure.

    Parameters
    ----------
    x : str or None
        Input string.
    default : float
        Fallback value if conversion fails.

    Returns
    -------
    float
    """
    try:
        return float(x)  # type: ignore[arg-type]
    except Exception:
        return default


def _safe_int(x: Optional[str], default: int = 0) -> int:
    """
    Convert a string to int, returning a default value on failure.

    The conversion is tolerant to inputs that are floats in string form.

    Parameters
    ----------
    x : str or None
        Input string.
    default : int
        Fallback value if conversion fails.

    Returns
    -------
    int
    """
    try:
        return int(float(x))  # type: ignore[arg-type]
    except Exception:
        return default


def _iterparse_clean_stream(xml_path: str) -> Iterator[ET.Element]:
    """
    Yield 'end' events from a minimally recovered XML stream.

    The only recovery applied is trimming any leading bytes before the
    first '<' character. No per-line filtering is performed, so tags
    split across multiple lines remain intact.

    Parameters
    ----------
    xml_path : str
        Path to the XML file.

    Yields
    ------
    xml.etree.ElementTree.Element
        Elements as they complete ('end' events).
    """
    import io

    with open(xml_path, "r", encoding="utf-8-sig") as f:
        data = f.read()

    i = data.find("<")
    if i > 0:
        data = data[i:]

    stream = io.StringIO(data)
    context = ET.iterparse(stream, events=("end",))
    for _event, elem in context:
        yield elem


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #

def fragility_to_xml(
    fragility_collection,
    xml_file: str,
    ndl: float = 0.1,
) -> None:
    """
    Write a fragility collection to NRML 0.5 XML.

    Parameters
    ----------
    fragility_collection
        Collection exposing an iterable 'model' of fragility models.
        Supported types: FragilityModelParametric, FragilityModelDiscrete.
    xml_file : str
        Output XML path.
    ndl : float, default 0.1
        No-damage limit for the 'imls' element.
    """
    nrml = ET.Element(
        "nrml",
        {"xmlns": _NRML_NS, "xmlns:gml": _GML_NS},
    )

    fm = ET.SubElement(
        nrml,
        "fragilityModel",
        {
            "assetCategory": "buildings",
            "id": "buildings",
            "lossCategory": "structural",
        },
    )

    ET.SubElement(fm, "description").text = "Fragility model"
    ET.SubElement(fm, "limitStates").text = "D1 D2 D3 D4 D5"

    for m in fragility_collection.model:
        if isinstance(m, FragilityModelParametric):
            ff = ET.SubElement(
                fm,
                "fragilityFunction",
                {"id": f"{m.id}", "format": "continuous", "shape": "lognormal"},
            )
            ET.SubElement(
                ff,
                "imls",
                {
                    "imt": f"{m.gmt}",
                    "noDamageLimit": f"{ndl}",
                    "minIML": f"{m.bounds[0]}",
                    "maxIML": f"{m.bounds[1]}",
                },
            )
            for dsl, (mean, stdv) in m.damage_state.items():
                ET.SubElement(
                    ff,
                    "params",
                    {"ls": f"{dsl}", "mean": f"{mean}", "stddev": f"{stdv}"},
                )

        elif isinstance(m, FragilityModelDiscrete):
            ff = ET.SubElement(
                fm,
                "fragilityFunction",
                {"id": f"{m.id}", "format": "discrete"},
            )
            im = ET.SubElement(
                ff,
                "imls",
                {"imt": f"{m.gmt}", "noDamageLimit": f"{ndl}"},
            )
            im.text = " ".join(str(n) for n in m.gmi)
            for dsl, poes in m.damage_state.items():
                poe = ET.SubElement(ff, "poes", {"ls": f"{dsl}"})
                poe.text = " ".join(str(n) for n in poes)

    indent(nrml)
    ET.ElementTree(nrml).write(
        xml_file, encoding="utf-8", xml_declaration=True
    )


def exposure_to_xml(exposure, taxonomy_tree, xml_file: str) -> None:
    """
    Write an exposure model (case-specific) to NRML 0.5 XML.

    WARNING
    -------
    This function is tailored to the Friuli case study and expects the
    'exposure' object to expose attributes like:
    - exposure.location (iterable)
    - for each location 'li': li.id, li.code, li.area, li.latitude,
      li.longitude, and li.taxonomy
    - each taxonomy entry exposes: id, number_of_buildings, occupants,
      and is resolved through 'taxonomy_tree.get_element(...)' which
      provides a '.branch' dict {branch_id: weight}.

    Parameters
    ----------
    exposure
        Project-specific data structure (see warning above).
    taxonomy_tree
        Structure exposing 'get_element(tax_id) -> element', with
        'element.branch' mapping to branch weights.
    xml_file : str
        Output XML path.
    """
    nrml = ET.Element("nrml", {"xmlns": _NRML_NS, "xmlns:gml": _GML_NS})

    em = ET.SubElement(
        nrml,
        "exposureModel",
        {"id": "buildings", "category": "buildings", "taxonomySource": "GEM"},
    )
    ET.SubElement(em, "description").text = "Buildings"

    con = ET.SubElement(em, "conversions")
    ctp = ET.SubElement(con, "costTypes")
    for name in (
        "structural",
        "nonstructural",
        "contents",
        "business_interruption",
    ):
        ET.SubElement(
            ctp,
            "costType",
            {"name": name, "unit": "EUR", "type": "per_asset"},
        )
    ET.SubElement(con, "area", {"type": "per_asset", "unit": "SQM"})
    ET.SubElement(em, "tagNames").text = "sezione comune"

    ass = ET.SubElement(em, "assets")
    for li in exposure.location:
        for tax in li.taxonomy:
            tte = taxonomy_tree.get_element(tax.id)
            for bri, brw in tte.branch.items():
                nob = int(tax.number_of_buildings * brw)
                if nob <= 0:
                    continue

                asset_id = "_".join([li.id, tax.id, bri])
                ast = ET.SubElement(
                    ass,
                    "asset",
                    {
                        "id": asset_id,
                        "name": li.code,
                        "area": f"{li.area}",
                        "number": f"{nob}",
                        "taxonomy": bri,
                    },
                )
                ET.SubElement(
                    ast,
                    "location",
                    {"lon": f"{li.longitude}", "lat": f"{li.latitude}"},
                )

                occ = ET.SubElement(ast, "occupancies")
                for period in ("day", "transit", "night"):
                    occupants = tax.occupants.get(period, 0)
                    ET.SubElement(
                        occ,
                        "occupancy",
                        {
                            "occupants": f"{int(occupants)}",
                            "period": period,
                        },
                    )

                cst = ET.SubElement(ast, "costs")
                for ctype in (
                    "structural",
                    "nonstructural",
                    "contents",
                    "business_interruption",
                ):
                    ET.SubElement(cst, "cost", {"type": ctype, "value": "0"})

                ET.SubElement(
                    ast,
                    "tags",
                    {"sezione": "all_included", "comune": li.code},
                )

    indent(nrml)
    ET.ElementTree(nrml).write(
        xml_file, encoding="utf-8", xml_declaration=True
    )


# --------------------------------------------------------------------------- #
# Readers
# --------------------------------------------------------------------------- #

def iter_exposure_assets(
    xml_path: str,
    recover: bool = False,
) -> Iterator[Dict]:
    """
    Stream NRML 0.5 exposure assets as plain dicts, with robust
    namespace handling. The iterator is memory-efficient and yields
    one asset at a time.

    This reader first tries exact lookups with the NRML namespace
    and then falls back to localname-only matches so it works even
    when the XML uses default namespaces or different prefixes.

    Parameters
    ----------
    xml_path : str
        Path to the NRML exposure XML file.
    recover : bool, default False
        If True, uses a tolerant stream that trims any garbage bytes
        preceding the first '<'. This is non-destructive and only
        affects the parsing stream.

    Yields
    ------
    dict
        A compact record for each <asset>, e.g.:
        {
          "id": str, "name": str, "taxonomy": str,
          "number": int, "area": float,
          "location": {"lat": float|None, "lon": float|None},
          "occupancies": {"day": int, "night": int, ...},
          "costs": {"structural": float, "nonstructural": float, ...},
          "tags": {"comune": str, "sezione": str, ...}
        }

    Notes
    -----
    * 'location' attributes are cast to float; missing values are None.
    * 'occupancies' values are summed if multiple entries per period
      exist in the source.
    * 'costs' values are summed if multiple entries per type exist.
    * Memory is released by clearing the processed <asset> element.
    """

    # ------- helpers (local to avoid polluting the module scope) ---------

    def _tag_localname(tag: str) -> str:
        """Return local part of a tag, stripping '{ns}' if present."""
        if not isinstance(tag, str):
            return ""
        if tag.startswith("{"):
            return tag.split("}", 1)[1]
        return tag

    def _find_any_ns(parent, localname: str):
        """Find first direct child by localname, ignoring namespace."""
        for child in list(parent):
            if _tag_localname(child.tag) == localname:
                return child
        return None

    def _findall_any_ns(parent, localname: str):
        """Yield direct children by localname, ignoring namespace."""
        for child in list(parent):
            if _tag_localname(child.tag) == localname:
                yield child

    # -------------------- parsing stream selection -----------------------

    events = (
        _iterparse_clean_stream(xml_path)
        if recover
        else (elem for _ev, elem in ET.iterparse(xml_path, events=("end",)))
    )

    n = _NRML_NS
    asset_tag = f"{{{n}}}asset"
    loc_tag = f"{{{n}}}location"
    occs_tag = f"{{{n}}}occupancies"
    occ_tag = f"{{{n}}}occupancy"
    costs_tag = f"{{{n}}}costs"
    cost_tag = f"{{{n}}}cost"
    tags_tag = f"{{{n}}}tags"

    for elem in events:
        is_asset = elem.tag == asset_tag or _tag_localname(elem.tag) == "asset"
        if not is_asset:
            continue

        at = elem.attrib

        # ---- location (try exact-NS then fallback by localname) ----
        loc_el = elem.find(loc_tag)
        if loc_el is None:
            loc_el = _find_any_ns(elem, "location")
        location = {
            "lat": (
                _safe_float(loc_el.attrib.get("lat"))
                if loc_el is not None else None
            ),
            "lon": (
                _safe_float(loc_el.attrib.get("lon"))
                if loc_el is not None else None
            ),
        }

        # ---- occupancies (sum by period) ----
        occs_el = elem.find(occs_tag)
        if occs_el is None:
            occs_el = _find_any_ns(elem, "occupancies")
        occs: Dict[str, int] = {}
        if occs_el is not None:
            occ_elems = list(occs_el.findall(occ_tag))
            if not occ_elems:
                occ_elems = list(_findall_any_ns(occs_el, "occupancy"))
            for o in occ_elems:
                period = o.attrib.get("period")
                if not period:
                    continue
                val = _safe_int(o.attrib.get("occupants"), 0)
                occs[period] = occs.get(period, 0) + val

        # ---- costs (sum by type) ----
        costs_el = elem.find(costs_tag)
        if costs_el is None:
            costs_el = _find_any_ns(elem, "costs")
        costs: Dict[str, float] = {}
        if costs_el is not None:
            cost_elems = list(costs_el.findall(cost_tag))
            if not cost_elems:
                cost_elems = list(_findall_any_ns(costs_el, "cost"))
            for c in cost_elems:
                ctype = c.attrib.get("type")
                if not ctype:
                    continue
                cval = _safe_float(c.attrib.get("value"), 0.0)
                costs[ctype] = costs.get(ctype, 0.0) + cval

        # ---- tags (attributes of <tags/>) ----
        tags_el = elem.find(tags_tag)
        if tags_el is None:
            tags_el = _find_any_ns(elem, "tags")
        tags = dict(tags_el.attrib) if tags_el is not None else {}

        yield {
            "id": at.get("id"),
            "name": at.get("name"),
            "taxonomy": at.get("taxonomy"),
            "number": _safe_int(at.get("number"), 1),
            "area": _safe_float(at.get("area"), 0.0),
            "location": location,
            "occupancies": occs,
            "costs": costs,
            "tags": tags,
        }

        # Free memory as we stream the tree.
        elem.clear()


def load_oq_exposure(xml_path: str, recover: bool = False) -> List[Dict]:
    """
    Load all NRML 0.5 exposure assets into memory.

    NOTE
    ----
    Indexing has been removed to keep the function lightweight and
    predictable. If you need indices for lookups (e.g., by 'id' or by
    'comune' in tags), build them externally to suit your use case.

    Parameters
    ----------
    xml_path : str
        Path to the NRML exposure XML file.
    recover : bool, default False
        If True, performs a minimal, non-destructive trimming of any
        garbage that may precede the first '<' in the file.

    Returns
    -------
    list of dict
        List of asset records as produced by `iter_exposure_assets`.
    """
    return [rec for rec in iter_exposure_assets(xml_path, recover=recover)]


# --------------------------------------------------------------------------- #
# Conversion to Shakelab exposure
# --------------------------------------------------------------------------- #

def convert_oq_exposure_strict(
    records: List[Dict[str, Any]],
    *,
    group_by: str = "none",
    metadata: Optional[Dict[str, Any]] = None,
) -> "Exposure":
    """
    Convert OpenQuake exposure records to a ShakeLab Exposure object.

    Strict policy:
    - No inference from taxonomy (no by-prefix mapping).
    - Do not fabricate missing attributes.
    - Do not generate geometry (Asset.geometry is left as None).
    - Use reference_location only (mandatory in ShakeLab).

    Parameters
    ----------
    records
        List of OpenQuake exposure records as dictionaries, typically
        returned by `load_oq_exposure()`.
    group_by
        Grouping strategy:
        - "none": one ShakeLab asset per OpenQuake record
        - "id": group by record["id"]
        - "name": group by record["name"]
        - "tag:<key>": group by record["tags"][<key>]
    metadata
        Metadata dictionary for the root Exposure object. Must include at
        least `name` and `date` for ShakeLab validation.

    Returns
    -------
    Exposure
        A validated ShakeLab Exposure object.

    Raises
    ------
    ValueError
        On missing mandatory fields (e.g., cannot compute reference point).
    """
    # Local import to avoid hard dependency at module import time.
    from shakelab.engineering.exposure.exposure import (  # noqa: WPS433
        Asset,
        Exposure,
        Typology,
    )

    def _is_number(x: Any) -> bool:
        return isinstance(x, (int, float)) and not isinstance(x, bool)

    def _parse_int_or_none(x: Any) -> Optional[int]:
        if x is None:
            return None
        if isinstance(x, bool):
            return None
        if isinstance(x, int):
            return x
        if isinstance(x, float):
            return int(x) if x.is_integer() else None
        if isinstance(x, str):
            s = x.strip()
            if not s:
                return None
            try:
                return int(s)
            except ValueError:
                return None
        return None

    def _parse_float_or_none(x: Any) -> Optional[float]:
        if x is None:
            return None
        if isinstance(x, bool):
            return None
        if isinstance(x, (int, float)):
            return float(x)
        if isinstance(x, str):
            s = x.strip()
            if not s:
                return None
            try:
                return float(s)
            except ValueError:
                return None
        return None

    def _group_key(rec: Dict[str, Any]) -> str:
        if group_by == "none":
            rid = rec.get("id")
            if isinstance(rid, str) and rid.strip():
                return rid
            return str(id(rec))

        if group_by == "id":
            rid = rec.get("id")
            if isinstance(rid, str) and rid.strip():
                return rid
            return str(id(rec))

        if group_by == "name":
            name = rec.get("name")
            if isinstance(name, str) and name.strip():
                return name
            rid = rec.get("id")
            if isinstance(rid, str) and rid.strip():
                return rid
            return str(id(rec))

        if group_by.startswith("tag:"):
            key = group_by.split(":", 1)[1].strip()
            tags = rec.get("tags") or {}
            if isinstance(tags, dict):
                val = tags.get(key)
                if val is not None:
                    sval = str(val).strip()
                    if sval:
                        return sval
            rid = rec.get("id")
            if isinstance(rid, str) and rid.strip():
                return rid
            return str(id(rec))

        raise ValueError(f"Unsupported group_by='{group_by}'")

    def _pick_lon_lat(group: List[Dict[str, Any]]) -> tuple[float, float]:
        for r in group:
            loc = r.get("location") or {}
            if not isinstance(loc, dict):
                continue
            lon = loc.get("lon")
            lat = loc.get("lat")
            if _is_number(lon) and _is_number(lat):
                return float(lon), float(lat)
        raise ValueError(
            "Cannot build reference_location: missing lon/lat in group."
        )

    def _pick_occupants(rec: Dict[str, Any]) -> Optional[Dict[str, float]]:
        occ = rec.get("occupancies") or {}
        if not isinstance(occ, dict) or not occ:
            return None

        day = occ.get("day")
        night = occ.get("night")
        if _is_number(day) and _is_number(night):
            if float(day) < 0 or float(night) < 0:
                return None
            return {"day": float(day), "night": float(night)}

        return None

    def _pick_replacement_cost(rec: Dict[str, Any]) -> Optional[float]:
        costs = rec.get("costs") or {}
        if not isinstance(costs, dict) or not costs:
            return None

        repl = costs.get("replacement")
        if _is_number(repl) and float(repl) >= 0:
            return float(repl)

        numeric = [v for v in costs.values() if _is_number(v) and float(v) >= 0]
        if len(numeric) == 1:
            return float(numeric[0])

        return None

    def _pick_period(rec: Dict[str, Any]) -> Optional[Dict[str, Optional[int]]]:
        tags = rec.get("tags") or {}
        if not isinstance(tags, dict) or not tags:
            return None

        start = _parse_int_or_none(tags.get("period_start"))
        end = _parse_int_or_none(tags.get("period_end"))

        if start is None and end is None:
            return None

        return {"start": start, "end": end}

    def _pick_code_level(rec: Dict[str, Any]) -> Optional[str]:
        tags = rec.get("tags") or {}
        if not isinstance(tags, dict):
            return None
        val = tags.get("code_level")
        if val is None:
            return None
        s = str(val).strip()
        if not s:
            return None
        # Keep only values accepted by the ShakeLab model.
        allowed = {"none", "low", "moderate", "high", "retrofit"}
        return s if s in allowed else None

    def _pick_str_tag(rec: Dict[str, Any], key: str) -> Optional[str]:
        tags = rec.get("tags") or {}
        if not isinstance(tags, dict):
            return None
        val = tags.get(key)
        if val is None:
            return None
        s = str(val).strip()
        return s if s else None

    def _pick_stories(rec: Dict[str, Any]) -> Optional[int]:
        tags = rec.get("tags") or {}
        if not isinstance(tags, dict):
            return None
        return _parse_int_or_none(tags.get("stories"))

    if metadata is None:
        metadata = {}

    groups: Dict[str, List[Dict[str, Any]]] = {}
    for r in records:
        if not isinstance(r, dict):
            continue
        key = _group_key(r)
        groups.setdefault(key, []).append(r)

    assets: List[Asset] = []
    for gid, group in groups.items():
        lon, lat = _pick_lon_lat(group)

        # aggregated is informational now (geometry is optional).
        # Make it True when grouping creates multi-record assets, or when an
        # OQ record represents multiple buildings (number > 1).
        aggregated = False
        if len(group) > 1:
            aggregated = True
        else:
            number = group[0].get("number")
            if isinstance(number, int) and number > 1:
                aggregated = True

        typologies: List[Typology] = []
        for r in group:
            taxonomy = r.get("taxonomy")
            if not isinstance(taxonomy, str) or not taxonomy.strip():
                continue

            number = r.get("number")
            if not isinstance(number, int) or number < 1:
                continue

            typologies.append(
                Typology(
                    taxonomy=taxonomy.strip(),
                    count=number,
                    usage=_pick_str_tag(r, "usage"),
                    building_type=_pick_str_tag(r, "building_type"),
                    code_level=_pick_code_level(r),
                    occupants=_pick_occupants(r),
                    period=_pick_period(r),
                    replacement_cost=_pick_replacement_cost(r),
                    stories=_pick_stories(r),
                    damage_state=_pick_str_tag(r, "damage_state"),
                )
            )

        if not typologies:
            continue

        assets.append(
            Asset(
                id=str(gid),
                name=str(gid),
                aggregated=aggregated,
                aggregation_area=_parse_float_or_none(
                    group[0].get("area") if len(group) == 1 else None
                ),
                critical=False,
                geometry=None,
                reference_location={"longitude": lon, "latitude": lat},
                reference_geology={},
                typologies=typologies,
            )
        )

    exposure = Exposure(
        type="ShakeLabExposure",
        schema_version="1.0.0",
        metadata=metadata,
        assets=assets,
    )

    exposure.validate()
    return exposure
