# ****************************************************************************
#
# Copyright (C) 2019-2020, ShakeLab Developers.
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

from typing import Dict, Iterator, List, Optional

from xml.etree import ElementTree as ET

from shakelab.engineering.fragility_old import (
    FragilityModelParametric,
    FragilityModelDiscrete,
)


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
    Stream NRML 0.5 exposure assets as plain dicts.

    This iterator is memory-efficient and yields one asset at a time.

    Parameters
    ----------
    xml_path : str
        Path to the NRML exposure XML file.
    recover : bool, default False
        If True, performs a minimal, non-destructive trimming of any
        garbage that may precede the first '<' in the file.

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
    """
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
        if elem.tag != asset_tag:
            continue

        at = elem.attrib

        loc_el = elem.find(loc_tag)
        location = {
            "lat": (
                _safe_float(loc_el.attrib.get("lat")) if loc_el else None
            ),
            "lon": (
                _safe_float(loc_el.attrib.get("lon")) if loc_el else None
            ),
        }

        occs: Dict[str, int] = {}
        occs_el = elem.find(occs_tag)
        if occs_el is not None:
            for o in occs_el.findall(occ_tag):
                period = o.attrib.get("period")
                if period:
                    occs[period] = _safe_int(o.attrib.get("occupants"))

        costs: Dict[str, float] = {}
        costs_el = elem.find(costs_tag)
        if costs_el is not None:
            for c in costs_el.findall(cost_tag):
                ctype = c.attrib.get("type")
                if ctype:
                    costs[ctype] = _safe_float(c.attrib.get("value"))

        tags_el = elem.find(tags_tag)
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

        # Free memory as we stream.
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
