# ****************************************************************************
# Copyright (C) 2019-2026, ShakeLab Developers.
# This file is part of ShakeLab.
#
# ShakeLab is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ShakeLab is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
# ****************************************************************************
"""
Reading XML documents with the standard library.

The module collects the XML utilities shared by the format readers of
ShakeLab (StationXML, QuakeML): reading a document from any source and
converting an element tree into plain Python containers.

Design choices
--------------
- Based only on :mod:`xml.etree.ElementTree`.
- Namespaces are removed from tags and attribute names, since ShakeLab
  readers address elements by their local names
  (``"{http://www.fdsn.org/xml/station/1}Channel"`` becomes
  ``"Channel"``).
- Text is never converted implicitly: codes such as ``"00"`` or
  ``"0123"`` must stay strings. Readers convert values explicitly.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from typing import Any, Optional

__all__ = [
    "read_xml",
    "strip_namespaces",
    "to_dict",
    "find_text",
    "find_float",
]


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------


def read_xml(source: Any) -> ET.Element:
    """
    Read an XML document and return its root element.

    Parameters
    ----------
    source : str, path-like, bytes, file object or Element
        File path (``str`` or path-like), XML text (``str`` starting
        with ``<``, with or without XML declaration), ``bytes``, an
        open file object, or an existing element (returned unchanged).

    Returns
    -------
    xml.etree.ElementTree.Element
        Root element, with namespaces removed.

    Raises
    ------
    FileNotFoundError
        If a path does not exist.
    TypeError
        If the source type is not supported.
    xml.etree.ElementTree.ParseError
        If the document is not well-formed XML.
    """
    if isinstance(source, ET.Element):
        return source

    if isinstance(source, (bytes, bytearray)):
        root = ET.fromstring(bytes(source))
    elif isinstance(source, str) and source.lstrip().startswith("<"):
        root = ET.fromstring(source)
    elif isinstance(source, (str, os.PathLike)):
        if not os.path.isfile(source):
            raise FileNotFoundError(f"File not found: {source}")
        root = ET.parse(source).getroot()
    elif hasattr(source, "read"):
        root = ET.parse(source).getroot()
    else:
        raise TypeError(
            "source must be a path, XML text, bytes or a file object"
        )

    return strip_namespaces(root)


def strip_namespaces(element: ET.Element) -> ET.Element:
    """
    Remove namespaces from tags and attribute names, in place.

    Parameters
    ----------
    element : xml.etree.ElementTree.Element
        Root of the tree to process.

    Returns
    -------
    xml.etree.ElementTree.Element
        The same element.
    """
    for node in element.iter():
        if isinstance(node.tag, str) and "}" in node.tag:
            node.tag = node.tag.split("}", 1)[1]
        if any("}" in key for key in node.attrib):
            node.attrib = {
                key.split("}", 1)[-1]: value
                for key, value in node.attrib.items()
            }
    return element


# ----------------------------------------------------------------------------
# Conversion
# ----------------------------------------------------------------------------


def to_dict(element: ET.Element) -> dict:
    """
    Convert an element tree into nested dictionaries.

    Parameters
    ----------
    element : xml.etree.ElementTree.Element
        Element to convert.

    Returns
    -------
    dict
        ``{"tag": str, "attrib": dict, "text": str or None,
        "children": list}``, where ``text`` is the stripped text of the
        element (None if empty) and ``children`` holds the converted
        child elements in document order. Text is kept as a string.
    """
    text = element.text.strip() if element.text else ""
    return {
        "tag": element.tag,
        "attrib": dict(element.attrib),
        "text": text or None,
        "children": [to_dict(child) for child in element],
    }


def find_text(
    element: Optional[ET.Element],
    path: str,
    default: Optional[str] = None,
) -> Optional[str]:
    """
    Return the stripped text of a sub-element.

    Parameters
    ----------
    element : xml.etree.ElementTree.Element or None
        Element to search (None gives ``default``).
    path : str
        ElementTree path of the sub-element.
    default : str, optional
        Value returned if the sub-element or its text are missing.

    Returns
    -------
    str or None
        Stripped text, or ``default``.
    """
    if element is None:
        return default
    text = element.findtext(path)
    if text is None or not text.strip():
        return default
    return text.strip()


def find_float(
    element: Optional[ET.Element],
    path: str,
    default: Optional[float] = None,
) -> Optional[float]:
    """
    Return the text of a sub-element as a float.

    Parameters
    ----------
    element : xml.etree.ElementTree.Element or None
        Element to search (None gives ``default``).
    path : str
        ElementTree path of the sub-element.
    default : float, optional
        Value returned if the sub-element or its text are missing.

    Returns
    -------
    float or None
        Value of the text, or ``default``.

    Raises
    ------
    ValueError
        If the text is not a number.
    """
    text = find_text(element, path)
    return default if text is None else float(text)
