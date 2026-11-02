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
Unit tests for :mod:`shakelab.libutils.xml`.

The tests cover reading documents from every supported source, the
removal of namespaces, the conversion to dictionaries (text kept as
strings) and the helpers for text and numeric values.
"""

from __future__ import annotations

import io
import xml.etree.ElementTree as ET

import pytest

from shakelab.libutils.xml import (
    find_float,
    find_text,
    read_xml,
    strip_namespaces,
    to_dict,
)

DOC = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<Root xmlns="http://example.org/ns" xmlns:x="http://example.org/x">\n'
    '  <Item code="00" x:flag="yes">\n'
    '    <Value> 1.5 </Value>\n'
    '    <Name>first line\n    second line</Name>\n'
    '  </Item>\n'
    '  <Empty/>\n'
    '</Root>\n'
)


def _check(root):
    """The document is read with local names and intact text."""
    assert root.tag == "Root"
    item = root.find("Item")
    assert item.attrib == {"code": "00", "flag": "yes"}
    assert find_float(item, "Value") == 1.5
    assert "first line" in item.findtext("Name")
    assert "second line" in item.findtext("Name")


def test_read_from_every_source(tmp_path):
    """Paths, text, bytes and file objects give the same tree."""
    path = tmp_path / "doc.xml"
    path.write_text(DOC, encoding="utf-8")

    _check(read_xml(str(path)))
    _check(read_xml(path))
    _check(read_xml(DOC))
    _check(read_xml(DOC.split("\n", 1)[1]))
    _check(read_xml(DOC.encode("utf-8")))
    _check(read_xml(io.BytesIO(DOC.encode("utf-8"))))

    root = read_xml(DOC)
    assert read_xml(root) is root


def test_read_errors(tmp_path):
    """Missing files, wrong types and malformed XML are rejected."""
    with pytest.raises(FileNotFoundError):
        read_xml(str(tmp_path / "missing.xml"))

    with pytest.raises(TypeError):
        read_xml(12)

    with pytest.raises(ET.ParseError):
        read_xml("<Root><Open></Root>")


def test_strip_namespaces_in_place():
    """Namespaces are removed from tags and attributes."""
    root = ET.fromstring('<a:R xmlns:a="u"><a:C a:k="1"/></a:R>')
    assert strip_namespaces(root) is root
    assert root.tag == "R"
    assert root.find("C").attrib == {"k": "1"}


def test_to_dict_keeps_strings():
    """Text is stripped but not converted; children keep their order."""
    result = to_dict(read_xml(DOC))

    assert result["tag"] == "Root"
    assert [child["tag"] for child in result["children"]] == [
        "Item", "Empty"
    ]
    item = result["children"][0]
    assert item["attrib"]["code"] == "00"
    assert item["children"][0] == {
        "tag": "Value", "attrib": {}, "text": "1.5", "children": []
    }
    assert result["children"][1]["text"] is None


def test_find_helpers():
    """Missing elements give the default; invalid numbers raise."""
    root = read_xml("<R><A> x </A><B/><C>abc</C></R>")

    assert find_text(root, "A") == "x"
    assert find_text(root, "B") is None
    assert find_text(root, "Z", "d") == "d"
    assert find_text(None, "A", "d") == "d"
    assert find_float(root, "Z", 2.0) == 2.0

    with pytest.raises(ValueError):
        find_float(root, "C")
