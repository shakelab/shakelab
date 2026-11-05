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
Unit tests for :mod:`shakelab.signals.io`.

They cover the sources of :func:`reader` and :func:`metadata_reader`
(files, patterns, directories, lists), the format table, the options
of the format readers and writers, the outputs of :func:`writer` (file,
directory, name template), the reading of channel metadata and their
assignment while reading waveforms, and the pickle format. Real data
come from ``tests/signals/data`` (see the README in that folder).
"""

from __future__ import annotations

import os
import pickle
import shutil
import warnings
from pathlib import Path

import numpy as np
import pytest

from shakelab.signals import io
from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, Stream, StreamCollection
from shakelab.signals.codes import parse_sid_pattern
from shakelab.signals.io import (
    IoWarning,
    formats,
    metadata_reader,
    reader,
    writer,
)
from shakelab.signals.libio.dyna import is_dyna
from shakelab.signals.libio.mseed import is_mseed
from shakelab.signals.libio.sac import is_sac
from shakelab.signals.libio.simpleascii import is_simpleascii
from shakelab.signals.libio.tdms import is_tdms
from shakelab.signals.metadata import MetadataCollection, MetadataWarning

DATA = Path(__file__).parent / "data"
MSEED = str(DATA / "MODE.mseed")
SXML = str(DATA / "MODE.xml")
OTHER = str(DATA / "TT_T01.xml")


# ----------------------------------------------------------------------------
# metadata_reader
# ----------------------------------------------------------------------------


def test_metadata_reader_single_file():
    """One file gives the same content as the collection constructor."""
    meta = metadata_reader(SXML)

    assert isinstance(meta, MetadataCollection)
    assert sorted(meta.sid) == sorted(MetadataCollection(SXML).sid)
    assert len(meta) == 3


def test_metadata_reader_list():
    """The files of a list are joined in one collection."""
    meta = metadata_reader([SXML, Path(OTHER)])

    assert len(meta) == 3 + 30
    assert len(meta.stations()) == 2


def test_metadata_reader_pattern_and_directory(tmp_path):
    """Wildcard patterns and directories are expanded as in reader."""
    for path in (SXML, OTHER):
        shutil.copy(path, tmp_path)

    by_pattern = metadata_reader(str(tmp_path / "*.xml"))
    by_directory = metadata_reader(tmp_path)

    assert len(by_pattern) == len(by_directory) == 33


def test_metadata_reader_repeated_file():
    """Epochs equal to those already read are ignored, silently."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        meta = metadata_reader([SXML, SXML])

    assert len(meta) == 3


def test_metadata_reader_format(tmp_path):
    """The format can be given explicitly for unusual names."""
    path = tmp_path / "mode.station"
    shutil.copy(SXML, path)

    with pytest.raises(ValueError, match="Unknown metadata format"):
        metadata_reader(path)

    assert len(metadata_reader(path, format="xml")) == 3


def test_metadata_reader_type_error():
    """Sources other than paths are rejected."""
    with pytest.raises(TypeError):
        metadata_reader(42)

    with pytest.raises(TypeError):
        metadata_reader(b"<FDSNStationXML/>")


# ----------------------------------------------------------------------------
# reader(..., metadata=...)
# ----------------------------------------------------------------------------


def test_reader_metadata_from_file():
    """Records read with a metadata file carry their channel epoch."""
    collection = reader(MSEED, metadata=SXML)
    meta = MetadataCollection(SXML)

    for stream in collection:
        record = stream[0]
        assert not record.head.meta.empty
        assert record.head.meta == meta.get(record.sid, record.time)

    assert collection["IV.MODE..HNE"][0].orientation.azimuth == 90.0


def test_reader_metadata_same_as_set_metadata():
    """Reading with metadata equals reading and then assigning them."""
    meta = MetadataCollection(SXML)

    direct = reader(MSEED, metadata=meta)
    later = reader(MSEED)
    later.set_metadata(meta)

    raw = direct["IV.MODE..HNZ"][0].copy()
    raw.deconvolve_response()
    expected = later["IV.MODE..HNZ"][0].copy()
    expected.deconvolve_response()

    np.testing.assert_allclose(raw.data, expected.data)


def test_reader_metadata_list_of_files():
    """A list of metadata files is read with metadata_reader."""
    collection = reader(MSEED, metadata=[SXML, OTHER])

    assert all(not stream[0].head.meta.empty for stream in collection)


def test_reader_missing_channels_one_warning():
    """Channels without metadata are reported in a single warning."""
    with pytest.warns(MetadataWarning) as record:
        collection = reader(MSEED, metadata=OTHER)

    messages = [
        str(item.message) for item in record
        if issubclass(item.category, MetadataWarning)
    ]
    assert len(messages) == 1
    assert "3 channel(s)" in messages[0]
    for cha in ("HNE", "HNN", "HNZ"):
        assert f"IV.MODE..{cha}" in messages[0]

    assert all(stream[0].head.meta.empty for stream in collection)


def test_reader_existing_collection_untouched():
    """Records already in the collection keep their metadata."""
    collection = StreamCollection()
    collection.append(
        Record(data=[1.0, 2.0], delta=0.01, time="2026-01-01T00:00:00",
               sid="XX.STA..HHZ")
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", MetadataWarning)
        returned = reader(MSEED, stream_collection=collection,
                          metadata=SXML)

    assert returned is collection
    assert collection["XX.STA..HHZ"][0].head.meta.empty
    assert not collection["IV.MODE..HNZ"][0].head.meta.empty


def test_collection_read_metadata():
    """StreamCollection.read forwards the metadata to reader."""
    collection = StreamCollection()
    returned = collection.read(MSEED, metadata=SXML)

    assert returned is collection
    assert collection["IV.MODE..HNN"][0].orientation.azimuth == 0.0


def test_reader_metadata_checked_first(tmp_path):
    """Invalid metadata are reported before reading the data."""
    missing = str(tmp_path / "missing.mseed")

    with pytest.raises(TypeError, match="metadata must be"):
        reader(missing, metadata=42)

    with pytest.raises(FileNotFoundError):
        reader(MSEED, metadata=str(tmp_path / "missing.xml"))


# ----------------------------------------------------------------------------
# pickle
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["data.pkl", "data.pickle"])
def test_pickle_round_trip(tmp_path, name):
    """Collections keep data and metadata through a pickle file."""
    collection = reader(MSEED, metadata=SXML)
    path = str(tmp_path / name)
    writer(collection, path)

    again = reader(path)
    for stream in collection:
        rec, other = stream[0], again[stream.sid][0]
        assert other.time == rec.time
        assert other.delta == rec.delta
        np.testing.assert_array_equal(other.data, rec.data)
        assert other.location == rec.location
        assert other.head.meta.response is not None


def test_pickle_record_and_existing_file(tmp_path):
    """A Record is written; files are not overwritten unless asked."""
    rec = reader(MSEED)["IV.MODE..HNZ"][0]
    path = str(tmp_path / "rec.pkl")
    writer(rec, path)

    with pytest.raises(FileExistsError):
        writer(rec, path)
    writer(rec, path, owrite=True)
    again = reader(path, format="pickle")
    assert list(again)[0].sid == "IV.MODE..HNZ"


def test_pickle_of_other_objects(tmp_path):
    """Pickle files of other objects are rejected."""
    path = tmp_path / "other.pkl"
    path.write_bytes(pickle.dumps({"a": 1}))

    with pytest.raises(TypeError, match="not a Record"):
        reader(str(path))


# ----------------------------------------------------------------------------
# Sources
# ----------------------------------------------------------------------------


def _collection(nsamp=10):
    """Three channels of one station, one record each."""
    collection = StreamCollection()
    for cha in ("HHE", "HHN", "HHZ"):
        collection.append(
            Record(data=np.arange(float(nsamp)), delta=0.01,
                   time="2026-01-01T00:00:00", sid=f"XX.AA..{cha}",
                   units="m/s**2")
        )
    return collection


def test_reader_directory_skips_unknown_files(tmp_path):
    """Unknown files of a directory are skipped with one warning."""
    shutil.copy(MSEED, tmp_path)
    (tmp_path / "notes.txt").write_text("not a waveform")
    (tmp_path / "readme").write_text("not a waveform")
    (tmp_path / ".hidden.txt").write_text("ignored")
    (tmp_path / "sub").mkdir()

    with pytest.warns(IoWarning, match="2 file") as caught:
        collection = reader(tmp_path)

    assert len(collection) == 3
    assert "notes.txt" in str(caught[0].message)
    assert ".hidden" not in str(caught[0].message)
    assert caught[0].filename == __file__


def test_reader_named_file_of_unknown_format(tmp_path):
    """A file given by name must have a recognizable format."""
    path = tmp_path / "notes"
    path.write_text("not a waveform")

    with pytest.raises(ValueError, match="set format"):
        reader(path)

    # Pickle files are never recognized from their content.
    pickled = tmp_path / "saved"
    writer(reader(MSEED), str(tmp_path / "saved.pkl"))
    (tmp_path / "saved.pkl").rename(pickled)
    with pytest.raises(ValueError, match="set format"):
        reader(pickled)
    assert len(reader(pickled, format="pkl")) == 3


# ----------------------------------------------------------------------------
# Formats recognized from the content
# ----------------------------------------------------------------------------


LIBIO = Path(__file__).parent / "libio" / "data"

DETECTORS = {
    "mseed": is_mseed,
    "sac": is_sac,
    "dyna": is_dyna,
    "simpleascii": is_simpleascii,
    "tdms": is_tdms,
}


def _expected_format(path):
    """Format of a test data file, or None for other files."""
    if path.suffix in (".json", ".npz", ".xml", ".md", ".txt") and (
        "simpleascii" not in path.parts
    ):
        return None
    if path.name.startswith("LICENSE"):
        return None
    if "libmseed" in path.parts or "MSEED" in path.parts:
        return "mseed"
    if path.name == "MODE.mseed":
        return "mseed"
    if "sac" in path.parts:
        return None if path.name == "testxy.sac" else "sac"   # text SAC
    for name, folder in (("simpleascii", "simpleascii"), ("tdms", "tdms"),
                         ("dyna", "ASCII")):
        if folder in path.parts:
            return name
    return None


SAMPLES = sorted(
    path for folder in (DATA, LIBIO) for path in folder.rglob("*")
    if path.is_file()
)


@pytest.mark.parametrize(
    "path", SAMPLES, ids=[path.name for path in SAMPLES]
)
def test_detectors_recognize_only_their_format(path):
    """Each test file is recognized by its own format only."""
    head = path.read_bytes()[:1024]
    found = [name for name, detect in DETECTORS.items() if detect(head)]

    expected = _expected_format(path)
    assert found == ([expected] if expected else [])


def test_detectors_short_or_wrong_input():
    """Short heads are not recognized; text is not accepted."""
    for detect in DETECTORS.values():
        assert not detect(b"")
        assert not detect(b"\x00" * 16)
        with pytest.raises(TypeError):
            detect("TIMESERIES")
    assert formats()["mseed"]["detect"]
    assert not formats()["pickle"]["detect"]


def test_reader_by_content(tmp_path):
    """Files without a known extension are read after their content."""
    sds = tmp_path / "2012" / "IV" / "MODE" / "HNZ.D"
    sds.mkdir(parents=True)
    shutil.copy(MSEED, sds / "IV.MODE..HNZ.D.2012.153")
    shutil.copy(LIBIO / "simpleascii" / "int32.slist.txt", sds / "slist")

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        collection = reader(sds)
    assert len(collection) == 3 + 1

    dyna = tmp_path / "acc"
    shutil.copy(DATA / "ESM" / "ASCII" /
                "IV.MODE..HNE.D.IT-2012-0008.ACC.MP.ASC", dyna)
    assert reader(dyna).sid == ["IV.MODE..HNE"]

    simple = tmp_path / "waveform.asc"    # .asc, but not DYNA
    shutil.copy(LIBIO / "simpleascii" / "int32.slist.txt", simple)
    assert len(reader(simple)) == 1

    pattern = reader(str(LIBIO / "mseed" / "libmseed" / "*-int32.mseed2"))
    assert len(pattern) == 1


def test_reader_missing_files(tmp_path):
    """Missing files and patterns without matches are errors."""
    with pytest.raises(FileNotFoundError, match="No such file"):
        reader(tmp_path / "missing.mseed")

    with pytest.raises(FileNotFoundError, match="No file matches"):
        reader(str(tmp_path / "*.mseed"))

    with pytest.raises(ValueError, match="source or list_file"):
        reader()

    with pytest.raises(TypeError, match="source must be"):
        reader(42)


def test_reader_list_file(tmp_path):
    """List files hold paths, patterns and comments; both sources add."""
    folder = tmp_path / "data"
    folder.mkdir()
    shutil.copy(MSEED, folder / "a.mseed")
    writer(_collection(), str(folder / "b.mseed"))
    writer(Record(data=[1.0, 2.0], delta=0.5, time="2026-01-01T00:00:00",
                  sid="XX.BB..HHZ"), str(tmp_path / "c.sac"))

    listing = tmp_path / "files.txt"
    listing.write_text("# waveforms\n\ndata/a.mseed\n  data/b*  \n")

    by_list = reader(list_file=listing)
    assert sorted(by_list.sid) == [
        "IV.MODE..HNE", "IV.MODE..HNN", "IV.MODE..HNZ",
        "XX.AA..HHE", "XX.AA..HHN", "XX.AA..HHZ",
    ]

    both = reader(str(tmp_path / "c.sac"), list_file=str(listing))
    assert len(both) == 7
    assert "XX.BB..HHZ" in both.sid


def test_reader_options_checked():
    """Options must be accepted by the reader of the format."""
    with pytest.raises(TypeError, match="unknown option"):
        reader(MSEED, mseed_backend="python", encoding=11)

    with pytest.raises(ValueError, match="mseed_backend"):
        reader(MSEED, mseed_backend=True)

    collection = reader(MSEED, mseed_backend="python", byte_order="be")
    assert len(collection) == 3


def test_format_names():
    """Formats are listed, aliased, rejected or not supported yet."""
    table = formats()
    assert table["sac"]["per_record"]
    assert table["mseed"]["write"] and not table["tdms"]["write"]
    assert ".asc" in table["dyna"]["extensions"]

    with pytest.raises(ValueError, match="unknown format"):
        reader(MSEED, format="xyz")
    with pytest.raises(NotImplementedError):
        reader(MSEED, format="seisan")
    with pytest.raises(TypeError):
        reader(MSEED, format=5)


def test_metadata_reader_skips_unknown_files(tmp_path):
    """Metadata directories may hold other files."""
    shutil.copy(SXML, tmp_path)
    (tmp_path / "notes.txt").write_text("not metadata")

    with pytest.warns(IoWarning, match="notes.txt"):
        meta = metadata_reader(tmp_path)

    assert len(meta) == 3


def test_metadata_warning_points_to_the_caller():
    """Metadata warnings while reading point to the user's code."""
    with pytest.warns(MetadataWarning) as caught:
        reader(MSEED, metadata=OTHER)
    assert caught[0].filename == __file__

    with pytest.warns(MetadataWarning) as caught:
        StreamCollection().read(MSEED, metadata=OTHER)
    assert caught[0].filename == __file__


# ----------------------------------------------------------------------------
# writer
# ----------------------------------------------------------------------------


def test_writer_one_record_per_file(tmp_path):
    """SAC and DYNA files hold one record; a file name gets one only."""
    collection = _collection()
    record = collection[0][0]

    writer(record, str(tmp_path / "rec.sac"))
    assert (tmp_path / "rec.sac").is_file()
    writer(record, str(tmp_path / "rec.asc"))
    assert reader(str(tmp_path / "rec.asc"))[0][0].nsamp == 10

    with pytest.raises(ValueError, match="one record per file"):
        writer(collection, str(tmp_path / "all.sac"))
    assert not (tmp_path / "all.sac").exists()


def test_writer_directory(tmp_path):
    """A directory gets one file per record or per channel."""
    collection = _collection()

    writer(collection, str(tmp_path / "sac") + "/", format="sac")
    assert sorted(path.name for path in (tmp_path / "sac").iterdir()) == [
        f"XX.AA..{cha}.20260101T000000.sac" for cha in ("HHE", "HHN", "HHZ")
    ]

    writer(collection, tmp_path / "sac", format="mseed")
    assert (tmp_path / "sac" / "XX.AA..HHZ.mseed").is_file()

    with pytest.raises(ValueError, match="set format"):
        writer(collection, tmp_path / "sac")


def test_writer_template(tmp_path):
    """Templates name the files and create the folders."""
    collection = _collection()

    writer(collection, str(tmp_path / "{sta}" / "{cha}.{start}.sac"))
    assert len(list((tmp_path / "AA").glob("HH?.20260101T000000.sac"))) == 3

    writer(collection, str(tmp_path / "{net}.{sta}.mseed"))
    again = reader(str(tmp_path / "XX.AA.mseed"))
    assert again.sid == collection.sid

    with pytest.raises(ValueError, match="same name"):
        writer(collection, str(tmp_path / "{sta}.sac"))

    with pytest.raises(ValueError, match="fields"):
        writer(collection, str(tmp_path / "{month}.mseed"))


def test_writer_existing_files(tmp_path):
    """No file is written if one of them exists, unless owrite."""
    collection = _collection()
    template = str(tmp_path / "{cha}.mseed")
    (tmp_path / "HHN.mseed").write_bytes(b"")

    with pytest.raises(FileExistsError):
        writer(collection, template)
    assert not (tmp_path / "HHE.mseed").exists()

    writer(collection, template, owrite=True)
    assert len(reader(str(tmp_path / "HHN.mseed"))) == 1

    with pytest.raises(FileExistsError):
        writer(collection, str(tmp_path / "HHN.mseed"), mseed_backend="python")


def test_writer_errors(tmp_path):
    """Invalid data, targets, formats and options are rejected."""
    collection = _collection()

    with pytest.raises(TypeError, match="data must be"):
        writer([1, 2], str(tmp_path / "a.mseed"))
    with pytest.raises(TypeError, match="target"):
        writer(collection, 3)
    with pytest.raises(ValueError, match="set format"):
        writer(collection, str(tmp_path / "a.xyz"))
    with pytest.raises(NotImplementedError, match="writing"):
        writer(collection, str(tmp_path / "a.tdms"))
    with pytest.raises(TypeError, match="unknown option"):
        writer(collection[0][0], str(tmp_path / "a.sac"), reclen=512)
    assert not list(tmp_path.iterdir())


def test_writer_template_time_fields(tmp_path):
    """Templates take the fields of the SeedLink archives."""
    collection = _collection()

    writer(collection, str(tmp_path / "{year}" / "{jday}" /
                           "{net}.{sta}.{loc}.{cha}.mseed"))
    assert sorted(path.name for path in (tmp_path / "2026" / "001").iterdir()
                  ) == [f"XX.AA..{cha}.mseed" for cha in ("HHE", "HHN",
                                                          "HHZ")]

    writer(collection, str(tmp_path / "{time:%Y-%m-%d}.{sid}.sac"))
    assert (tmp_path / "2026-01-01.XX.AA..HHZ.sac").is_file()


# ----------------------------------------------------------------------------
# Selection
# ----------------------------------------------------------------------------


def test_reader_select():
    """Channels are selected with the patterns of signals.codes."""
    assert reader(MSEED, select="IV.MODE..HN[EN]").sid == [
        "IV.MODE..HNE", "IV.MODE..HNN"
    ]
    assert reader(MSEED, select=["IV.MODE..HNZ", "XX"]).sid == [
        "IV.MODE..HNZ"
    ]
    assert len(reader(MSEED, select="IV.{MODE,MIRA}")) == 3
    assert len(reader(MSEED, select="IV.MIRA")) == 0

    with pytest.raises(ValueError):
        reader(MSEED, select="IV..")
    with pytest.raises(TypeError, match="select"):
        reader(MSEED, select=5)


def test_reader_time_window():
    """Records are cut to the window, bounds included."""
    full = reader(MSEED)["IV.MODE..HNZ"][0]
    start = full.time + 10.0
    end = full.time + 20.0

    cut = reader(MSEED, starttime=start, endtime=end)["IV.MODE..HNZ"][0]
    assert cut.time == start
    assert cut.nsamp == round(10.0 / full.delta) + 1
    np.testing.assert_array_equal(
        cut.data, full.data[1000:1000 + cut.nsamp]
    )

    assert reader(MSEED, starttime=start.iso8601).sid == reader(MSEED).sid
    assert len(reader(MSEED, endtime=full.time - 1.0)) == 0

    with pytest.raises(ValueError, match="precedes"):
        reader(MSEED, starttime=end, endtime=start)
    with pytest.raises(TypeError, match="ISO 8601"):
        reader(MSEED, starttime=10.0)


def test_reader_selection_with_collection_and_metadata():
    """The metadata go to the selected records only; others untouched."""
    collection = StreamCollection()
    collection.append(Record(data=[1.0, 2.0], delta=0.01,
                             time="2026-01-01T00:00:00", sid="XX.STA..HHZ"))

    with pytest.warns(MetadataWarning) as caught:
        reader(MSEED, select="IV.MODE..HNZ", metadata=OTHER,
               stream_collection=collection)
    messages = [str(item.message) for item in caught]
    assert len(messages) == 1 and "1 channel(s)" in messages[0]
    assert collection.sid == ["IV.MODE..HNZ", "XX.STA..HHZ"]
    assert collection["XX.STA..HHZ"][0].nsamp == 2


# ----------------------------------------------------------------------------
# SDS archives
# ----------------------------------------------------------------------------


SDS = "{year}/{net}/{sta}/{cha}.D/{net}.{sta}.{loc}.{cha}.D.{year}.{jday}"


@pytest.fixture
def archive(tmp_path):
    """
    SDS archive of two stations at 1 Hz: per channel, one record per
    day from 2026-05-02 to 2026-05-04 (00:00:00 to 23:58:59) and one
    record from 23:59:00 to 00:00:59, stored in the file of the day it
    begins.
    """
    root = tmp_path / "sds"
    collection = StreamCollection()
    for sta in ("AAA", "BBB"):
        for cha in ("HHZ", "HHN"):
            sid = f"XX.{sta}..{cha}"
            for day in ("02", "03", "04"):
                collection.append(Record(
                    data=np.arange(86340.0), delta=1.0,
                    time=f"2026-05-{day}T00:00:00", sid=sid))
                collection.append(Record(
                    data=np.arange(120.0) + 1e6, delta=1.0,
                    time=f"2026-05-{day}T23:59:00", sid=sid))
    # Records kept separate (Stream.append would join them).
    records = StreamCollection()
    for record in collection.records():
        records.stream.append(Stream(sid=record.sid))
        records.stream[-1].record.append(record)
    writer(records, str(root / SDS), format="mseed",
           mseed_backend="python")
    return root


def test_sds_reads_only_the_files_needed(archive):
    """Only the requested channels and days are read."""
    files = sorted(
        path.relative_to(archive).as_posix()
        for path in archive.rglob("*") if path.is_file()
    )
    assert len(files) == 2 * 2 * 3

    # The day before the window and the day of the window.
    needed = io._sds_files(
        str(archive), parse_sid_pattern("XX.AAA..HHZ"),
        Date("2026-05-03T10:00:00"), Date("2026-05-03T10:00:09"),
    )
    assert [os.path.basename(path) for path in needed] == [
        "XX.AAA..HHZ.D.2026.122", "XX.AAA..HHZ.D.2026.123"
    ]

    collection = reader(archive, layout="sds", select="XX.AAA..HHZ",
                        starttime="2026-05-03T10:00:00",
                        endtime="2026-05-03T10:00:09")
    assert collection.sid == ["XX.AAA..HHZ"]
    record = collection[0][0]
    assert record.nsamp == 10
    assert record.data[0] == 36000.0


def test_sds_previous_day(archive):
    """Data just after midnight are in the file of the day before."""
    collection = reader(archive, layout="sds", select="XX.*..HHN",
                        starttime="2026-05-04T00:00:10",
                        endtime="2026-05-04T00:00:19")
    assert collection.sid == ["XX.AAA..HHN", "XX.BBB..HHN"]
    for stream in collection:
        assert stream[0].nsamp == 10
        assert stream[0].data[0] == 1e6 + 70


def test_sds_whole_days_and_patterns(archive):
    """Patterns select stations and channels; missing days are fine."""
    collection = reader(archive, layout="sds", select="XX.{AAA,BBB}",
                        starttime="2026-05-01T00:00:00",
                        endtime="2026-05-07T00:00:00")
    assert len(collection) == 4
    for stream in collection:
        assert sum(record.nsamp for record in stream) == 3 * 86460


def test_sds_errors(archive, tmp_path):
    """Invalid requests, missing archives and empty results."""
    with pytest.raises(ValueError, match="requires starttime and endtime"):
        reader(archive, layout="sds", starttime="2026-05-03T00:00:00")
    with pytest.raises(ValueError, match="unknown layout"):
        reader(archive, layout="bud", starttime="2026-05-03T00:00:00",
               endtime="2026-05-03T01:00:00")
    with pytest.raises(ValueError, match="list_file"):
        reader(archive, layout="sds", list_file="files.txt",
               starttime="2026-05-03T00:00:00",
               endtime="2026-05-03T01:00:00")
    with pytest.raises(FileNotFoundError, match="archive"):
        reader(tmp_path / "missing", layout="sds",
               starttime="2026-05-03T00:00:00",
               endtime="2026-05-03T01:00:00")

    with pytest.warns(IoWarning, match="no data"):
        collection = reader(archive, layout="sds", select="YY",
                            starttime="2026-05-03T00:00:00",
                            endtime="2026-05-03T01:00:00")
    assert len(collection) == 0
