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
# ****************************************************************************
"""
Regression tests for the pure-Python MiniSEED reader and writer.
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from shakelab.signals.libio.mseed import (
    ENCODING_ASCII,
    ENCODING_FLOAT32,
    ENCODING_INT16,
    ENCODING_INT32,
    ENCODING_STEIM1,
    ENCODING_STEIM2,
    MSRecord,
    _decode_steim,
    _encode_steim1,
    _encode_steim2,
    _prepare_steim_data,
    msrawread,
    msrawwrite,
)


class TestMiniSeed(unittest.TestCase):
    """Test MiniSEED encoding, decoding and file round trips."""

    def setUp(self):
        """Create a temporary directory for output files."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)

    def tearDown(self):
        """Remove all temporary files."""
        self.temp_dir.cleanup()

    @staticmethod
    def _make_record(data, encoding):
        """
        Create a minimal valid MiniSEED record.

        Parameters
        ----------
        data : array-like or str
            Data assigned to the record.
        encoding : int
            MiniSEED encoding code.

        Returns
        -------
        MSRecord
            Record ready to be written.
        """
        record = MSRecord()

        record.header.update({
            "SEQUENCE_NUMBER": "000001",
            "DATA_HEADER_QUALITY_INDICATOR": "D",
            "RESERVED_BYTE": " ",
            "STATION_CODE": "TEST ",
            "LOCATION_IDENTIFIER": "00",
            "CHANNEL_IDENTIFIER": "BHZ",
            "NETWORK_CODE": "XX",
            "YEAR": 2026,
            "DAY": 200,
            "HOURS": 12,
            "MINUTES": 0,
            "SECONDS": 0,
            "UNUSED": 0,
            "MSECONDS": 0,
            "NUMBER_OF_SAMPLES": len(data),
            "SAMPLE_RATE_FACTOR": 100,
            "SAMPLE_RATE_MULTIPLIER": 1,
            "ACTIVITY_FLAGS": 0,
            "IO_FLAGS": 0,
            "DATA_QUALITY_FLAGS": 0,
            "NUMBER_OF_BLOCKETTES_TO_FOLLOW": 1,
            "TIME_CORRECTION": 0,
            "OFFSET_TO_BEGINNING_OF_DATA": 56,
            "OFFSET_TO_BEGINNING_OF_BLOCKETTE": 48,
        })

        record.blockette = {
            1000: {
                "OFFSET_NEXT": 0,
                "ENCODING_FORMAT": encoding,
                "WORD_ORDER": 1,
                "DATA_RECORD_LENGTH": 12,
                "RESERVED": 0,
            }
        }

        if isinstance(data, str):
            record.data = data
        else:
            record.data = np.asarray(data)

        return record

    def _round_trip(
        self,
        data,
        encoding,
        record_length=4096,
    ):
        """
        Write and read one record and return the decoded data.

        Parameters
        ----------
        data : array-like or str
            Input record data.
        encoding : int
            MiniSEED encoding code.
        record_length : int, optional
            Output record length in bytes.

        Returns
        -------
        numpy.ndarray or str
            Data read from the generated file.
        """
        record = self._make_record(data, encoding)

        output_file = self.output_dir / (
            "encoding_{0}.mseed".format(encoding)
        )

        msrawwrite(
            [record],
            output_file,
            record_length=record_length,
            encoding=encoding,
        )

        restored = msrawread(output_file)

        self.assertEqual(len(restored), 1)
        self.assertEqual(restored[0].nsamp, len(data))

        return restored[0].data

    def test_ascii_round_trip(self):
        """Test ASCII encoding round trip."""
        original = "MiniSEED regression test"

        restored = self._round_trip(
            original,
            ENCODING_ASCII,
        )

        self.assertEqual(restored, original)

    def test_int16_round_trip(self):
        """Test signed 16-bit integer round trip."""
        original = np.array(
            [
                -32768,
                -1000,
                -1,
                0,
                1,
                1000,
                32767,
            ],
            dtype=np.int16,
        )

        restored = self._round_trip(
            original,
            ENCODING_INT16,
        )

        np.testing.assert_array_equal(restored, original)

    def test_int32_round_trip(self):
        """Test signed 32-bit integer round trip."""
        original = np.array(
            [
                -2000000000,
                -1000000,
                -1,
                0,
                1,
                1000000,
                2000000000,
            ],
            dtype=np.int32,
        )

        restored = self._round_trip(
            original,
            ENCODING_INT32,
        )

        np.testing.assert_array_equal(restored, original)

    def test_float32_round_trip(self):
        """Test IEEE float32 round trip."""
        time = np.arange(200, dtype=np.float32) / 100.0

        original = (
            np.exp(-0.5 * time)
            * np.sin(2.0 * np.pi * 5.0 * time)
        ).astype(np.float32)

        restored = self._round_trip(
            original,
            ENCODING_FLOAT32,
        )

        np.testing.assert_array_equal(restored, original)

    def test_steim1_payload_round_trip(self):
        """Test STEIM1 payload encoding and decoding."""
        original = self._integer_signal()

        payload, sample_count, frame_count = _encode_steim1(
            original,
            frame_count=20,
        )

        restored = _decode_steim(
            payload,
            encoding=ENCODING_STEIM1,
            nsamp=sample_count,
        )

        self.assertGreater(frame_count, 0)
        self.assertEqual(sample_count, len(original))
        np.testing.assert_array_equal(restored, original)

    def test_steim2_payload_round_trip(self):
        """Test STEIM2 payload encoding and decoding."""
        original = self._integer_signal()

        payload, sample_count, frame_count = _encode_steim2(
            original,
            frame_count=20,
        )

        restored = _decode_steim(
            payload,
            encoding=ENCODING_STEIM2,
            nsamp=sample_count,
        )

        self.assertGreater(frame_count, 0)
        self.assertEqual(sample_count, len(original))
        np.testing.assert_array_equal(restored, original)

    def test_steim1_file_round_trip(self):
        """Test complete STEIM1 file round trip."""
        original = self._integer_signal()

        restored = self._round_trip(
            original,
            ENCODING_STEIM1,
            record_length=4096,
        )

        np.testing.assert_array_equal(restored, original)

    def test_steim2_file_round_trip(self):
        """Test complete STEIM2 file round trip."""
        original = self._integer_signal()

        restored = self._round_trip(
            original,
            ENCODING_STEIM2,
            record_length=4096,
        )

        np.testing.assert_array_equal(restored, original)

    def test_supported_record_lengths(self):
        """Test representative MiniSEED record lengths."""
        original = np.arange(32, dtype=np.int32)

        for record_length in (256, 512, 1024, 2048, 4096):
            with self.subTest(record_length=record_length):
                restored = self._round_trip(
                    original,
                    ENCODING_INT32,
                    record_length=record_length,
                )

                np.testing.assert_array_equal(
                    restored,
                    original,
                )

    def test_steim_rejects_float_data(self):
        """Test that STEIM encoding rejects floating-point data."""
        data = np.array(
            [0.0, 1.5, 2.0],
            dtype=np.float32,
        )

        with self.assertRaises(ValueError):
            _prepare_steim_data(data)

    def test_steim2_rejects_large_difference(self):
        """Test the signed 30-bit STEIM2 difference limit."""
        data = np.array(
            [0, 536870912],
            dtype=np.int32,
        )

        with self.assertRaises(ValueError):
            _encode_steim2(
                data,
                frame_count=1,
            )

    def test_record_too_small(self):
        """Test that silent record truncation is prevented."""
        original = np.arange(1000, dtype=np.int32)
        record = self._make_record(
            original,
            ENCODING_INT32,
        )

        output_file = self.output_dir / "too_small.mseed"

        with self.assertRaises(ValueError):
            msrawwrite(
                [record],
                output_file,
                record_length=256,
                encoding=ENCODING_INT32,
            )

    @staticmethod
    def _integer_signal():
        """
        Return a deterministic integer waveform.

        The signal combines a damped harmonic component with a linear
        trend, producing differences suitable for several STEIM packing
        modes.

        Returns
        -------
        numpy.ndarray
            Synthetic signed 32-bit signal.
        """
        time = np.arange(500, dtype=np.float64) / 100.0

        signal = (
            20000.0
            * np.exp(-0.3 * time)
            * np.sin(2.0 * np.pi * 3.0 * time)
        )

        trend = 10.0 * np.arange(500)

        return np.rint(signal + trend).astype(np.int32)


if __name__ == "__main__":
    unittest.main()