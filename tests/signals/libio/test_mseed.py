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
from unittest.mock import patch

import numpy as np

from shakelab.libutils.timeN import Date
from shakelab.signals.base import Record, StreamCollection

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
    _encode_steim_payload,
    _prepare_steim_data,
    _sampling_rate_factors,
    msread,
    msrawread,
    msrawwrite,
    mswrite,
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

    @staticmethod
    def _make_shakelab_record(data, delta=0.01):
        """
        Create a ShakeLab waveform record for high-level writer tests.
        """
        record = Record()

        record.head.sid = "XX.TEST.00.BHZ"
        record.head.time = Date("2026-200T12:00:00.0000")
        record.head.delta = delta
        record.data = np.asarray(data)

        return record

    def test_mswrite_steim2_segmented_round_trip(self):
        """Test high-level STEIM2 writing across multiple records."""
        original = self._integer_signal()

        collection = StreamCollection()
        collection.append(
            self._make_shakelab_record(original)
        )

        output_file = self.output_dir / "segmented_steim2.mseed"

        mswrite(
            output_file,
            collection,
            encoding=ENCODING_STEIM2,
            reclen=256,
        )

        raw_records = msrawread(output_file)

        self.assertGreater(len(raw_records), 1)

        restored = np.concatenate([
            record.data
            for record in raw_records
        ])

        np.testing.assert_array_equal(restored, original)

    def test_mswrite_int32_segmented_round_trip(self):
        """Test high-level INT32 writing across multiple records."""
        original = np.arange(
            -500,
            500,
            dtype=np.int32,
        )

        collection = StreamCollection()
        collection.append(
            self._make_shakelab_record(original)
        )

        output_file = self.output_dir / "segmented_int32.mseed"

        mswrite(
            output_file,
            collection,
            encoding=ENCODING_INT32,
            reclen=256,
        )

        raw_records = msrawread(output_file)

        self.assertGreater(len(raw_records), 1)

        restored = np.concatenate([
            record.data
            for record in raw_records
        ])

        np.testing.assert_array_equal(restored, original)

    def test_mswrite_record_metadata(self):
        """Test sequence numbers and segmented record start times."""
        original = np.arange(
            1000,
            dtype=np.int32,
        )

        delta = 0.01
        record = self._make_shakelab_record(
            original,
            delta=delta,
        )

        collection = StreamCollection()
        collection.append(record)

        output_file = self.output_dir / "metadata.mseed"

        mswrite(
            output_file,
            collection,
            encoding=ENCODING_INT32,
            reclen=256,
        )

        raw_records = msrawread(output_file)

        sample_index = 0

        for index, raw_record in enumerate(
            raw_records,
            start=1,
        ):
            self.assertEqual(raw_record.seqn, index)

            expected_time = (
                record.head.time
                + sample_index * delta
            )

            self.assertAlmostEqual(
                raw_record.time - expected_time,
                0.0,
                places=4,
            )

            sample_index += raw_record.nsamp

        self.assertEqual(sample_index, len(original))

    def test_sampling_rate_factors_slow_rate(self):
        """Test representation of a 0.00875 Hz sampling rate."""
        sampling_rate = 0.00875
        delta = 1.0 / sampling_rate

        factor, multiplier = _sampling_rate_factors(delta)

        self.assertEqual(factor, 7)
        self.assertEqual(multiplier, -800)

        reconstructed_rate = factor / abs(multiplier)

        self.assertAlmostEqual(
            reconstructed_rate,
            sampling_rate,
        )

    def test_mswrite_rejects_mseed3(self):
        """Test rejection of unsupported MiniSEED 3 output."""
        collection = StreamCollection()
        collection.append(
            self._make_shakelab_record(
                np.arange(10, dtype=np.int32)
            )
        )

        output_file = self.output_dir / "mseed3.mseed"

        with self.assertRaises(ValueError):
            mswrite(
                output_file,
                collection,
                msformat=3,
            )

    def test_mswrite_rejects_invalid_sid(self):
        """Test rejection of an invalid waveform source identifier."""
        record = self._make_shakelab_record(
            np.arange(10, dtype=np.int32)
        )
        record.head.sid = "INVALID"

        collection = StreamCollection()
        collection.append(record)

        output_file = self.output_dir / "invalid_sid.mseed"

        with self.assertRaises(ValueError):
            mswrite(
                output_file,
                collection,
                encoding=ENCODING_INT32,
            )

    def test_mswrite_msread_round_trip(self):
        """Test the complete high-level MiniSEED API round trip."""
        original = self._integer_signal()

        collection = StreamCollection()
        collection.append(
            self._make_shakelab_record(original)
        )

        output_file = self.output_dir / "high_level.mseed"

        mswrite(
            output_file,
            collection,
            encoding=ENCODING_STEIM2,
            reclen=256,
        )

        restored = msread(output_file)

        self.assertEqual(len(restored), 1)
        self.assertEqual(len(restored[0]), 1)

        output_record = restored[0][0]

        self.assertEqual(
            output_record.head.sid,
            "XX.TEST.00.BHZ",
        )
        self.assertAlmostEqual(
            output_record.head.delta,
            0.01,
        )
        self.assertAlmostEqual(
            output_record.head.time
            - Date("2026-200T12:00:00.0000"),
            0.0,
            places=4,
        )

        np.testing.assert_array_equal(
            output_record.data,
            original,
        )

    def test_mswrite_steim1_segmented_round_trip(self):
        """Test high-level STEIM1 writing across multiple records."""
        original = self._integer_signal()

        collection = StreamCollection()
        collection.append(
            self._make_shakelab_record(original)
        )

        output_file = self.output_dir / "segmented_steim1.mseed"

        mswrite(
            output_file,
            collection,
            encoding=ENCODING_STEIM1,
            reclen=256,
        )

        raw_records = msrawread(output_file)

        self.assertGreater(len(raw_records), 1)

        restored = np.concatenate([
            record.data
            for record in raw_records
        ])

        np.testing.assert_array_equal(restored, original)

    def test_partial_steim_limits_prepared_input(self):
        """Test that partial STEIM encoding limits input preparation."""
        original = np.arange(
            100000,
            dtype=np.int32,
        )

        prepared_sizes = []

        original_prepare = _prepare_steim_data

        def tracked_prepare(data):
            prepared_sizes.append(len(data))
            return original_prepare(data)

        with patch(
            "shakelab.signals.libio.mseed._prepare_steim_data",
            side_effect=tracked_prepare,
        ):
            _, _, sample_count, _ = _encode_steim_payload(
                original,
                encoding=ENCODING_STEIM2,
                record_length=256,
                data_offset=64,
                require_all=False,
            )

        max_frame_count = (256 - 64) // 64
        data_word_count = 13 + 15 * (max_frame_count - 1)
        expected_limit = data_word_count * 7

        self.assertEqual(
            prepared_sizes,
            [expected_limit],
        )
        self.assertLessEqual(
            sample_count,
            expected_limit,
        )
        self.assertLess(
            expected_limit,
            len(original),
        )


if __name__ == "__main__":
    unittest.main()