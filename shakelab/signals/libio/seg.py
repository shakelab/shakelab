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
\import numpy as np
from struct import pack, unpack

class SEG1:
    def __init__(self, file_path=None, byte_order='>'):
        """
        Classe base per la lettura e la scrittura di file SEG (SEG-1).

        :param file_path: Percorso del file SEG (SEG-1)
        :param byte_order: Ordine dei byte ('>' per big-endian, '<' per little-endian)
        """
        self.header = {}
        self.data = []
        self.byte_order = byte_order

        if file_path:
            self.read(file_path)

    def read(self, file_path):
        """
        Legge un file SEG (SEG-1) dal disco.

        :param file_path: Percorso del file SEG (SEG-1)
        """
        with open(file_path, 'rb') as f:
            self.header = self._read_header(f)
            self.data = self._read_data(f)

    def _read_header(self, file):
        """
        Legge l'header del file SEG (SEG-1).

        :param file: File object
        :return: Dizionario contenente le informazioni dell'header
        """
        header = {}
        # Campi dell'header (esempio basato su una struttura tipica SEG-1)
        header['job_id'] = unpack(self.byte_order + 'I', file.read(4))[0]
        header['line_number'] = unpack(self.byte_order + 'I', file.read(4))[0]
        header['reel_number'] = unpack(self.byte_order + 'I', file.read(4))[0]
        header['traces_per_ensemble'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['aux_traces_per_ensemble'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['sample_interval'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['sample_interval_orig'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['samples_per_trace'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['samples_per_trace_orig'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['data_sample_format'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['ensemble_fold'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['trace_sorting_code'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['vertical_sum_code'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['sweep_frequency_start'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['sweep_frequency_end'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['sweep_length'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['sweep_type'] = unpack(self.byte_order + 'H', file.read(2))[0]
        header['trace_number'] = unpack(self.byte_order + 'I', file.read(4))[0]
        # Aggiungi altri campi specifici dell'header SEG-1...

        return header

    def _read_data(self, file):
        """
        Legge i dati sismici dal file SEG (SEG-1).

        :param file: File object
        :return: Lista contenente i dati sismici
        """
        data = []
        while True:
            trace_header = file.read(240)  # Legge l'header della traccia (240 byte)
            if not trace_header:
                break
            # Esempio: leggere dati di traccia (modificare in base alla struttura del file SEG-1)
            trace_data = np.frombuffer(file.read(1024), dtype=self.byte_order + 'f')
            data.append(trace_data)
        return data

    def write(self, file_path):
        """
        Scrive un file SEG (SEG-1) sul disco.

        :param file_path: Percorso del file di output SEG (SEG-1)
        """
        with open(file_path, 'wb') as f:
            # Scrittura dell'header
            self._write_header(f)
            # Scrittura dei dati
            for trace in self.data:
                f.write(trace.tobytes())

    def _write_header(self, file):
        """
        Scrive l'header del file SEG (SEG-1).

        :param file: File object
        """
        file.write(pack(self.byte_order + 'I', self.header['job_id']))
        file.write(pack(self.byte_order + 'I', self.header['line_number']))
        file.write(pack(self.byte_order + 'I', self.header['reel_number']))
        file.write(pack(self.byte_order + 'H', self.header['traces_per_ensemble']))
        file.write(pack(self.byte_order + 'H', self.header['aux_traces_per_ensemble']))
        file.write(pack(self.byte_order + 'H', self.header['sample_interval']))
        file.write(pack(self.byte_order + 'H', self.header['sample_interval_orig']))
        file.write(pack(self.byte_order + 'H', self.header['samples_per_trace']))
        file.write(pack(self.byte_order + 'H', self.header['samples_per_trace_orig']))
        file.write(pack(self.byte_order + 'H', self.header['data_sample_format']))
        file.write(pack(self.byte_order + 'H', self.header['ensemble_fold']))
        file.write(pack(self.byte_order + 'H', self.header['trace_sorting_code']))
        file.write(pack(self.byte_order + 'H', self.header['vertical_sum_code']))
        file.write(pack(self.byte_order + 'H', self.header['sweep_frequency_start']))
        file.write(pack(self.byte_order + 'H', self.header['sweep_frequency_end']))
        file.write(pack(self.byte_order + 'H', self.header['sweep_length']))
        file.write(pack(self.byte_order + 'H', self.header['sweep_type']))
        file.write(pack(self.byte_order + 'I', self.header['trace_number']))
        # Aggiungi altri campi specifici dell'header SEG-1...

# Esempio di utilizzo
# seg_file = SEG1('percorso/del/file.seg')
# seg_file.read('percorso/del/file.seg')
# print(seg_file.header)
# print(seg_file.data)
