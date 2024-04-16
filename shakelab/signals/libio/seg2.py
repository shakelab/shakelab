import struct

class SEG2Parser:
    def __init__(self, filename):
        self.filename = filename
        self.header_info = {}
        self.trace_data = []

    def parse(self):
        with open(self.filename, 'rb') as seg2_file:
            # Read and parse the text header (ASCII header).
            text_header = seg2_file.read(3200).decode('utf-8').strip()
            self.header_info['TextHeader'] = text_header

            # Read and parse the binary header (400 bytes).
            binary_header = seg2_file.read(400)
            self.header_info['BinaryHeader'] = binary_header

            # Define format strings for extracting header values from the binary header.
            header_format = '<10s10s20s10s10s10s10s20s4s12s'
            header_size = struct.calcsize(header_format)
            header_data = struct.unpack(header_format, binary_header)

            # Store specific header values in the dictionary.
            self.header_info['FieldRecord'] = header_data[0].strip()
            self.header_info['TraceNumber'] = int(header_data[1].strip())
            self.header_info['ReceiverPoint'] = int(header_data[2].strip())
            self.header_info['ReceiverPointIndex'] = int(header_data[3].strip())
            self.header_info['SourcePoint'] = int(header_data[4].strip())
            self.header_info['SourcePointIndex'] = int(header_data[5].strip())
            self.header_info['DataUse'] = header_data[6].strip()
            self.header_info['Distance'] = float(header_data[7].strip())
            self.header_info['Elevation'] = float(header_data[8].strip())
            self.header_info['CoordinateUnits'] = header_data[9].strip()

            # Read and parse trace data.
            while True:
                trace_header = seg2_file.read(240)
                if not trace_header:
                    break

                # Define format string for extracting trace header values.
                trace_header_format = '<ihihhhhhhhhhh'
                trace_header_data = struct.unpack(trace_header_format, trace_header)

                # Store trace header information in a dictionary.
                trace_info = {
                    'TraceSequenceLine': trace_header_data[0],
                    'TraceNumber': trace_header_data[1],
                    'DataUse': trace_header_data[2],
                    'Offset': trace_header_data[3],
                    'ReceiverElevation': trace_header_data[4],
                    'SourceElevation': trace_header_data[5],
                    'SourceDepth': trace_header_data[6],
                    'DatumElevationReceiverGroup': trace_header_data[7],
                    'DatumElevationSource': trace_header_data[8],
                    'WaterDepthSource': trace_header_data[9],
                    'WaterDepthGroup': trace_header_data[10],
                    'ScalerToBeAppliedToAllCoordinates': trace_header_data[11],
                    'SeismicData': []  # List to store seismic data samples
                }

                # Determine the data format based on header information.
                data_use = trace_info['DataUse']
                num_samples = trace_info['TraceNumber']

                # Define format string based on the data use code.
                if data_use == '1':
                    sample_format = '<' + 'h' * num_samples  # 16-bit signed integers
                elif data_use == '3':
                    sample_format = '<' + 'f' * num_samples  # 32-bit floating-point
                elif data_use == '4':
                    sample_format = '<' + 'd' * num_samples  # 64-bit floating-point
                else:
                    # Handle additional data formats as needed.
                    raise NotImplementedError(f"Data use code {data_use} not supported.")

                # Read and parse the seismic data samples.
                seismic_data = struct.unpack(sample_format, seg2_file.read(struct.calcsize(sample_format)))

                trace_info['SeismicData'] = seismic_data
                self.trace_data.append(trace_info)
    
    def get_header_info(self):
        return self.header_info

    def get_trace_data(self):
        return self.trace_data

# Example usage:
if __name__ == "__main__":
    seg2_filename = "your_seg2_file.seg2"
    seg2_parser = SEG2Parser(seg2_filename)
    seg2_parser.parse()
    
    header_info = seg2_parser.get_header_info()
    trace_data = seg2_parser.get_trace_data()
    
    # Print the parsed header information and trace information as needed.
    # ...
