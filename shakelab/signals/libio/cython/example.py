import libmseed

# Import binary buffer from an example file
#file_path = 'libmseed/libmseed-3.1.3/test/data/reference-testdata-steim2.mseed3'
file_path = 'test_data/ANMO.IU.mseed'

with open(file_path, 'rb') as f:
    byte_stream = f.read()

# Initialise a Miniseed object
ms = libmseed.MiniSeed()

# Read the buffer
ms.read(byte_stream, verbose=1)

# Import data from a dictionary
import numpy as np
nsamp = 100000
data = 1000*np.random.randn(nsamp) * np.exp(-np.arange(nsamp)/(nsamp/10))

record_dict = {
    'network' : 'OX',
    'station' : 'ABCD',
    'location' : '00',
    'channel' : 'HHZ',
    'starttime' : '2010-01-10T08:23:45.019538Z',
    'rate' : 200,
    'nsamp' : nsamp,
    'data' : data
}

ms.import_record(record_dict)

# Export data to a list of dictionaries.
rec_list = ms.export_records()

# Save binary buffer to a file
if 1:
    byte_stream = ms.write()

    # Write bytes to file
    with open("my_file.mseed", "wb") as binary_file:
        binary_file.write(byte_stream)