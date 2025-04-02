# ****************************************************************************
#
# Copyright (C) 2019-2024, ShakeLab Developers.
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
"""
import pickle
import logging

# Set up logging
logging.basicConfig(level=logging.WARNING)


def write_to_pickle(stream_collection, filename: str) -> None:
    """
    Writes a StreamCollection object to a pickle file.
    
    :param stream_collection: The StreamCollection object to write.
    :param filename: The filename (with path if necessary) to write the
        object to.
    :raises IOError: If there is an issue writing the file.
    """
    try:
        with open(filename, 'wb') as f:
            pickle.dump(stream_collection, f)
        logging.info(f"StreamCollection written to {filename}")
    except IOError as e:
        logging.error(f"Failed to write to file {filename}: {e}")
        raise


def read_from_pickle(filename: str):
    """
    Reads a StreamCollection object from a pickle file.
    
    :param filename: The filename (with path if necessary) to read the 
        object from.
    :return: The StreamCollection object.
    :raises IOError: If there is an issue reading the file.
    :raises pickle.UnpicklingError: If there is an issue unpickling the data.
    """
    try:
        with open(filename, 'rb') as f:
            stream_collection = pickle.load(f)
        logging.info(f"StreamCollection read from {filename}")
        return stream_collection
    except (IOError, pickle.UnpicklingError) as e:
        logging.error(f"Failed to read from file {filename}: {e}")
        raise


