# ****************************************************************************
#
# Copyright (C) 2019-2025, ShakeLab Developers.
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
Module to handle completeness information of seismic catalogues
"""

import csv

class CompletenessTable:
    """
    A class to represent a completeness table for a seismic catalog.
    """

    def __init__(self):
        """
        Initializes the CompletenessTable with an empty list.
        """
        self.table = []

    def add_entry(self, magnitude, year_start, year_end, completeness):
        """
        Adds an entry to the completeness table.

        Parameters
        ----------
        magnitude : float
            The magnitude of the earthquake.
        year_start : int
            The start year of the completeness period.
        year_end : int
            The end year of the completeness period.
        completeness : float
            The completeness percentage.
        """
        entry = {
            'magnitude': magnitude,
            'year_start': year_start,
            'year_end': year_end,
            'completeness': completeness
        }
        self.table.append(entry)
        self.sort_by_magnitude()

    def remove_entry(self, index):
        """
        Removes an entry from the completeness table by index.

        Parameters
        ----------
        index : int
            The index of the entry to remove.
        """
        if 0 <= index < len(self.table):
            del self.table[index]

    def update_entry(self, index, magnitude,
                     year_start, year_end, completeness):
        """
        Updates an existing entry in the completeness table.

        Parameters
        ----------
        index : int
            The index of the entry to update.
        magnitude : float
            The magnitude of the earthquake.
        year_start : int
            The start year of the completeness period.
        year_end : int
            The end year of the completeness period.
        completeness : float
            The completeness percentage.
        """
        if 0 <= index < len(self.table):
            self.table[index] = {
                'magnitude': magnitude,
                'year_start': year_start,
                'year_end': year_end,
                'completeness': completeness
            }
            self.sort_by_magnitude()

    def find_entry_by_magnitude(self, magnitude):
        """
        Finds entries in the completeness table by magnitude.

        Parameters
        ----------
        magnitude : float
            The magnitude to search for.

        Returns
        -------
        list
            A list of entries with the specified magnitude.
        """
        return [e for e in self.table if entry['magnitude'] == magnitude]

    def find_entry_by_year(self, year):
        """
        Finds entries in the completeness table by year.

        Parameters
        ----------
        year : int
            The year to search for.

        Returns
        -------
        list
            A list of entries that include the specified year.
        """
        return [e for e in self.table if entry['year_start'] <= year <= entry['year_end']]

    def sort_by_magnitude(self):
        """
        Sorts the completeness table by magnitude.
        """
        self.table.sort(key=lambda x: x['magnitude'])

    def sort_by_year(self):
        """
        Sorts the completeness table by year start.
        """
        self.table.sort(key=lambda x: x['year_start'])

    def calculate_average_completeness(self):
        """
        Calculates the average completeness of the table.

        Returns
        -------
        float
            The average completeness.
        """
        if not self.table:
            return 0.0
        total_completeness = sum(entry['completeness'] for entry in self.table)
        return total_completeness / len(self.table)

    def export_to_csv(self, filename):
        """
        Exports the completeness table to a CSV file.

        Parameters
        ----------
        filename : str
            The name of the CSV file.
        """
        with open(filename, 'w', newline='') as csvfile:
            fieldnames = ['magnitude', 'year_start',
                          'year_end', 'completeness']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for entry in self.table:
                writer.writerow(entry)

    def import_from_csv(self, filename):
        """
        Imports data from a CSV file into the completeness table.

        Parameters
        ----------
        filename : str
            The name of the CSV file.
        """
        with open(filename, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            self.table = [row for row in reader]
            self.sort_by_magnitude()

    def clear_table(self):
        """
        Clears all entries from the completeness table.
        """
        self.table = []

    def display(self):
        """
        Displays the completeness table in a formatted manner.
        """
        print(f"{'Magnitude':<10} {'Year Start':<10} {'Year End':<10} "
              f"{'Completeness':<15}")
        for entry in self.table:
            print(f"{entry['magnitude']:<10} {entry['year_start']:<10} "
                  f"{entry['year_end']:<10} {entry['completeness']:<15}")