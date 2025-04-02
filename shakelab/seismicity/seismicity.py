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
"""
import numpy as np


def poisson_probability(lambda_rate, investigation_time=50.):
    """
    Probability of “at-least” one occurrence of a Poisson process
    with a known constant mean rate (λ) and time interval (t).
    NOTE: clearly, rate is the inverse of the return period....
    """

    return 1 - np.exp(-lambda_rate * investigation_time) 

def poisson_rate(probability, investigation_time=50.):
    """
    Occurrence rate of “at-least” one event of a Poisson process with 
    known probability (P) time interval (t).
    NOTE: clearly, rate is the inverse of the return period....
    """

    return - np.log(1 - probability)/investigation_time

def sample_magnitude_distribution(aval, bval, mmin, mmax, duration=1.):
    """
    """
    # Computing cumulative annual rate for events larger than Mmin
    mmin_rate = (10.**aval)*(10.**(-bval*mmin)-10.**(-bval*mmax))

    # Magnitude distribution:
    # Derived as the inverse of the 
    # double-truncated cumulative GR
    # See Inverse transform sampling (also known as inversion sampling,
    # the inverse probability integral transform)

    Um = np.random.rand(int(mmin_rate * duration))
    C = 1.-(10.**(-bval*(mmax-mmin)))
    magnitude_samples = mmin-np.log10(1.-(Um*C))/bval

    return magnitude_samples

def sample_time_distribution(duration):
    """
    """
    # Seconds per year (Not YET counting leap years)
    secyear = duration*365.*24.*60.*60.

    timesec = np.sort(np.random.rand(len(magnitudes)) * secyear)

    return timesec

def generate_synthetic_catalogue(aval, bval, mmin, mmax, duration=1.):
    """
    Must define an initial date!

    also, it is probably worth to generalize the MFD as input
    instead of aval, bval.
    """

    magnitudes = sample_magnitude_distribution(aval, bval,
                                               mmin, mmax, duration)

    # Average interval time obtained
    # from the inverse of the cumulative
    # of the exponential distribution
    # (Poisson assumption)
    # TOCHECK!

    #Ut = np.random.rand(catlen)
    #dT = -np.log(1-Ut)/catlen
    #T = np.cumsum(dT)
    #secyear = duration*365.*24.*60.*60.
    #timesec = T*secyear



    return timesec, magnitudes
