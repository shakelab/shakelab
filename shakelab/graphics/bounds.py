# ****************************************************************************
#
# Copyright (C) 2019-2023, ShakeLab Developers.
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
(External) utilities to create axis labels
"""

import math
import numpy as np


def logb(value, base=10):

    #return np.log(value) / np.log(base)
    return math.log(value, base)

def to_scientific_notation(value, base=10):
    """
    Decompose a number into a mantissa and an exponent of given base
    as for scientific notation.
    """
    exponent = math.floor(logb(value, base))
    mantissa = value / base ** exponent

    return mantissa, exponent

def exponent_range(axis_start, axis_end, base=10):
    """
    """
    m0, e0 = to_scientific_notation(axis_start, base)
    m1, e1 = to_scientific_notation(axis_end, base)

    if m0 > 1:
        e0 += 1

    return list(range(e0, e1 + 1))

def logb_ticks(axis_start, axis_end, base=10):
    """
    """
    exp_list = exponent_range(axis_start, axis_end, base)

    return [base**exp for exp in exp_list]

def lin_ticks(axis_start, axis_end, tick_num, inside_only=True):

    if tick_num < 2:
        tick_num = 2

    tick_min, tick_max, step = nice_bounds(axis_start, axis_end, tick_num)
    ticks = np.arange(tick_min, tick_max + step, step)

    if inside_only:
        ticks = ticks[ticks >= axis_start]
        ticks = ticks[ticks <= axis_end]

    return ticks

def nice_number(value, round_=False):
    """
    nice_number(value, round_=False) -> float
    
    Original code from "The Shmoo" of Stack overflow:
    https://stackoverflow.com/users/2226708/the-shmoo
    """
    mantissa, exponent = to_scientific_notation(value, 10)

    if round_:
        if mantissa < 1.5:
            fraction = 1.
        elif mantissa < 3.:
            fraction = 2.
        elif mantissa < 7.:
            fraction = 5.
        else:
            fraction = 10.
    else:
        if mantissa <= 1:
            fraction = 1.
        elif mantissa <= 2:
            fraction = 2.
        elif mantissa <= 5:
            fraction = 5.
        else:
            fraction = 10.

    return fraction * 10 ** exponent

def nice_bounds(axis_start, axis_end, num_ticks=10):
    """
    nice_bounds(axis_start, axis_end, num_ticks=10) -> tuple
    @return: tuple as (nice_axis_start, nice_axis_end, nice_tick_width)
    """
    axis_width = axis_end - axis_start

    if axis_width == 0:
        nice_tick = 0

    else:
        nice_range = nice_number(axis_width)
        nice_tick = nice_number(nice_range / (num_ticks - 1), round_=True)

        axis_start = math.ceil(axis_start / nice_tick) * nice_tick
        axis_end = math.floor(axis_end / nice_tick) * nice_tick

    return axis_start, axis_end, nice_tick

def auto_ticks(data_range, max_tick=10, tf_inside=False):
    """
    Tool function that automatically calculate optimal ticks based
    on range and the max number of ticks

    :param data_range:   range of data, e.g. [-0.1, 0.5]
    :param max_tick:     max number of ticks, an integer, default to 10
    :param tf_inside:    True to force ticks only inside the data range
    :return:             list of ticks

    Original code from "Shaobo Guan" of Stack overflow:
    https://stackoverflow.com/users/7470085/shaobo-guan
    """
    data_span = data_range[1] - data_range[0]

    # Scale of data as the order of 10, e.g. 1, 10, 100, 0.1, 0.01, ...
    scale = 10.0**np.floor(np.log10(data_span))

    # Possible tick sizes for normalized data in range [1, 10]
    list_tick_size_nmlz = [5.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01] 

    # Initial tick size for normalized data
    tick_size_nmlz = 1.0

    for i in range(len(list_tick_size_nmlz)):
        num_tick = data_span / scale / list_tick_size_nmlz[i]

        if num_tick > max_tick:
            tick_size_nmlz = list_tick_size_nmlz[i-1]
            break

    # Tick size for the original data
    tick_size = tick_size_nmlz * scale

    # List of ticks
    dr0 = data_range[0] / tick_size
    dr1 = data_range[1] / tick_size
    ticks = np.unique(np.arange(dr0, dr1).round()) * tick_size

    # Ticks allowed only within the given range
    if tf_inside:
        ticks = ticks[(ticks>=data_range[0]) * (ticks<=data_range[1])]

    return ticks