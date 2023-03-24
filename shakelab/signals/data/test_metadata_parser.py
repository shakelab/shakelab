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

import shakelab.signals.metadata as slmd

# example file n.1: contains info on station T01 of network TT 
# (all codes are dummies). It has many channels, each with different location 
# codes and working intervals (corresponding to different instrumental 
# responses)
file_in = 'TT_T01.xml'


# load all the information contained in T01.xml at level 'response'
#------------------------------------------------------------------
md_t01 = slmd.Metadata(mdfile=file_in, level='response')


# create new Metadata() objects using different selections:
#----------------------------------------------------------

# information for HHE channels
md_t01_HHE = md_t01.select('TT.T01..HHE')

# information for HHE channels, at a given time
time = '2019-10-08T18:30:00.000000Z'
md_t01_HHE_time = md_t01.select('TT.T01..HHE', time=time)

# information for HHE channels with locationcode 01
md_t01_HHE = md_t01.select('TT.T01.01.HHE')

# the same dictionary stored in md_t01_HHE_time.data can also be accessed using 
# select_nslc
md_t01_HHE_time.data == slmd.select_nslc(md_t01.data, 'TT.T01..HHE', time=time)

# to access the data at level 'channel' now we can use __getitem__
chan_info = md_t01_HHE_time['Channel']
print(chan_info)