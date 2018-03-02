# msutil.py: Utility MeasurementSet functions
# Copyright (C) 2011
# Associated Universities, Inc. Washington DC, USA.
# Copyright (C) 2016
# Joint Institute for VLBI ERIC, Dwingeloo, The Netherlands
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
from casac import casac

def addCorrectedData(msname):
    """ Add CORRECTED_DATA column to an MS.

    """
    tb = casac.table()

    # Open the MS
    tb.open(msname, nomodify=False)
    cnames = tb.colnames()
    # Get the description of the DATA column.
    try:
        cdesc = tb.getcoldesc('DATA')
    except:
        raise ValueError('Column DATA does not exist')
    hasTiled = False
    # Determine if the DATA storage specification is tiled.
    hasTiled = False
    try:
        dminfo = None
        for rec in tb.getdminfo("DATA"):
            if recs[rec]['COLUMNS'][0] == 'DATA':
                dminfo = recs[rec]
        if dminfo['TYPE'][:5] == 'Tiled':
            hasTiled = True
    except:
        hasTiled = False
    # Use TiledShapeStMan if needed.
    if not hasTiled:
        dminfo = {'TYPE': 'TiledShapeStMan', 'SPEC': {'DEFAULTTILESHAPE':[4,32,128]}}
    # Add the columns (if not existing). Use the description of the DATA column.
    if 'CORRECTED_DATA' in cnames:
        print("Column CORRECTED_DATA not added; it already exists")
    else:
        dminfo['NAME'] = 'correcteddata'
        cdesc['comment'] = 'The corrected data column'
        tb.addcols({'CORRECTED_DATA': cdesc}, dminfo)
    # Flush the table to make sure it is written.
    tb.flush()
