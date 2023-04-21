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
WX Widget to draw 2D data
"""

import wx
import numpy as np
from copy import deepcopy
from numpy.random import randn

from shakelab.gui.bounds import lin_ticks, logb_ticks


DEFAULT_PARAMS = {
    'top_margin' : 25,
    'bottom_margin' : 60,
    'left_margin' : 60,
    'right_margin' : 25,
    'background_colour' : 'white',
    'logx' : False,
    'logy' : False,
    'xlim' : 'auto',
    'ylim' : 'auto',
    'xaxis' : 'bottom',
    'yaxis' : 'left',
    'axis_line_width' : 1,
    'axis_line_colour' : 'black',
    'axis_bg_colour' : None,
    'xticks' : 'auto',
    'yticks' : 'auto',
    'xtick_labels' : 'auto',
    'ytick_labels' : 'auto',
    'xtick_spacing' : 60,
    'ytick_spacing' : 30,
    'xtick_length' : 5,
    'ytick_length' : 5,
    'font_size' : 14,
    'xgrid' : True,
    'ygrid' : True,
    'grid_colour' : 'light grey',
    'grid_line_style' : 'dot',
    'grid_line_width' : 2,
    'box' : True,
    'xlabel' : None,
    'ylabel' : None,
    'label_size' : 16,
    'ticks_outside' : True}

DATA_STRUCT = {
    'id' : None,
    'x' : None,
    'y' : None,
    'z' : None,
    'colour' : 'black',
    'width' : 1,
    'style' : 'solid',
    'label' : None}

WXSTYLE = {
    'solid' : wx.SOLID,
    'dot' : wx.DOT,
    '*' : wx.DOT,
    'long_dash' : wx.LONG_DASH,
    '--' : wx.LONG_DASH,
    'short_dash' : wx.SHORT_DASH,
    '-' : wx.SHORT_DASH,
    'dot_dash' : wx.DOT_DASH,
    '*-' : wx.DOT_DASH}


class DataStore():
    """
    Create an internal archive to store data
    and related plot settings.
    Note: data are plotted following the storage order.
    """
    def __init__(self):
        self.db = []

    def __iter__(self):
        self._i = 0
        return self

    def __next__(self):
        try:
            item = self.db[self._i]
            self._i += 1
        except IndexError:
            raise StopIteration
        return item

    def __getitem__(self, idx):
        if idx < len(self.db):
            return self.db[idx]
        else:
            raise ValueError('No data set found')

    def __len__ (self):
        return len(self.db)

    def add(self, x, y, **kwargs):

        item = deepcopy(DATA_STRUCT)

        id = -1
        while True:
            id += 1
            if id not in self.ids:
                item['id'] = id
                break

        item['x'] = np.array(x)
        item['y'] = np.array(y)

        for k, v in kwargs.items():
            if k in DATA_STRUCT.keys():
                item[k] = v

        self.db.append(item)

    def remove(self, id):
        pass

    @property
    def ids(self):
        return [i['id'] for i in self.db]

    @property
    def xmin(self):
        if not self.db:
            return 0
        else:
            return min([min(i['x']) for i in self.db])

    @property
    def xmax(self):
        if not self.db:
            return 1
        else:
            return max([max(i['x']) for i in self.db])

    @property
    def ymin(self):
        if not self.db:
            return 0
        else:
            return min([min(i['y']) for i in self.db])

    @property
    def ymax(self):
        if not self.db:
            return 1
        else:
            return max([max(i['y']) for i in self.db])

    def xlim(self):
        return self.xmin, self.xmax

    def ylim(self):
        return self.ymin, self.ymax


class BasePlot(wx.Panel):

    def __init__(self, parent, tile_grid=(1,1), params={}, **kwargs):

        # Platform specific settings
        if wx.Platform == '__WXGTK__':
            pass
        elif wx.Platform == '__WXMAC__':
            pass
        elif wx.Platform == '__WXMSW__':
            pass

        # Initialise grid
        rows, cols = tile_grid
        self.graph = [[Graph()]*cols]*rows

        self.CTRL = False

        wx.Panel.__init__(self, parent, size=(-1, -1))
        self.parent = parent

        self.Bind(wx.EVT_PAINT, self.OnPaint)
        # self.Bind(wx.EVT_SIZE, self.OnSize)

        self.Bind(wx.EVT_MOUSE_EVENTS, self.OnMouseEvent)
        self.Bind(wx.EVT_MOTION, self.OnMouseEvent)

        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseLeftClick)

        self.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)
        self.Bind(wx.EVT_KEY_UP, self.OnKeyUp)

    def Plot(self, x, y, tile=(0,0), params={}, **kwargs):
        """
        """
        self.data[tile[0]][tile[1]].add(x, y, *args, **kwargs)

    def SetProperty(self, tile=(0,0), params={}, **kwargs):
        """
        """
        params = {**params, **kwargs}

        for key, value in params.items():
            if key in DEFAULT_PARAM.keys():
                self.params[tile[0]][tile[1]][key] = value

    def ExportImage(self, file_name):
        """
        """
        if self._bitmap is not None:
            self._bitmap.SaveFile(file_name, wx.BITMAP_TYPE_PNG)

    def DrawGraph(self):

        size = self.GetClientSize()

        panel_width = max(size[0], 1)
        panel_height = max(size[1], 1)

        i_margin = self.params['left_margin'] + self.params['right_margin']
        j_margin = self.params['top_margin'] + self.params['bottom_margin']

        # Plot area width and height
        axis_width = max(panel_width - i_margin, 0)
        axis_height = max(panel_height - j_margin, 0)

        i0 = self.params['left_margin'] + 1
        j0 = self.params['top_margin'] + axis_height

        # Background canvas
        self._bitmap = wx.Bitmap(panel_width, panel_height, 32)

        with wx.BufferedPaintDC(self, self._bitmap) as dc:

            dc.SetBackground(wx.Brush(self.params['background_colour']))
            dc.Clear()

            #gcdc = wx.GCDC(dc)
            #gc = gcdc.GetGraphicsContext()
            #gc = wx.GraphicsContext.Create(dc)

            #dc.SetBackgroundMode(wx.TRANSPARENT)

            graph.Draw()

    def OnPaint(self, event):

        self.DrawGraph()
        self.Refresh()

    def OnSize(self, event):

        self.DrawGraph()
        self.Refresh()

    def OnMouseEvent(self, event):

        if event.LeftDown():
            if self.CTRL:
                print(pos.x, pos.y)
            else:
                print('Not active')
        else:
            pass

    def OnMouseLeftClick(self, event):

        if self.CTRL:
            pos = event.GetPosition()
            print(pos.x, pos.y)
        else:
            print('Not active')

    def OnKeyDown(self, event):

        #key = event.GetKeyCode()
        #if chr(key).upper() == 'A':
        #    self.CTRL = True

        if event.ControlDown():
            self.CTRL = True

        event.Skip()

    def OnKeyUp(self, event):

        self.CTRL = False
        event.Skip()


# ---------------------------------------------------------------------------

class Graph():
    """
    """
    def __init__(self, origin, params={}, **kwargs):
        """
        """
        self.data = DataStore()

        self.i0 = origin[0]       # axis origin's i index
        self.j0 = origin[1]       # axis origin's j index

        # Override default properties
        self.params = deepcopy(DEFAULT_PARAMS)
        self.SetProperty(params, **kwargs)

        self.xmin = 0
        self.xmax = 1

        self.ymin = 0
        self.ymax = 1

        self.xtick = np.array([])
        self.ytick = np.array([])

        self.lx = 0
        self.ly = 0

    def Draw(self, dc):
        """
        """
        self.SetDataBounds()

        self.SetAxisTicks()

        self.DrawBackground(dc)

        self.DrawGrid(dc)

        self.DrawData(dc)

        self.DrawAxis(dc)

        self.DrawLabels(dc)

        self.DrawBox(dc)

    def SetProperty(self, params={}, **kwargs):
        """
        Override default settings
        with user-defined properties
        """
        params = {**params, **kwargs}

        for key, value in params.items():
            if key in self.params.keys():
                self.params[key] = value

    def GetProperty(self, key):
        """
        """
        if key in self.params.keys():
            return self.params[key]
        else:
            print('Property not found')

    def AddData(self, x, y, **kwargs):
        """
        """
        self.data.add(x, y, *args, **kwargs)

    def SetCanvasSize(self, size):
        """
        """
        self.aw = max(size[0], 1)  # axis width
        self.ah = max(size[1], 1)  # axis height

        self.i1 = self.i0 + self.aw - 1
        self.j1 = self.j0 - self.ah + 1

    def SetDataBounds(self):
        """
        """
        if self.params['xlim'] == 'auto':
            xmin, xmax = self.data.xlim()
        else:
            xmin, xmax = self.params['xlim']

        if self.params['ylim'] == 'auto':
            ymin, ymax = self.data.ylim()
        else:
            ymin, ymax = self.params['ylim']

        self.xmin = xmin
        self.xmax = xmax

        self.ymin = ymin
        self.ymax = ymax

    def SetAxisTicks(self):
        """
        """
        xnum = round(self.aw / self.params['xtick_spacing'])
        ynum = round(self.ah / self.params['ytick_spacing'])

        if xnum > 0:
            if self.params['logx']:
                self.xtick = logb_ticks(self.xmin, self.xmax)
            else:
                self.xtick = lin_ticks(self.xmin, self.xmax, xnum)

        if ynum > 0:
            if self.params['logy']:
                self.ytick = logb_ticks(self.ymin, self.ymax)
            else:
                self.ytick = lin_ticks(self.ymin, self.ymax, ynum)

    def DrawBackground(self, dc):
        """
        """
        bg_colour = self.params['axis_bg_colour']

        if bg_colour is not None:

            dc.SetPen(wx.Pen(bg_colour, 1))
            dc.SetBrush(wx.Brush(bg_colour))
            dc.DrawRectangle(self.i0, self.j1, self.aw, self.ah)

    def DrawGrid(self, dc):
        """
        """
        self.dc.SetPen(wx.Pen(self.params['grid_colour'],
                              self.params['grid_line_width'],
                              WXSTYLE[self.params['grid_line_style']]))

        if self.params['xgrid']:
            for xtick in self.xtick:

                xi = self.XToPix(xtick)

                if self.i0 < xi < self.i1:
                    dc.DrawLine(xi, self.j0, xi, self.j1)

        if self.params['ygrid']:
            for ytick in self.ytick:

                yi = self.YToPix(ytick)

                if self.j1 < yi < self.j0:
                    dc.DrawLine(self.i0, yi, self.i1, yi)

    def DrawData(self, dc):
        """
        """
        if (self.aw * self.ah) > 0:

            polygon = wx.Rect(self.i0, self.j1, self.aw, self.ah)
            dc.SetClippingRegion(polygon)

            for d in self.data:

                style = WXSTYLE.get(d['style'], 'solid')
                dc.SetPen(wx.Pen(d['colour'], d['width'], style))

                xi = self.XToPix(d['x'])
                yi = self.YToPix(d['y'])

                dc.DrawLines([xy for xy in zip(xi, yi)])

            dc.DestroyClippingRegion()

    def DrawAxis(self, dc):
        """
        """
        dc.SetPen(wx.Pen(self.params['axis_line_colour'],
                         self.params['axis_line_width']))

        font = dc.GetFont()
        font.SetPointSize(self.params['font_size'])
        dc.SetFont(font)

        tdir = -1 if self.params['ticks_outside'] else 1

        for xtick in self.xtick:

            xi = self.XToPix(xtick)
            yi = self.YToPix(self.ymin)
            dyi = tdir * self.params['xtick_length']

            dc.DrawLine(xi, yi, xi, yi - dyi)

            xtick_label = '{}'.format(round(xtick, 15))
            label_width, label_height = dc.GetTextExtent(xtick_label)

            dxi = rint(label_width / 2)
            dyi = 5 - (dyi if tdir < 0 else 0)

            self.ly = max(self.ly, label_height + dyi)

            dc.DrawText(xtick_label, xi - dxi, yi + dyi)

        for ytick in self.ytick:

            xi = self.XToPix(self.xmin)
            yi = self.YToPix(ytick)
            dxi = tdir * self.params['ytick_length']

            dc.DrawLine(xi, yi, xi + dxi, yi)

            ytick_label = '{}'.format(round(ytick, 15))
            label_width, label_height = dc.GetTextExtent(ytick_label)

            dxi = rint(label_width) + 5 - (dxi if tdir < 0 else 0)
            dyi = rint(label_height / 2)

            self.lx = max(self.lx, dxi)

            dc.DrawText(ytick_label, xi - dxi, yi - dyi)

    def DrawLabels(self, dc):
        """
        """
        #font = dc.GetFont()
        #font.SetPointSize(self.params['label_size'])

        font = wx.Font(pointSize=self.params['label_size'],
                       family=wx.FONTFAMILY_DEFAULT,
                       style=wx.FONTSTYLE_NORMAL,
                       weight = wx.FONTWEIGHT_NORMAL,
                       faceName = '',
                       encoding = wx.FONTENCODING_DEFAULT)

        dc.SetFont(font)

        xlabel = self.params['xlabel']
        ylabel = self.params['ylabel']

        if xlabel is not None:

            lw, lh = dc.GetTextExtent(xlabel)
            xl = rint(self.i0 + self.aw/2 - lw/2)
            yl = rint(self.j0 + self.ly)

            dc.DrawRotatedText(xlabel, xl, yl, 0)

        if ylabel is not None:

            lw, lh = dc.GetTextExtent(ylabel)
            xl = rint(self.i0 - self.lx - lh)
            yl = rint(self.j0 - self.ah/2 + lw/2)

            dc.DrawRotatedText(ylabel, xl, yl, 90)

    def DrawBox(self, dc):
        """
        """
        dc.SetPen(wx.Pen(self.params['axis_line_colour'],
                         self.params['axis_line_width']))

        if self.params['box']:
            dc.DrawLine(self.i0, self.j0, self.i1, self.j0)
            dc.DrawLine(self.i0, self.j1, self.i1, self.j1)
            dc.DrawLine(self.i0, self.j0, self.i0, self.j1)
            dc.DrawLine(self.i1, self.j0, self.i1, self.j1)

    def XToPix(self, x):
        """
        """
        if self.params['logx']:
            xn = np.log10(x / self.xmin) / np.log10(self.xmax / self.xmin)
        else:
            xn = (x - self.xmin) / (self.xmax - self.xmin)

        return rint(xn * (self.aw - 1)) + self.i0

    def YToPix(self, y):
        """
        """
        if self.params['logy']:
            yn = np.log10(self.ymax / y) / np.log10(self.ymax / self.ymin)
        else:
            yn = (self.ymax - y) / (self.ymax - self.ymin)

        return rint(yn * (self.ah - 1)) + self.j1

def rint(value):
    """
    Round and return integer
    """
    return np.rint(value).astype(int)

def unique(x, y):
    """
    x and y must be numpy arrays
    """
    i = [0] + [n for n in range(1, len(x)) if x[n] != x[n-1]]
    return x[i], y[i]
