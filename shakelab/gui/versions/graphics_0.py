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
WX Widget to draw various data
"""

import wx
import numpy as np
from copy import deepcopy
from numpy.random import randn

from bounds import nice_ticks, auto_ticks

PARAM = {'xlim' : 'auto',
         'ylim' : 'auto',
         'background_colour' : 'white',
         'top_margin' : 25,
         'bottom_margin' : 60,
         'left_margin' : 60,
         'right_margin' : 25,
         'xaxis' : 'bottom',
         'yaxis' : 'left',
         'axis_line_width' : 1,
         'axis_line_colour' : 'black',
         'axis_bg_colour' : None,
         'xtick_spacing' : 60,
         'ytick_spacing' : 30,
         'xtick_length' : 5,
         'ytick_length' : 5,
         'font_size' : 14,
         'xgrid' : True,
         'ygrid' : True,
         'grid_colour' : 'light grey',
         'grid_line_width' : 1,
         'box' : True,
         'xlabel' : None,
         'ylabel' : None,
         'label_size' : 16,
         'xlabel' : 'Time (s)',
         'ylabel' : 'Amplitude'}

DATA_STRUCT = {'id' : None,
               'x' : None,
               'y' : None,
               'z' : None,
               'colour' : 'black',
               'width' : 1,
               'style' : 'solid'}

WXSTYLE = {'solid' : wx.SOLID,
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
    and plotting properties.
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

    def add(self, x, y, *args, **kwargs):

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


class DataPlot(wx.Panel):

    def __init__(self, parent, params={}, **kwargs):

        self.data = DataStore()
        self.CTRL = False
        self.buffer = None

        self.SetProperty(params, **kwargs)

        wx.Panel.__init__(self, parent, size=(-1, -1))
        self.parent = parent
        self.SetBackgroundColour(PARAM['background_colour'])

        self.Bind(wx.EVT_PAINT, self.OnPaint)
        # self.Bind(wx.EVT_SIZE, self.OnSize)

        self.Bind(wx.EVT_MOUSE_EVENTS, self.OnMouseEvent)
        self.Bind(wx.EVT_MOTION, self.OnMouseEvent)

        self.Bind(wx.EVT_LEFT_DOWN, self.OnMouseLeftClick)

        self.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown)
        self.Bind(wx.EVT_KEY_UP, self.OnKeyUp)

    def Plot(self, x, y, *args, **kwargs):

        self.data.add(x, y, *args, **kwargs)

    def SetProperty(self, params={}, **kwargs):

        params = {**params, **kwargs}

        for k, v in params.items():
            if k in PARAM.keys():
                PARAM[k] = v

    def GetProperty(self, key):

        if key in selg.cfg.keys():
            return PARAM(key)
        else:
            print('Property not found')

    def ExportImage(self, file_name):

        if self.buffer is not None:
            self.buffer.SaveFile(file_name, wx.BITMAP_TYPE_PNG)

    def _SetSize(self):

        (w, h) = self.GetClientSize()

        # Canvas width and height
        self._CW = max(w, 1)
        self._CH = max(h, 1)

        # Margins to the plot area
        self._TM = PARAM['top_margin']
        self._BM = PARAM['bottom_margin']
        self._LM = PARAM['left_margin']
        self._RM = PARAM['right_margin']

        # Plot area width and height
        self._AW = max(self._CW - self._LM - self._RM, 0)
        self._AH = max(self._CH - self._TM - self._BM, 0)

        self._I0 = self._LM + 1
        self._I1 = self._LM + self._AW
        self._J0 = self._TM + self._AH
        self._J1 = self._TM + 1

        self._XL = 0
        self._YL = 0

    def _SetDataBounds(self):

        if PARAM['xlim'] == 'auto':
            xmin, xmax = self.data.xlim()
        else:
            xmin, xmax = PARAM['xlim']

        if PARAM['ylim'] == 'auto':
            ymin, ymax = self.data.ylim()
        else:
            ymin, ymax = PARAM['ylim']

        self._XMIN = xmin
        self._XMAX = xmax

        self._YMIN = ymin
        self._YMAX = ymax

    def _ToPixel(self, x, y):

        xn = (x - self._XMIN) / (self._XMAX - self._XMIN)
        yn = (self._YMAX - y) / (self._YMAX - self._YMIN)

        xi = np.rint(xn * (self._AW - 1)).astype(int) + self._LM + 1
        yi = np.rint(yn * (self._AH - 1)).astype(int) + self._TM + 1

        return xi, yi

    def _GenTicks(self):

        self._xtick = np.array([])
        self._ytick = np.array([])

        xnum = round(self._AW / PARAM['xtick_spacing'])
        ynum = round(self._AH / PARAM['ytick_spacing'])

        if xnum > 1:
            self._xtick = nice_ticks(self._XMIN, self._XMAX, xnum, True)

        if ynum > 1:
            self._ytick = nice_ticks(self._YMIN, self._YMAX, ynum, True)

    def _Uniques(self, x, y):
        """
        x and y must be numpy arrays
        """
        i = [0] + [n for n in range(1, len(x)) if x[n] != x[n-1]]
        return x[i], y[i]

    def OnPaint(self, event):

        self._SetSize()
        self._SetDataBounds()

        # Background canvas
        self.buffer = wx.Bitmap(self._CW, self._CH, 32)

        with wx.BufferedPaintDC(self, self.buffer) as dc:

            #dc = wx.GCDC(dc)
            #dc.SetBackgroundMode(wx.TRANSPARENT)

            self._DrawBackground(dc)

            self._GenTicks()
            self._DrawGrid(dc)

            if self._AW > 0 and self._AH > 0:
                self._DrawData(dc)

            self._DrawTicks(dc)
            self._DrawBox(dc)
            self._DrawLabels(dc)

            self.Refresh()

    def _DrawBackground(self, dc):
        """
        """
        dc.SetBackground(wx.Brush(PARAM['background_colour']))
        dc.Clear()

        if PARAM['axis_bg_colour'] is not None:

            dc.SetPen(wx.Pen(PARAM['axis_bg_colour'], 1))
            dc.SetBrush(wx.Brush(PARAM['axis_bg_colour']))

            dc.DrawRectangle(self._I0, self._J1,
                             self._AW, self._AH)

    def _DrawGrid(self, dc):
        """
        """
        dc.SetPen(wx.Pen(PARAM['grid_colour'],
                         PARAM['grid_line_width'],
                         WXSTYLE['*']))

        for xtick in self._xtick:
            xi, yi = self._ToPixel(xtick, self._YMIN)

            if PARAM['xgrid'] and xi > self._I0 and xi < self._I1:
                dc.DrawLine(xi, self._J0, xi, self._J1)

        for ytick in self._ytick:
            xi, yi = self._ToPixel(self._XMIN, ytick)

            if PARAM['ygrid'] and yi > self._J1 and yi < self._J0:
                dc.DrawLine(self._I0, yi, self._I1, yi)

    def _DrawData(self, dc):
        """
        """
        dc.SetClippingRegion(wx.Rect(self._I0, self._J1,
                                     self._AW, self._AH))

        for d in self.data:

            xi, yi = self._ToPixel(d['x'], d['y'])

            style = WXSTYLE.get(d['style'], 'solid')
            dc.SetPen(wx.Pen(d['colour'], d['width'], style))

            dc.DrawLines([xy for xy in zip(xi, yi)])

        dc.DestroyClippingRegion()

    def _DrawTicks(self, dc):
        """
        """
        dc.SetPen(wx.Pen(PARAM['axis_line_colour'],
                         PARAM['axis_line_width']))

        font = dc.GetFont()
        font.SetPointSize(PARAM['font_size'])
        dc.SetFont(font)

        for xtick in self._xtick:

            xi, yi = self._ToPixel(xtick, self._YMIN)

            dc.DrawLine(xi, yi, xi, yi - PARAM['xtick_length'])

            #xtick_label = '{:.2E}'.format(xtick)
            xtick_label = '{}'.format(round(xtick, 15))

            label_width, label_height = dc.GetTextExtent(xtick_label)

            dxi = rint(label_width / 2)
            dyi = 5

            self._YL = max(self._YL, label_height + dyi)

            dc.DrawText(xtick_label, xi - dxi, yi + dyi)

        for ytick in self._ytick:

            xi, yi = self._ToPixel(self._XMIN, ytick)

            dc.DrawLine(xi, yi, xi + PARAM['ytick_length'], yi)

            #ytick_label = '{:.2E}'.format(ytick)
            ytick_label = '{}'.format(round(ytick, 15))

            label_width, label_height = dc.GetTextExtent(ytick_label)

            dxi = rint(label_width) + 5
            dyi = rint(label_height / 2)

            self._XL = max(self._XL, dxi)

            dc.DrawText(ytick_label, xi - dxi, yi - dyi)

    def _DrawLabels(self, dc):
        """
        """
        #font = dc.GetFont()
        #font.SetPointSize(PARAM['label_size'])

        font = wx.Font(pointSize=16,
                       family=wx.FONTFAMILY_DEFAULT,
                       style=wx.FONTSTYLE_NORMAL,
                       weight = wx.FONTWEIGHT_NORMAL,
                       faceName = '',
                       encoding = wx.FONTENCODING_DEFAULT)

        dc.SetFont(font)

        xlabel = PARAM['xlabel']
        ylabel = PARAM['ylabel']

        if xlabel is not None:

            lw, lh = dc.GetTextExtent(xlabel)
            xl = rint(self._I0 + self._AW/2 - lw/2)
            yl = rint(self._J0 + self._YL)

            dc.DrawRotatedText(xlabel, xl, yl, 0)

        if ylabel is not None:

            lw, lh = dc.GetTextExtent(ylabel)
            xl = rint(self._I0 - self._XL - lh)
            yl = rint(self._J0 - self._AH/2 + lw/2)

            dc.DrawRotatedText(ylabel, xl, yl, 90)

    def _DrawBox(self, dc):

        dc.SetPen(wx.Pen(PARAM['axis_line_colour'],
                         PARAM['axis_line_width']))

        if PARAM['box']:
            dc.DrawLine(self._I0, self._J0, self._I1, self._J0)
            dc.DrawLine(self._I0, self._J1, self._I1, self._J1)
            dc.DrawLine(self._I0, self._J0, self._I0, self._J1)
            dc.DrawLine(self._I1, self._J0, self._I1, self._J1)

    def OnSize(self, event):

        self.OnPaint(event)
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

def rint(value):
    """
    """
    return np.rint(value).astype(int)


def generate_data(data_length):

        if 1:
            x = np.linspace(0, 1, data_length)
            y = np.random.randn(1, data_length)[0]

        if 0:
            x = np.array([1,2,3,4])
            y = np.array([0, 10, -30, 0])

        if 0:
            freq = 10
            x = np.linspace(0, 1, data_length)
            y = np.sin(2*np.pi*x*freq)
        return x, y


class MainWindow(wx.Frame):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs, size=(800, 400))

        panel = wx.Panel(self)

        sizer = wx.BoxSizer(wx.VERTICAL)

        tp = TracePlot(panel)
        sizer.Add(tp, 1, wx.EXPAND|wx.ALL, 10)
        panel.SetSizer(sizer)

        # Example data
        x, y = generate_data(400)
        tp.Plot(x, y, colour='blue', width=2)

        x, y = generate_data(20)
        x += 0
        y += 0
        tp.Plot(x, y, colour='black', width=2, style='solid')

        tp.SetProperty(axis_bg_colour='#F0F0F0')
        #tp.SetProperty(top_margin=5)
        #tp.SetProperty(bottom_margin=30)
        #tp.SetProperty(left_margin=15)
        #tp.SetProperty(right_margin=5)

        tp.SetProperty(xlim='auto')
        #tp.SetProperty(ylim=(-5, 5))

        self.SetTitle("Trace Plot Example")
        self.Centre()
        self.Show()

    def OnQuit(self, event):

        self.Close(True)


def main():

    app = wx.App(False)
    ex = MainWindow(None)
    app.MainLoop()


if __name__ == '__main__':
    main()
