# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:19:33 2018

This module includes classes and functions relating string coloring and format.

@author: Guest Group
"""
from colorpalette import Color, Palette

__all__ = ('CS', 'color_scheme')

# %% Classes for coloring

ansicolor = Color.ansicolor
CS = color_scheme = Palette(
        exception=ansicolor('\x1b[31m\x1b[1m', 'brightred'),
        info=ansicolor('\x1b[34m\x1b[1m', 'brightblue'),
        note=ansicolor('\x1b[36m\x1b[1m', 'brightcyan'))

CS.dark_warning = Color('darkgoldenrod')
CS.warning = Color('gold')
CS.request = Color('peachpuff')
CS.dim = ansicolor('\x1b[37m\x1b[22m', ID='dim')
CS.reset = ansicolor('\x1b[0m', ID='reset')
CS.ul = Color(style='underline')
CS.num = Color('rednum', '#f74f4f')
CS.str = Color('greenstr', '#89e189')
CS.next = CS.info
CS.bold = ansicolor('\033[1m', ID='bold')
CS.italic = ansicolor('\x1b[3m', ID='italic')