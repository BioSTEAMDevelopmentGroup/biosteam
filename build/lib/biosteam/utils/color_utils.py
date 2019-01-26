# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 17:19:33 2018

This module includes classes and functions relating string coloring and format.

@author: Guest Group
"""

from colorama import Fore, Style
from colorpalette import Color, Palette

__all__ = ('CS', 'color_scheme')

# %% Classes for coloring

ansicolor = Color.ansicolor
CS = color_scheme = Palette(
        exception=ansicolor(Fore.RED + Style.BRIGHT, 'brightred'),
        info=ansicolor(Fore.BLUE + Style.BRIGHT, 'brightblue'),
        note=ansicolor(Fore.CYAN + Style.BRIGHT, 'brightcyan'))

CS.dark_warning = Color('darkgoldenrod')
CS.warning = Color('gold')
CS.request = Color('peachpuff')
CS.dim = ansicolor(Fore.WHITE + Style.NORMAL, ID='dim')
CS.reset = ansicolor(Style.RESET_ALL, ID='reset')
CS.ul = Color(style='underline')
CS.num = Color('rednum', '#f74f4f')
CS.str = Color('greenstr', '#89e189')

