# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam import Stream
from thermosteam.units_of_measure import parse_units_notation
import yaml
import os

__all__ = ('preferences', 'TemporaryPreferences')

class DisplayPreferences:
    """
    All preferences for BioSTEAM diagram and results display.

    Examples
    --------
    >>> from biosteam import preferences
    >>> preferences.show()
    DisplayPreferences:
    label_streams: True
    autodisplay: True
    minimal_nodes: False
    number_path: False
    profile: False
    raise_exception: False
    background_color: 'transparent'
    stream_color: '#90918e'
    label_color: '#90918e'
    depth_colors: ['#f98f609f']
    stream_width: 'F_mass'
    unit_color: '#555f69'
    unit_label_color: 'white'
    unit_periphery_color: '#90918e'
    fill_cluster: False
    graphviz_format: 'svg'
    tooltips_full_results: False
    graphviz_html_height: {'system': ('400px', '600px'), 'unit': ('225px', '400px')}
    flow: 'kmol/hr:.3g'
    T: 'K:.5g'
    P: 'Pa:.6g'
    composition: False
    N: 7
    sort: False
    
    """
    __slots__ = ('label_streams', 'autodisplay', 'minimal_nodes', 'number_path',
                 'profile', 'raise_exception', 'background_color', 'stream_color',
                 'label_color', 'label_color', 'depth_colors', 'stream_width',
                 'unit_color', 'unit_label_color', 'unit_periphery_color',
                 'fill_cluster', 'graphviz_format', 'tooltips_full_results',
                 'graphviz_html_height')
    
    def __init__(self):
        #: Whether to label the ID of streams with sources and sinks in process 
        #: flow diagrams.
        self.label_streams: bool = True
        
        #: Whether to automatically generate diagrams when displaying an object in the
        #: IPython console.
        self.autodisplay: bool = True
        
        #: Whether to ignore unit graphics and display unit nodes as dots in process
        #: flow diagrams.
        self.minimal_nodes: bool = False
        
        #: Whether to number unit operations according to their order in the system path.
        self.number_path: bool = False
        
        #: Whether to clock the simulation time of unit operations in diagrams.
        self.profile: bool = False
        
        #: Whether to raise exception regarding problems displaying graphviz diagrams.
        self.raise_exception: bool = False
        
        #: Background color in graphviz diagrams.
        self.background_color: str = 'transparent'
        
        #: Color of streams in graphviz diagrams.
        self.stream_color: str = '#90918e'
        
        #: Color of stream labels in graphviz diagrams.
        self.label_color: str = '#90918e'
        
        #: Color of subsystem clusters in BioSTEAM graphviz diagrams.
        self.depth_colors: list[str] = ['#f98f609f']
        
        #: Property to scale stream widths in BioSTEAM graphviz diagrams.
        self.stream_width: str = 'F_mass'
        
        #: Unit node fill color in BioSTEAM graphviz diagrams.
        self.unit_color: str = '#555f69'
        
        #: Unit node label color in BioSTEAM graphviz diagrams.
        self.unit_label_color: str = 'white'
        
        #: Unit node periphery color in BioSTEAM graphviz diagrams.
        self.unit_periphery_color: str = '#90918e'
        
        #: Whether to fill subsystem boxes in BioSTEAM 'cluster' diagrams.
        self.fill_cluster: bool = False
        
        #: Image format of BioSTEAM graphviz diagrams.
        self.graphviz_format: str = 'svg'
        
        #: Whether to add full results in tooltips by inserting java script into graphviz html outputs.
        self.tooltips_full_results: bool = False
        
        #: Displayed height of graphviz html diagrams without and with full results.
        self.graphviz_html_height: dict[str, tuple[str, str]] = {
            'system': ('400px', '600px'),
            'unit': ('225px', '400px'),
        }
        
    def temporary(self):
        """Return a TemporaryPreferences object that will revert back to original
        preferences after context management."""
        return TemporaryPreferences()
        
    def reset(self, save=False):
        """Reset to BioSTEAM defaults."""
        self.__init__()
        self.flow = 'kmol/hr'
        self.T = 'K'
        self.P = 'Pa'
        self.composition = False
        self.N = 7
        self.sort = False
        if save: self.save()
        
    @property
    def flow(self) -> str:
        """Flow rate units and notation."""
        return ":".join([Stream.display_units.flow, Stream.display_notation.flow])
    @flow.setter
    def flow(self, units_notation):
        units, notation = parse_units_notation(units_notation)
        if units is not None:
            Stream.display_units.flow = units
        if notation is not None:
            Stream.display_notation.flow = notation
        
    @property
    def T(self) -> str:
        """Temperature units and notation."""
        return ":".join([Stream.display_units.T, Stream.display_notation.T])
    @T.setter
    def T(self, units_notation):
        units, notation = parse_units_notation(units_notation)
        if units is not None:
            Stream.display_units.T = units
        if notation is not None:
            Stream.display_notation.T = notation

    @property
    def P(self) -> str:
        """Pressure units and notation."""
        return ":".join([Stream.display_units.P, Stream.display_notation.P])
    @P.setter
    def P(self, units_notation):
        units, notation = parse_units_notation(units_notation)
        if units is not None:
            Stream.display_units.P = units
        if notation is not None:
            Stream.display_notation.P = notation

    @property
    def composition(self) -> bool:
        """Whether to show composition.""" 
        return Stream.display_units.composition
    @composition.setter
    def composition(self, composition):
        Stream.display_units.composition = composition

    @property
    def sort(self) -> bool:
        """Whether to sort flows in decreasing order.""" 
        try: return Stream.display_units.sort
        except: return False
    @sort.setter
    def sort(self, sort):
        try: Stream.display_units.sort = sort
        except: pass

    @property
    def N(self) -> int:
        """Number of compounds to display."""
        return Stream.display_units.N
    @N.setter
    def N(self, N):
        Stream.display_units.N = N
        
    def update(self, *, save=False, **kwargs):
        for i, j in kwargs.items(): setattr(self, i, j)
        if save: self.save()
        
    def _set_mode(self, stream, label, bg, cluster, unit_color, 
                  unit_label_color, unit_periphery_color, fill_cluster, save):
        self.background_color = bg
        self.stream_color = stream
        self.label_color = label
        self.depth_colors = cluster
        self.unit_color = unit_color
        self.unit_label_color = unit_label_color
        self.unit_periphery_color = unit_periphery_color
        self.fill_cluster = fill_cluster
        if save: self.save()
    
    def classic_mode(self, 
                     stream='#90918e', 
                     label='#90918e', 
                     bg='transparent',
                     cluster=('#f98f609f',),
                     unit_color='#555f69',
                     unit_label_color='white',
                     unit_periphery_color='none',
                     fill_cluster=False,
                     save=False):
        """Set diagram display colors to classic mode."""
        self._set_mode(stream, label, bg, cluster, unit_color, 
                       unit_label_color, unit_periphery_color,
                       fill_cluster, save)
    
    def dark_mode(self, stream='#98a2ad', label='#e5e5e5', bg='transparent',
                  cluster=['#5172512f'], unit_color='#555f69', 
                  unit_label_color='white', unit_periphery_color='none',
                  fill_cluster=True, save=False):
        """Set diagram display colors to dark mode."""
        self._set_mode(stream, label, bg, cluster, unit_color, unit_label_color,
                       unit_periphery_color, fill_cluster, save)
    
    def light_mode(self, stream='#4e4e4e', label='#4e4e4e', bg='#ffffffff',
                   cluster=['#7ac0832f'], unit_color='white:#CDCDCD', 
                   unit_label_color='black', unit_periphery_color='#4e4e4e',
                   fill_cluster=True, save=False):
        """Set diagram display colors to light mode."""
        self._set_mode(stream, label, bg, cluster, unit_color, unit_label_color, 
                       unit_periphery_color, fill_cluster, save)
    
    night_mode = dark_mode
    day_mode = light_mode
    
    def autoload(self):
        folder = os.path.dirname(__file__)
        file = os.path.join(folder, 'preferences.yaml')
        with open(file, 'r') as stream: 
            data = yaml.full_load(stream)
            assert isinstance(data, dict), 'yaml file must return a dict' 
        self.update(**data)
        
    def to_dict(self):
        """Return dictionary of all preferences."""
        dct = {i: getattr(self, i) for i in preferences.__slots__}
        dct['flow'] = self.flow
        dct['T'] = self.T
        dct['P'] = self.P
        dct['composition'] = self.composition
        dct['N'] = self.N
        dct['sort'] = self.sort
        
        return dct
        
    def save(self):
        """Save preferences."""
        folder = os.path.dirname(__file__)
        file = os.path.join(folder, 'preferences.yaml')
        with open(file, 'w') as file:
            dct = self.to_dict()
            yaml.dump(dct, file)

    def show(self):
        """Print all specifications."""
        dct = self.to_dict()
        print(f'{type(self).__name__}:\n' + '\n'.join([f"{i}: {repr(j)}" for i, j in dct.items()])) 
    _ipython_display_ = show


class TemporaryPreferences:
    
    def __enter__(self):
        dct = self.__dict__
        dct.update({i: getattr(preferences, i) for i in preferences.__slots__})
        dct['flow'] = preferences.flow
        dct['T'] = preferences.T
        dct['P'] = preferences.P
        dct['composition'] = preferences.composition
        dct['N'] = preferences.N
        dct['sort'] = preferences.sort
        return preferences
        
    def __exit__(self, type, exception, traceback):
        preferences.update(**self.__dict__)
        if exception: raise exception

#: 
preferences: DisplayPreferences = DisplayPreferences()

if os.environ.get("FILTER_WARNINGS"):
    from warnings import filterwarnings; filterwarnings('ignore')
if not os.environ.get("DISABLE_PREFERENCES") == "1":
    try: preferences.autoload()
    except: pass 
