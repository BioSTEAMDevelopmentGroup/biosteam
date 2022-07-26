# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2022, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam import Stream
import yaml
import os

__all__ = ('preferences',)

class BioSTEAMDisplayPreferences:
    """Class containing preferences for BioSTEAM diagram and results display."""
    __slots__ = ('label_streams', 'autodisplay', 'minimal_nodes', 'number_path',
                 'profile', 'raise_exception', 'background_color', 'stream_color',
                 'label_color', 'label_color', 'depth_colors', 'stream_width',
                 'unit_color', 'unit_label_color', 'unit_periphery_color',
                 'fill_cluster')
    
    def __init__(self):
        #: Whether to label the ID of streams with sources and sinks in process 
        #: flow diagrams.
        self.label_streams = True
        
        #: Whether to automatically generate diagrams when displaying an object in the
        #: IPython console.
        self.autodisplay = True
        
        #: Whether to ignore unit graphics and display unit nodes as dots in process
        #: flow diagrams.
        self.minimal_nodes = False
        
        #: Whether to number unit operations according to their order in the system path.
        self.number_path = False
        
        #: Whether to clock the simulation time of unit operations in diagrams.
        self.profile = False
        
        #: Whether to raise exception regarding problems displaying graphviz diagrams.
        self.raise_exception = False
        
        #: Background color in graphviz diagrams.
        self.background_color = 'transparent'
        
        #: Color of streams in graphviz diagrams.
        self.stream_color = '#90918e'
        
        #: Color of stream labels in graphviz diagrams.
        self.label_color = '#90918e'
        
        #: Color of subsystem clusters in BioSTEAM graphviz diagrams.
        self.depth_colors = ['#7ac0836f']
        
        #: Property to scale stream widths in BioSTEAM graphviz diagrams.
        self.stream_width = 'F_mass'
        
        #: Unit node fill color in BioSTEAM graphviz diagrams.
        self.unit_color = '#555f69'
        
        #: Unit node label color in BioSTEAM graphviz diagrams.
        self.unit_label_color = 'white'
        
        #: Unit node periphery color in BioSTEAM graphviz diagrams.
        self.unit_periphery_color = '#90918e'
        
        #: Whether to fill subsystem boxes in BioSTEAM 'cluster' diagrams.
        self.fill_cluster = False
        
    def reset(self, save=False):
        """Reset to BioSTEAM defaults."""
        self.__init__()
        self.flow = 'kmol/hr'
        self.T = 'K'
        self.P = 'Pa'
        self.composition = False
        self.N = 7
        if save: self.save()
        
    @property
    def flow(self):
        """[str] Flow rate units."""
        return Stream.display_units.flow
    @flow.setter
    def flow(self, units):
        Stream.display_units.flow = units
        
    @property
    def T(self):
        """[str] Temperature units."""
        return Stream.display_units.T
    @T.setter
    def T(self, units):
        Stream.display_units.T = units

    @property
    def P(self):
        """[str] Pressure units."""
        return Stream.display_units.P
    @P.setter
    def P(self, units):
        Stream.display_units.P = units

    @property
    def composition(self):
        """[bool] Whether to show composition.""" 
        return Stream.display_units.composition
    @composition.setter
    def composition(self, composition):
        Stream.display_units.composition = composition

    @property
    def N(self):
        """[int] Number of compounds to display."""
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
                     cluster=('#7ac0836f',),
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
        
        self.label_streams = data['label_streams']
        self.autodisplay = data['autodisplay']
        self.minimal_nodes = data['minimal_nodes']
        self.number_path = data['number_path']
        self.profile = data['profile']
        self.raise_exception = data['raise_exception']
        self.background_color = data['background_color']
        self.stream_color = data['stream_color']
        self.label_color = data['label_color']
        self.depth_colors = data['depth_colors']
        self.flow = data['flow']
        self.T = data['T']
        self.P = data['P']
        self.composition = data['composition']
        self.N = data['N']
        self.stream_width = data.get('stream_width', 'F_mass')
        self.unit_color = data.get('unit_color', '#555f69')
        self.unit_label_color = data.get('unit_label_color', 'white')
        self.unit_periphery_color = data.get('unit_periphery_color', '#90918e')
        self.fill_cluster = data.get('fill_cluster', False)
        
    def save(self):
        folder = os.path.dirname(__file__)
        file = os.path.join(folder, 'preferences.yaml')
        with open(file, 'w') as file:
            dct = {i: getattr(self, i) for i in BioSTEAMDisplayPreferences.__slots__}
            dct['flow'] = self.flow
            dct['T'] = self.T
            dct['P'] = self.P
            dct['composition'] = self.composition
            dct['N'] = self.N
            yaml.dump(dct, file)

#: [BioSTEAMDisplayPreferences] All preferences for diagram and results display.
preferences = BioSTEAMDisplayPreferences()
# from warnings import filterwarnings; filterwarnings('ignore')
if not os.environ.get("DISABLE_PREFERENCES") == "1":
    try: preferences.autoload()
    except: pass 
