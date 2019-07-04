# -*- coding: utf-8 -*-
"""
Created on Tue May 14 14:35:55 2019

@author: yoelr
"""
import re
import biosteam as bst
__all__ = ('elementname',)

def elementname(element):
    if element:
        if isinstance(element, type):
            return re.sub(r"\B([A-Z])", r" \1", element.__name__.replace('_', ' ')).capitalize()
        elif isinstance(element, str):
            return element.replace('_', ' ')
        else:
            if isinstance(element, bst.Stream):
                return 'Stream-' + element.ID.replace('_', ' ')
            else:
                return element.line + '-' + element.ID.replace('_', ' ')
    else:
        return 'None'