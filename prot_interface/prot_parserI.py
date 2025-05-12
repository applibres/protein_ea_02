#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22-11-2024

@author: Rolando Armas
    Phage Therapy Group
    Yachay Tech
"""

from optparse import OptionParser


def parse_value(string):
    ''' convert string to int, float or string '''
    try:
        return int(string)
    except ValueError:
        try:
            return float(string)
        except ValueError:
            return string         

def parse_params(param_str):
    ''' parse a param string to a dict '''
    dict_params = {}
    if param_str:
        for key_value_str in param_str.split(','):
            key, value = key_value_str.split('=')
            if key in ['pmops']:
                dict_params[key] = np.fromstring(str(value), dtype=float, sep=':')
            else:
                dict_params[key] = parse_value(value)
    return dict_params
