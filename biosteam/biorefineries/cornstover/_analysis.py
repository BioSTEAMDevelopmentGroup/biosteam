# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 22:30:47 2019

@author: yoelr
"""
from biosteam.biorefineries.cornstover.model import cornstover_model as model_cs

N_samples = 5000
rule = 'L'

samples = model_cs.sample(N_samples, rule)
model_cs.load_samples(samples)
model_cs.evaluate()
model_cs.table.to_excel('Monte Carlo cornstover.xlsx')
metric_names = ('Minimum ethanol selling price',)
spearman = model_cs.spearman(metric_names)
spearman.to_excel("Spearman correlation cornstover.xlsx")