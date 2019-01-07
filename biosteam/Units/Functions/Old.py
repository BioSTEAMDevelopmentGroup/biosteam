#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 22:45:26 2017

@author: Yoel Rene Cortes-Pena
"""


def constant_solute_solvent_ratio(self):
    """
         top_Rspecies = list of species that correspond to Rs
         top_Rs = list of molar ratios of solute to solvent in top
         bot_Rspecies = list of species that correspond to bot_Rs
         bot_Rs =  list of molar ratios of solute to solvent in bot
         top_solvents = list of species that correspond to top_split
         top_split = list of splits for top
         bot_solvents = list of species that correspond to bot_split
         bot_split = list of splits for bot
    """
    feed = self.ins[0]
    top, bot = self.outs
    top_Rspecies, top_Rs, bot_Rspecies, bot_Rs, top_solvents, top_split, bot_solvents, bot_split = (
        (self.kwargs[i] for i in ('top_Rspecies', 'top_Rs', 'bot_Rspecies', 'bot_Rs', 'top_solvents', 'top_split', 'bot_solvents', 'bot_split')))
    top_index = [feed.Settings.ID_index_dictionary[specie]
                for specie in top_solvents]
    bot_index = [feed.Settings.ID_index_dictionary[specie]
                for specie in bot_solvents]
    top.mol[top_index] = feed.mol[top_index]*top_split
    bot.mol[top_index] = feed.mol[top_index]-top.mol[top_index]
    bot.mol[bot_index] = feed.mol[bot_index]*bot_split
    top.mol[bot_index] = feed.mol[bot_index]-bot.mol[bot_index]
    topnetsolv = sum(top.mol[top_index])
    botnetsolv = sum(bot.mol[bot_index])
    top_Rindex = [feed.Settings.ID_index_dictionary[specie]
                 for specie in top_Rspecies]
    bot_Rindex = [feed.Settings.ID_index_dictionary[specie]
                 for specie in bot_Rspecies]
    l1 = len(top_Rindex)
    for i in range(l1):
        index = top_Rindex[i]
        R = top_Rs[index]
        top.mol[index] = R*topnetsolv
    l2 = len(bot_Rindex)
    for i in range(l2):
        index = bot_Rindex[i]
        R = bot_Rs[index]
        bot.mol[index] = R*botnetsolv
