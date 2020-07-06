# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# LipidLynxX is Dual-licensed
#   For academic and non-commercial use: GPLv2 License:
#   For commercial use: please contact the SysMedOs team by email.
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: lipid annotations converter for large scale lipidomics and epilipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from dataclasses import dataclass
from typing import List


@dataclass
class Lipid:
    code: str
    formula: str
    elements: dict
    exact_mass: float
    smiles: str
    level: int
    modified: bool
    modifications: dict
    adducts: dict
    spectra: dict
    lipid_class: str


@dataclass
class FA(Lipid):
    lipid_class: str = "FA"


@dataclass
class PL(Lipid):
    lipid_sub_class: str
    fa_list: List[FA]
    fa1: FA
    fa2: FA
    position: bool


@dataclass
class GL(Lipid):
    lipid_sub_class: str


@dataclass
class TG(GL):
    fa_list: List[FA]
    fa1: FA
    fa2: FA
    fa3: FA


@dataclass
class DG(GL):
    fa_list: List[FA]
    fa1: FA
    fa2: FA


@dataclass
class MG(GL):
    fa_list: List[FA]
    fa1: FA
