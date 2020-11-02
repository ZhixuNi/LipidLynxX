# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
#
# LipidLynxX is using GPL V3 License
#
# Please cite our publication in an appropriate form.
#   LipidLynxX preprint on bioRxiv.org
#   Zhixu Ni, Maria Fedorova.
#   "LipidLynxX: a data transfer hub to support integration of large scale lipidomics datasets"
#   DOI: 10.1101/2020.04.09.033894
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import asyncio
import pytest
from typing import List

from lynx.controllers.linker import link_lipids


lipid_name_lst = ["HexCer(18:1;O2/24:0)"]
lipid_names_lst = [["palmitic acid", "SM d18:1/24:0", "SM(d18:1/24:0)"]]


@pytest.mark.parametrize("lipid_name", lipid_name_lst)
def test_link_str(lipid_name: str):

    loop = asyncio.get_event_loop()
    results_df = loop.run_until_complete(link_lipids([lipid_name]))
    print(results_df.T)


@pytest.mark.parametrize("lipid_names", lipid_names_lst)
def test_link_list(lipid_names: List[str]):

    loop = asyncio.get_event_loop()
    results_df = loop.run_until_complete(link_lipids(lipid_names))
    print(results_df.T)
