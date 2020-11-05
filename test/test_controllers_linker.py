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
from datetime import datetime
from typing import List

import pandas as pd

from lynx.controllers.linker import link_lipid, link_lipids
from lynx.utils.file_handler import get_abs_path


# lipid_name_lst = ["PLPC"]
lipid_name_lst = ["PC(16:0/18:2)"]
lipid_names_lst = [["palmitic acid", "SM d18:1/24:0", "SM(d18:1/24:0)"]]
file_info_lst = [[r"test/test_input/SearchLipids.csv", "SEARCH_NAME"]]


@pytest.mark.parametrize("lipid_name", lipid_name_lst)
def test_link_str(lipid_name: str):

    loop = asyncio.get_event_loop()
    results_dct = loop.run_until_complete(link_lipid(lipid_name, direct_search=True))
    print(results_dct)


@pytest.mark.parametrize("lipid_names", lipid_names_lst)
def test_link_list(lipid_names: List[str]):

    loop = asyncio.get_event_loop()
    results_df = loop.run_until_complete(link_lipids(lipid_names))
    print(results_df.T)


@pytest.mark.parametrize("file_info", file_info_lst)
def test_link_file(file_info: List[str]):

    abs_file_path_str = get_abs_path(file_info[0])
    search_col = file_info[1]

    output_name = (
        f"LipidLynxX-Linker_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.csv"
    )

    df_dct = pd.read_csv(abs_file_path_str).to_dict(orient="list")
    lipid_names = df_dct.get(search_col)

    loop = asyncio.get_event_loop()
    results_df = loop.run_until_complete(link_lipids(lipid_names, direct_search=True))
    print(results_df.T)
    results_df.to_csv(output_name)
