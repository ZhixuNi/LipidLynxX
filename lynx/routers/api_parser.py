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

from fastapi import APIRouter, status

from lynx.controllers.parser import parse_lipid
from lynx.models.api_models import (
    InputStrData,
    LevelsData,
)

router = APIRouter()

default_levels = LevelsData(levels=["B1", "D1"])


# Get APIs
@router.get("/lipid/")
async def parse_lipid(lipid_name: str = "PLPC"):
    """
    Parse one lipid name from path parameter
    """
    return parse_lipid(lipid_name)


# Post APIs
@router.post("/str/", status_code=status.HTTP_201_CREATED)
async def parse_str(data: InputStrData):
    """
    Parse one lipid name from data
    """

    return parse_lipid(data.data)
