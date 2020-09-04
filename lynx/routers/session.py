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

from fastapi import APIRouter, HTTPException, Response

from lynx.models.defaults import (
    default_temp_folder,
    default_temp_max_days,
    default_temp_max_files,
)
from lynx.utils.log import app_logger
from lynx.utils.file_handler import clean_temp_folder
from lynx.utils.session_tools import get_session_output, is_session_finished

router = APIRouter()

removed_files = clean_temp_folder(
    default_temp_folder, default_temp_max_days, default_temp_max_files
)
if removed_files:
    app_logger.info(
        f"Remove temporary output files older than {default_temp_max_days} days..."
    )
    for removed_file in removed_files:
        app_logger.info(f"File removed: {removed_file}")


# Get APIs
@router.get("/status/{session_id}", status_code=201, include_in_schema=False)
async def is_finished(session_id: str, response: Response) -> bool:
    """
    check if a session is finished or not
    curl -X GET "http://127.0.0.1:1399/session/status/{session_id}" -H  "accept: application/json"
    """
    session_finished = is_session_finished(session_id)
    if session_finished:
        response.status_code = 200
        return True
    else:
        response.status_code = 404
        return False


@router.get("/output/{session_id}", status_code=201, include_in_schema=False)
async def get_output(session_id: str, response: Response) -> dict:
    """
    get the output of a session
    curl -X GET "http://127.0.0.1:1399/session/output/{session_id}" -H  "accept: application/json"
    """
    session_output = await get_session_output(session_id)
    if session_output:
        response.status_code = 200
    else:
        response.status_code = 404
    return session_output
