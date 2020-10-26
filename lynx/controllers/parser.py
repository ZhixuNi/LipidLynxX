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

import re

from pydantic import parse_obj_as

from lynx.controllers.converter import Converter
from lynx.controllers.decoder import Decoder
from lynx.controllers.encoder import Encoder
from lynx.models.lipid import LipidType
from lynx.utils.log import app_logger


def parse_lipid(lipid_name: str):
    lynx_converter = Converter()
    converted_results = lynx_converter.convert_str(input_str=lipid_name)
    converted_name = converted_results.output
    decoder = Decoder()
    encoder = Encoder(
        logger=app_logger,
    )
    extracted_info = decoder.extract(converted_name)
    best_input_rule = encoder.get_best_rule(lipid_name)
    best_id = encoder.convert(lipid_name)
    parsed_info = {"id": best_id, "input_style": best_input_rule}
    if extracted_info:
        for p in extracted_info:
            p_info = extracted_info[p]
            for in_r in p_info:
                if re.search("#LipidLynxX", in_r, re.IGNORECASE):
                    parsed_info["id"] = converted_name
                    r_info = p_info[in_r]
                    segments = r_info.get("segments", {})
                    parsed_info["lipid_class"] = segments.get("CLASS", [""])[0]
                    parsed_info["residues"] = r_info.get("residues", {})

                    residues_separator_level = r_info.get("residues", {}).get(
                        "residues_separator_level"
                    )

                    if residues_separator_level == "S":
                        parsed_info["exact_sn_position"] = True
                    else:
                        parsed_info["exact_sn_position"] = False

    return parsed_info


def detect_style(lipid_name: str) -> str:

    encoder = Encoder(
        logger=app_logger,
    )
    input_style = encoder.get_best_rule(lipid_name)

    return input_style


if __name__ == "__main__":
    # ls = (
    #     "SM(18:1{4E};2OH{1R,3S}<2OH{8S,18R},oxo{10}>/18:2{8E,11E}<2OH{7S,18R},oxo{12}>)"
    # )
    ls = (
        "SM 16:0"
    )
    exp = parse_lipid(ls)
    print(exp)

    br = detect_style(ls)
    print(f"Auto detect rule: [{br}] for lipid: {ls}")
