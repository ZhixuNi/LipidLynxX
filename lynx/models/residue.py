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

import json
import os

import regex as re
from jsonschema import Draft7Validator, RefResolver

from lynx.models.defaults import default_output_rules, res_schema, res_schema_path
from lynx.models.modification import merge_mods, Modifications
from lynx.utils.log import app_logger
from lynx.utils.params_loader import load_output_rule
from lynx.utils.toolbox import check_json


class Residue(object):
    def __init__(
        self,
        residue_info: dict,
        schema: str = "lynx_residues",
        output_rules: dict = default_output_rules,
        nomenclature: str = "LipidLynxX",
        logger=app_logger,
    ):
        self.logger = logger
        self.export_rule = load_output_rule(output_rules, nomenclature)
        self.res_rule = self.export_rule.get("RESIDUE", None)
        self.res_rule_orders = self.res_rule.get("RESIDUE_INFO", {}).get("ORDER", [])
        self.res_separators = self.export_rule.get("SEPARATOR", [])
        self.res_info = residue_info
        self.res_level = residue_info.get("level", {})
        # transfer the m/d/t label in Cer/SM/SP from link section into O section
        self.__replace_mdt__()

        self.schema = schema
        self.type = "FattyAcid"
        resolver = RefResolver(
            referrer=res_schema, base_uri=f"file://{os.path.dirname(res_schema_path)}/"
        )
        self.validator = Draft7Validator(res_schema, resolver=resolver)

        self.mod_info = self.res_info.get("mod_info_sum", {}).get("info", {})
        self.mod_level = self.res_info.get("mod_info_sum", {}).get("level", 0)

        self.mod_obj = Modifications(
            self.res_info.get("mod_info_sum", {}),
            o_count=self.res_info.get("o_count", 0),
            nomenclature=nomenclature,
        )
        self.sum_mod_info = self.mod_obj.info

        if float(self.mod_level) > 0:
            self.is_modified = True
        else:
            self.is_modified = False

        if self.is_modified and self.sum_mod_info:
            self.linked_levels = self.sum_mod_info.get("linked_levels", ["0"])
        else:
            self.linked_levels = ["0"]

        self.linked_ids = self.__post_init__()

    def __replace_mdt__(self):
        link = self.res_info.get("link", "")
        if link:
            existed_o_count = self.res_info.get("o_count", 0)
            # convert m/d/t into corresponding O section
            if link.lower() == "m":
                self.res_info["link"] = ""
                self.res_info["o_count"] = 1 + existed_o_count
            elif link.lower() == "d":
                self.res_info["link"] = ""
                self.res_info["o_count"] = 2 + existed_o_count
            elif link.lower() == "t":
                self.res_info["link"] = ""
                self.res_info["o_count"] = 3 + existed_o_count
            else:
                pass
        else:
            pass

    def __post_init__(self):
        res_str_dct = {}
        num_o = self.res_info.get("o_count", 0)
        for lv in self.linked_levels:
            res_str = ""
            for o in self.res_rule_orders:
                if o in self.res_info or o in self.res_separators or o in ["SUM_MODS"]:
                    if o == "o_count":
                        if num_o > 0:
                            o_seg_rgx = self.res_rule.get("RESIDUE", {}).get("o_count")
                            # print("o_seg_rgx", o_seg_rgx)
                            if o_seg_rgx:
                                if num_o == 1:
                                    if re.match(o_seg_rgx, str(num_o)):
                                        res_str += str(num_o)
                                    elif re.match(o_seg_rgx, "1"):
                                        res_str += "1"
                                    elif re.match(o_seg_rgx, "O"):
                                        res_str += "O"
                                    else:
                                        res_str += "O"
                                else:
                                    if re.match(o_seg_rgx, str(num_o)):
                                        res_str += str(num_o)
                                    elif re.match(o_seg_rgx, f"{num_o}O"):
                                        res_str += f"{num_o}O"
                                    elif re.match(o_seg_rgx, f"O{num_o}"):
                                        res_str += f"O{num_o}"
                                    else:
                                        res_str += str(num_o)
                            else:
                                res_str += str(num_o)
                        else:
                            pass

                    elif o.upper().endswith("_SEPARATOR"):
                        res_str += self.res_separators.get(o, "")
                        if num_o == 0:
                            res_str = res_str.strip(
                                self.res_separators.get("O_SEPARATOR", "")
                            )
                    elif re.search("BRACKET", o.upper()):
                        res_str += self.res_separators.get(o, "")
                    elif o in ["MODS", "MOD"]:
                        res_str += self.sum_mod_info.get("linked_ids", {}).get(lv, "")
                    elif o in ["SUM_MODS"]:
                        sum_o_seg = self.sum_mod_info.get("linked_ids", {}).get(lv, "")
                        if num_o > 0 and re.match(r"\d?O\d?", sum_o_seg):
                            pass
                        else:
                            res_str += sum_o_seg
                    else:
                        res_str += str(self.res_info.get(o, ""))
            na_brackets_lst = [r"\<\>", r"\{\}", r"\[\]", r"\(\)"]
            for b in na_brackets_lst:
                res_str = re.sub(b, "", res_str)
            res_str_dct[lv] = res_str.strip(";")

        return res_str_dct

    def to_json(self):
        fa_lite_info_dct = self.fa_info_dct
        fa_lite_info_dct.pop("mod_obj", None)
        fa_json_str = json.dumps(fa_lite_info_dct)
        if check_json(self.validator, json.loads(fa_json_str, logger=self.logger)):
            return fa_json_str
        else:
            raise Exception(f"JSON Schema check FAILED. Schema {self.schema}")


def merge_residues(
    residues_order: list,
    residues_info: dict,
    schema: str = "lynx_residues",
    output_rules: dict = default_output_rules,
    nomenclature: str = "LipidLynxX",
) -> Residue:

    sum_res_dct = {}
    if isinstance(residues_info, dict):
        pass
    else:
        raise TypeError(
            f"Requires multiple Residues in dict, "
            f"got type: {type(residues_info)} for {residues_info}"
        )

    all_mod_lst = [
        residues_info[residue_name].get("info", {}).get("mod_info_sum", {})
        for residue_name in residues_info
    ]
    sum_mods_obj = merge_mods(all_mod_lst)

    for res in residues_order:
        res_info = residues_info.get(res, {})
        link = res_info.get("LINK")
        if res.upper().startswith("P-") and link == "P-":
            res_info["LINK"] = "O-"
            res_info["DB_COUNT"] = res_info.get("DB_COUNT") + 1
        for res_seg in res_info:
            if re.search(r"MOD", res_seg):
                pass
            else:
                if res_seg not in sum_res_dct:
                    sum_res_dct[res_seg] = res_info[res_seg]
                else:
                    existed_count = sum_res_dct.get(res_seg, None)
                    res_seg_count = res_info.get(res_seg, None)
                    if res_seg_count:
                        if isinstance(existed_count, int) and isinstance(
                            res_seg_count, int
                        ):
                            sum_res_dct[res_seg] = res_seg_count + existed_count
                        elif isinstance(existed_count, str) and isinstance(
                            res_seg_count, str
                        ):
                            sum_res_dct[res_seg] = res_seg_count + existed_count
                        else:
                            raise TypeError
                    else:
                        pass

    sum_res_dct["MOD"] = {
        "MOD_LEVEL": sum_mods_obj.level,
        "MOD_INFO": sum_mods_obj.mod_info,
    }

    sum_res_obj = Residue(sum_res_dct, schema, output_rules, nomenclature)

    return sum_res_obj


if __name__ == "__main__":

    # usr_res_info = {
    #     "link": "",
    #     "c_count": 18,
    #     "db_count": 1,
    #     "db_info_sum": {
    #         "level": 0.2,
    #         "info": {
    #             "0_DB": {
    #                 "count": 1,
    #                 "cv": "",
    #                 "level": 0.2,
    #                 "order": 0,
    #                 "site": ["4"],
    #                 "site_info": ["4E"],
    #             }
    #         },
    #     },
    #     "o_count": 2,
    #     "o_info_sum": {
    #         "level": 5,
    #         "info": {
    #             "0_O": {
    #                 "count": 2,
    #                 "cv": "OH",
    #                 "level": 5,
    #                 "order": 0,
    #                 "site": ["1", "3"],
    #                 "site_info": ["1R", "3S"],
    #             }
    #         },
    #     },
    #     "mod_info_sum": {
    #         "level": 5,
    #         "info": {
    #             "5.01_OH": {
    #                 "count": 2,
    #                 "cv": "OH",
    #                 "level": 5,
    #                 "order": 5.01,
    #                 "site": ["8", "18"],
    #                 "site_info": ["8S", "18R"],
    #                 "verbose": {"elements": {"O": 1}, "mass_shift": 16},
    #             },
    #             "5.02_oxo": {
    #                 "count": 1,
    #                 "cv": "oxo",
    #                 "level": 4,
    #                 "order": 5.02,
    #                 "site": ["10"],
    #                 "site_info": [],
    #                 "verbose": {"elements": {"H": -2, "O": 1}, "mass_shift": 14},
    #             },
    #         },
    #     },
    # }
    usr_res_info = {
        "link": "d",
        "c_count": 18,
        "db_count": 1,
        "db_info_sum": {
            "level": 0.2,
            "info": {
                "0_DB": {
                    "count": 1,
                    "cv": "",
                    "level": 0.2,
                    "order": 0,
                    "site": ["4"],
                    "site_info": ["4E"],
                }
            },
        },
        "o_count": 0,
        "o_info_sum": {},
        "mod_info_sum": {
            "level": 5,
            "info": {
                "5.01_OH": {
                    "count": 2,
                    "cv": "OH",
                    "level": 5,
                    "order": 5.01,
                    "site": ["8", "18"],
                    "site_info": ["8S", "18R"],
                    "verbose": {"elements": {"O": 1}, "mass_shift": 16},
                },
                "5.02_oxo": {
                    "count": 1,
                    "cv": "oxo",
                    "level": 4,
                    "order": 5.02,
                    "site": ["10"],
                    "site_info": [],
                    "verbose": {"elements": {"H": -2, "O": 1}, "mass_shift": 14},
                },
            },
        },
    }
    # for r in usr_res_info:
    #     usr_res_obj = Residue(usr_res_info[r])
    #     logger.debug(usr_res_obj.linked_ids)
    #     # usr_res_json = res_obj.to_json()

    # usr_res_obj = merge_residues(usr_res_info)
    usr_res_obj = Residue(usr_res_info)
    app_logger.debug(usr_res_obj.linked_ids)
    # usr_res_json = res_obj.to_json()
    app_logger.info("FINISHED")
