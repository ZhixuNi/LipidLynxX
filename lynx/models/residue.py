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
from typing import List

from natsort import natsorted
import regex as re
from jsonschema import Draft7Validator, RefResolver

from lynx.models.db import DB
from lynx.models.defaults import (
    db_level_lst,
    default_output_rules,
    mod_level_lst,
    res_schema,
    res_schema_path,
)
from lynx.models.modification import merge_mods, Modification
from lynx.models.sp_o import SP_O
from lynx.utils.cfg_reader import api_version
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
        self.nomenclature = nomenclature
        self.logger = logger
        self.export_rule = load_output_rule(output_rules, nomenclature)
        self.res_rule = self.export_rule.get("RESIDUE", None)
        self.res_rule_orders = self.res_rule.get("RESIDUE_INFO", {}).get("ORDER", [])
        self.res_separators = self.export_rule.get("SEPARATOR", [])
        self.info = residue_info.get("info", {})
        self.level = residue_info.get("level", "0.0")
        # transfer the m/d/t label in Cer/SM/SP from link section into O section
        self.__replace_mdt__()

        self.db = self.__init_db__()
        self.sp_o = self.__init_sp_o__()
        self.mod = self.__init_mod__()
        self.linked_levels = self.__init_linked_levels__()

        if float(self.mod.level) > 0:
            self.is_modified = True
        else:
            self.is_modified = False

        self.schema = schema
        self.type = "FattyAcid"
        resolver = RefResolver(
            referrer=res_schema, base_uri=f"file://{os.path.dirname(res_schema_path)}/"
        )
        self.validator = Draft7Validator(res_schema, resolver=resolver)

        self.linked_ids = self.to_all_levels()

    def __replace_mdt__(self):
        link = self.info.get("link", "")
        if link:
            existed_o_count = self.info.get("sp_o_count", 0)
            # convert m/d/t into corresponding O section
            if link.lower() == "m":
                self.info["link"] = ""
                self.info["sp_o_count"] = 1 + existed_o_count
            elif link.lower() == "d":
                self.info["link"] = ""
                self.info["sp_o_count"] = 2 + existed_o_count
            elif link.lower() == "t":
                self.info["link"] = ""
                self.info["sp_o_count"] = 3 + existed_o_count
            else:
                pass
        else:
            pass

    def __init_db__(self) -> DB:
        db_info_sum = self.info.get("db_info_sum", {})

        db = DB(db_info_sum, nomenclature=self.nomenclature, logger=self.logger)

        return db

    def __init_sp_o__(self) -> SP_O:
        sp_o_info_sum = self.info.get("sp_o_info_sum", {})

        sp_o = SP_O(sp_o_info_sum, nomenclature=self.nomenclature, logger=self.logger)

        return sp_o

    def __init_mod__(self) -> Modification:
        mod_info_sum = self.info.get("mod_info_sum", {})
        db_count = self.info.get("db_count", 0)
        sp_o_count = self.info.get("sp_o_count", 0)

        mod = Modification(
            mod_info_sum,
            db_count=db_count,
            sp_o_count=sp_o_count,
            nomenclature=self.nomenclature,
            logger=self.logger,
        )
        return mod

    def __init_linked_levels__(self) -> List[str]:
        db_lv = self.db.level
        db_lvs = self.db.details.get("linked_levels", [db_lv])
        sp_o_lv = self.sp_o.level
        sp_o_lvs = self.sp_o.details.get("linked_levels", [sp_o_lv])
        mod_lv = self.mod.level
        mod_lvs = self.mod.details.get("linked_levels", [mod_lv])

        max_level_f = max(float(db_lv), float(sp_o_lv))
        max_level_f += int(mod_lv)
        max_level = f"{max_level_f:.1f}"

        if max_level != self.level:
            max_lv = max(float(max_level), float(self.level))
            self.level = f"{max_lv:.1f}"
        else:
            pass
        sub_lvs = list(set(db_lvs + sp_o_lvs))
        linked_lvs = []
        for main_lv in mod_lvs:
            for sub_lv in sub_lvs:
                linked_lv = int(main_lv) + float(sub_lv)
                linked_lvs.append(f"{linked_lv:.1f}")
        linked_lvs = natsorted(list(set(linked_lvs)))
        return linked_lvs

    def collect_info(self) -> dict:
        formatted_dct = {
            "RESIDUE_SEPARATOR": self.res_separators.get("RESIDUE_SEPARATOR", "_"),
            "LINK": self.info.get("link", ""),
            "C_COUNT": self.info.get("c_count", 0),
            "DB_SEPARATOR": self.res_separators.get("DB_SEPARATOR", ":"),
            "DB_COUNT": self.db.count,
        }
        if self.db.count > 0:
            formatted_dct["DB_INFO_SUM"] = self.db.details.get("linked_ids", {})
        if self.sp_o.count > 0:
            formatted_dct["SP_O_SEPARATOR"] = self.res_separators.get(
                "SP_O_SEPARATOR", ";"
            )
            formatted_dct["SP_O_COUNT"] = self.sp_o.count
            formatted_dct["SP_O_INFO_SUM"] = self.sp_o.details.get("linked_ids", {})
        if self.is_modified:
            formatted_dct["MOD_INFO_SUM"] = self.mod.details.get("linked_ids", {})

        return formatted_dct

    def to_all_levels(self):
        res_str_dct = {}
        collected_info_dct = self.collect_info()

        mod_left = self.res_separators.get("MOD_BRACKET_LEFT", "<")
        mod_right = self.res_separators.get("MOD_BRACKET_RIGHT", ">")
        mod_left = re.sub(r"\\", "", mod_left)
        mod_right = re.sub(r"\\", "", mod_right)

        for lv in self.linked_levels:
            res_str = ""
            lv_seg_lst = lv.split(".")
            if len(lv_seg_lst) == 1 and lv_seg_lst[0] in mod_level_lst:
                main_lv = lv_seg_lst[0]
                sub_lv = "0.0"
            elif len(lv_seg_lst) == 2 and lv_seg_lst[0] in mod_level_lst:
                main_lv = lv_seg_lst[0]
                sub_lv = f"0.{lv_seg_lst[1]}"
                if sub_lv in db_level_lst:
                    pass
                else:
                    sub_lv = "0.0"
            else:
                main_lv = "0"
                sub_lv = "0.0"

            for o in self.res_rule_orders:
                if o in collected_info_dct:
                    o_info = collected_info_dct.get(o)
                    if isinstance(o_info, int) or isinstance(o_info, float):
                        o_info = str(o_info)
                    else:
                        pass
                    if isinstance(o_info, str):
                        if re.match("RESIDUE_SEPARATOR", o):
                            pass
                        else:
                            res_str += o_info
                    else:
                        if re.match(r".*_INFO_SUM$", o) and isinstance(o_info, dict):
                            if o == "MOD_INFO_SUM":
                                mod_seg = o_info.get(main_lv, "")
                                if mod_seg:
                                    res_str += f'{mod_left}{mod_seg}{mod_right}'
                            else:
                                res_str += o_info.get(sub_lv, "")

            na_brackets_lst = [r"\<\>", r"\{\}", r"\[\]", r"\(\)"]
            for b in na_brackets_lst:
                res_str = re.sub(b, "", res_str)
            res_str_dct[lv] = res_str.strip(";")

        return res_str_dct

    def __str__(self):
        return self.to_json()

    def __repr__(self):
        return self.to_json()

    def to_dict(self):
        res_id = self.linked_ids.get(self.level, "")
        if float(self.level) >= 0:
            sum_res_info_dct = {
                "api_version": api_version,
                "type": self.type,
                "id": res_id,
                "level": self.level,
                "linked_ids": self.linked_ids,
                "linked_levels": self.linked_levels,
                "info": self.info,
            }
        else:
            raise ValueError(
                f"Cannot format Residue information to level {self.level} "
                f"from input: {self.info}"
            )

        return sum_res_info_dct

    # def to_json(self):
    #     fa_lite_info_dct = self.fa_info_dct
    #     fa_lite_info_dct.pop("mod_obj", None)
    #     fa_json_str = json.dumps(fa_lite_info_dct)
    #     if check_json(self.validator, json.loads(fa_json_str, logger=self.logger)):
    #         return fa_json_str
    #     else:
    #         raise Exception(f"JSON Schema check FAILED. Schema {self.schema}")


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
        "MOD_INFO": sum_mods_obj.info,
    }

    sum_res_obj = Residue(sum_res_dct, schema, output_rules, nomenclature)

    return sum_res_obj


if __name__ == "__main__":

    usr_res_info = {
        "level": 5.2,
        "info": {
            "link": "",
            "c_count": 18,
            "db_count": 1,
            "db_info_sum": {
                "level": 0.2,
                "info": {
                    "0.01_DB": {
                        "count": 1,
                        "cv": "",
                        "level": 0.2,
                        "order": 0.01,
                        "site": ["4"],
                        "site_info": ["4E"],
                    }
                },
            },
            "sp_o_count": 2,
            "sp_o_info_sum": {
                "level": 0.2,
                "info": {
                    "0.02_SP_O": {
                        "count": 2,
                        "cv": "OH",
                        "level": 0.2,
                        "order": 0.02,
                        "site": ["1", "3"],
                        "site_info": ["1R", "3S"],
                    }
                },
            },
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
        },
    }

    # usr_res_info = {
    #     "level": 5.2,
    #     "info": {
    #         "link": "",
    #         "c_count": 18,
    #         "db_count": 1,
    #         "db_info_sum": {
    #             "level": 0.2,
    #             "info": {
    #                 "0.01_DB": {
    #                     "count": 1,
    #                     "cv": "",
    #                     "level": 0.2,
    #                     "order": 0.01,
    #                     "site": ["4"],
    #                     "site_info": ["4E"],
    #                 }
    #             },
    #         },
    #         "sp_o_count": 0,
    #         "sp_o_info_sum": {
    #             "level": 0,
    #             "info": {
    #                 "0.02_SP_O": {
    #                     "count": 0,
    #                     "cv": "O",
    #                     "level": 0,
    #                     "order": 0.02,
    #                     "site": [],
    #                     "site_info": [],
    #                 }
    #             },
    #         },
    #         "mod_info_sum": {
    #             "level": 5,
    #             "info": {
    #                 "5.01_OH": {
    #                     "count": 2,
    #                     "cv": "OH",
    #                     "level": 5,
    #                     "order": 5.01,
    #                     "site": ["8", "18"],
    #                     "site_info": ["8S", "18R"],
    #                     "verbose": {"elements": {"O": 1}, "mass_shift": 16},
    #                 },
    #                 "5.02_oxo": {
    #                     "count": 1,
    #                     "cv": "oxo",
    #                     "level": 4,
    #                     "order": 5.02,
    #                     "site": ["10"],
    #                     "site_info": [],
    #                     "verbose": {"elements": {"H": -2, "O": 1}, "mass_shift": 14},
    #                 },
    #             },
    #         },
    #     },
    # }

    # for r in usr_res_info:
    #     usr_res_obj = Residue(usr_res_info[r])
    #     logger.debug(usr_res_obj.linked_ids)
    #     # usr_res_json = res_obj.to_json()

    # usr_res_obj = merge_residues(usr_res_info)
    usr_res_obj = Residue(usr_res_info)
    app_logger.debug(usr_res_obj.linked_ids)
    # usr_res_json = res_obj.to_json()
    app_logger.info("FINISHED")
