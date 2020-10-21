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
import re
from typing import Dict, List, Union

from jsonschema import Draft7Validator, RefResolver
from natsort import natsorted

from lynx.models.defaults import (
    core_schema,
    core_schema_path,
    mod_level_lst,
    default_output_rules,
    lynx_schema_cfg,
)
from lynx.utils.basics import get_abs_path
from lynx.utils.cfg_reader import api_version
from lynx.utils.log import app_logger
from lynx.utils.params_loader import load_output_rule
from lynx.utils.toolbox import check_json


class SP_O(object):
    def __init__(
        self,
        sp_o_info_sum: dict,
        schema: str = "lynx_o",
        output_rules: dict = default_output_rules,
        nomenclature: str = "LipidLynxX",
        logger=app_logger,
    ):
        self.logger = logger
        self.nomenclature = nomenclature
        self.export_rule = load_output_rule(output_rules, nomenclature)
        self.sp_o_sites_rule = self.export_rule.get("SITE", None)
        self.sp_o_separators = self.export_rule.get("SEPARATOR", {})
        if not self.sp_o_sites_rule:
            raise ValueError(
                f"Cannot find output rule for 'O_SITES' from nomenclature: {nomenclature}."
            )
        self.info = sp_o_info_sum.get("info", {}).get("0.02_SP_O", {})
        self.schema = schema
        self.type = "SP_O"
        self.level = str(sp_o_info_sum.get("level", 0))
        if self.level == "0":
            self.level = "0.0"
        with open(get_abs_path(lynx_schema_cfg[self.schema]), "r") as s_obj:
            self.validator = Draft7Validator(
                json.load(s_obj),
                resolver=RefResolver(
                    f"file://{core_schema_path}", referrer=core_schema
                ),
            )

        self.count = self.info.get("count", 0)
        self.site = natsorted(self.info.get("site", []))
        self.site_info = natsorted(self.info.get("site_info", []))
        self.details = self.to_dict()

    def __str__(self):
        return self.to_json()

    def __repr__(self):
        return self.to_json()

    def to_sp_o_level(self, level: str = "0") -> str:
        o_str = ""
        site_left = re.sub(
            r"\\", "", self.sp_o_separators.get("SITE_BRACKET_LEFT", "{")
        )
        site_right = re.sub(
            r"\\", "", self.sp_o_separators.get("SITE_BRACKET_RIGHT", "{")
        )
        if not isinstance(level, str):
            level = str(level)
        if float(level) > float(self.level):
            raise ValueError(
                f'Cannot convert to higher level than the sp_o_level "{self.level}". Input:{level}'
            )
        try:
            level_i = int(level)
        except ValueError:
            self.logger.warning(f"Cannot process db.level: {level}")
            level_i = 0
        if level_i > 3:
            sp_o_sites_lst = natsorted(self.info.get("site", []))
            sp_o_sites_info_lst = natsorted(self.info.get("site_info", []))
            sp_o_seg_site_str = ""
            if level == "4":
                sp_o_seg_site_str += ",".join(sp_o_sites_lst)
            elif level == "5":
                if sp_o_sites_info_lst:
                    sp_o_seg_site_str += ",".join(sp_o_sites_info_lst)
                else:
                    sp_o_seg_site_str += ",".join(sp_o_sites_lst)
            else:
                pass
        else:
            sp_o_seg_site_str = ""
        if sp_o_seg_site_str:
            o_str += f"{site_left}{sp_o_seg_site_str}{site_right}"

        return o_str

    def to_all_levels(self, as_list: bool = False) -> Union[Dict[str, str], List[str]]:
        all_levels_dct = {}
        if self.level in mod_level_lst:
            sp_o_idx = mod_level_lst.index(self.level)
            output_levels_lst = mod_level_lst[: sp_o_idx + 1]
        else:
            raise ValueError(f"DB level not supported: {self.level}")

        for level in output_levels_lst:
            all_levels_dct[level] = self.to_sp_o_level(level)
        all_levels_info = all_levels_dct

        return all_levels_info

    def to_dict(self):
        linked_ids = self.to_all_levels()
        sp_o_id = linked_ids.get(self.level, "")
        if float(self.level) >= 0:
            sum_sp_o_info_dct = {
                "api_version": api_version,
                "type": self.type,
                "id": sp_o_id,
                "level": self.level,
                "linked_ids": linked_ids,
                "linked_levels": natsorted(list(linked_ids.keys())),
                "info": self.info,
            }
        else:
            raise ValueError(
                f"Cannot format SP_O code to level {self.level} "
                f"from input: {self.info}"
            )

        return sum_sp_o_info_dct

    def to_json(self):
        sp_o_json_str = json.dumps(self.details)

        if check_json(
            validator=self.validator,
            json_obj=json.loads(sp_o_json_str),
        ):
            return sp_o_json_str
        else:
            raise Exception(f"Schema test FAILED. Schema {self.schema}")


if __name__ == "__main__":

    usr_sp_o = {
        "level": 5,
        "info": {
            "0.02_SP_O": {
                "count": 2,
                "cv": "",
                "level": 5,
                "order": 0.02,
                "site": ["1", "3"],
                "site_info": ["1R", "3S"],
            }
        },
    }

    sp_o_obj = SP_O(sp_o_info_sum=usr_sp_o)
    print(sp_o_obj.to_json())
