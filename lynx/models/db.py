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
    db_level_lst,
    default_output_rules,
    lynx_schema_cfg,
)
from lynx.utils.basics import get_abs_path
from lynx.utils.cfg_reader import api_version
from lynx.utils.log import app_logger
from lynx.utils.params_loader import load_output_rule
from lynx.utils.toolbox import check_json


class DB(object):
    def __init__(
        self,
        db_info_sum: dict,
        schema: str = "lynx_db",
        output_rules: dict = default_output_rules,
        nomenclature: str = "LipidLynxX",
        logger=app_logger,
    ):
        self.logger = logger
        self.nomenclature = nomenclature
        self.export_rule = load_output_rule(output_rules, nomenclature)
        self.db_sites_rule = self.export_rule.get("SITE", None)
        self.db_separators = self.export_rule.get("SEPARATOR", {})
        if not self.db_sites_rule:
            self.logger.info(
                f"Cannot find output rule for 'DB_SITES' from nomenclature: {nomenclature}."
            )
        self.info = db_info_sum.get("info", {}).get("0.01_DB", {})
        if not self.info:
            self.info = {
                "count": 0,
                "cv": "",
                "level": 0,
                "order": 0.01,
                "site": [],
                "site_info": [],
            }
        self.schema = schema
        self.type = "DB"
        self.level = str(db_info_sum.get("level", 0))
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
        self.linked_ids = self.details.get("linked_ids", {})
        self.linked_levels = self.details.get("linked_levels", [])

    def __str__(self):
        return self.to_json()

    def __repr__(self):
        return self.to_json()

    def to_db_level(self, level: str = "0.0") -> str:
        db_str = ""
        site_left = self.db_separators.get("SITE_BRACKET_LEFT", "{")
        site_right = self.db_separators.get("SITE_BRACKET_RIGHT", "{")
        site_left = re.sub(r"\\", "", site_left)
        site_right = re.sub(r"\\", "", site_right)
        if not isinstance(level, str):
            level = str(level)
        if float(level) > float(self.level):
            raise ValueError(
                f'Cannot convert to higher level than the db_level "{self.level}". Input:{level}'
            )
        try:
            level_f = float(level)
        except ValueError:
            self.logger.warning(f"Cannot process db.level: {level}")
            level_f = 0.0
        if level_f > 0:
            db_sites_lst = natsorted(self.info.get("site", []))
            db_sites_info_lst = natsorted(self.info.get("site_info", []))
            db_seg_site_str = ""
            if level == "0.1":
                db_seg_site_str += ",".join(db_sites_lst)
            elif level == "0.2":
                if db_sites_info_lst:
                    db_seg_site_str += ",".join(db_sites_info_lst)
                else:
                    db_seg_site_str += ",".join(db_sites_lst)
            else:
                pass
        else:
            db_seg_site_str = ""
        if db_seg_site_str:
            db_str += f"{site_left}{db_seg_site_str}{site_right}"

        return db_str

    def to_all_levels(self, as_list: bool = False) -> Union[Dict[str, str], List[str]]:
        all_levels_dct = {}
        if self.level in db_level_lst:
            db_idx = db_level_lst.index(self.level)
            output_levels_lst = db_level_lst[: db_idx + 1]
        else:
            raise ValueError(f"DB level not supported: {self.level}")

        for level in output_levels_lst:
            all_levels_dct[level] = self.to_db_level(level)
        all_levels_info = all_levels_dct

        return all_levels_info

    def to_dict(self):
        linked_ids = self.to_all_levels()
        db_id = linked_ids.get(self.level, "")
        if float(self.level) >= 0:
            sum_db_info_dct = {
                "api_version": api_version,
                "type": self.type,
                "id": db_id,
                "level": self.level,
                "linked_ids": linked_ids,
                "linked_levels": natsorted(list(linked_ids.keys())),
                "info": self.info,
            }
        else:
            raise ValueError(
                f"Cannot format DB code to level {self.level} "
                f"from input: {self.info}"
            )

        return sum_db_info_dct

    def to_json(self):
        db_json_str = json.dumps(self.details)

        if check_json(
            validator=self.validator,
            json_obj=json.loads(db_json_str),
        ):
            return db_json_str
        else:
            raise Exception(f"Schema test FAILED. Schema {self.schema}")


if __name__ == "__main__":

    usr_db = {
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
    }

    # db_obj = DB(db_info_sum=usr_db)
    # print(db_obj.to_json())

    db_obj2 = DB(db_info_sum={})
    print(db_obj2.to_json())
