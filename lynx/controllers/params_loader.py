# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import configparser
import os
import re
from typing import Dict, List, Tuple

import pandas as pd
from natsort import natsorted

from lynx.controllers.general_functions import get_abs_path, load_folder
from lynx.models.rules import InputRules, OutputRules
from lynx.models.log import logger


def load_cfg_info(cfg_path: str = None) -> Dict[str, str]:
    cfg_path_dct = {}
    default_fields = [
        "cv",
        "rules",
        "input_rules",
        "output_rules",
        "mod_cfg",
        "abbr_cfg",
        "base_url",
    ]
    config = configparser.ConfigParser()
    if cfg_path and isinstance(cfg_path, str):
        config_path = get_abs_path(cfg_path)
    else:
        try:
            config_path = get_abs_path("config.ini")
        except FileNotFoundError:
            config_path = get_abs_path("configure.ini")

    config.read(config_path)
    if config.has_section("settings"):
        user_cfg = "settings"
    elif config.has_section("default"):
        user_cfg = "default"
    else:
        user_cfg = ""
        raise ValueError(f"Cannot load settings from file {config_path}")

    if len(user_cfg) > 2:
        options = config.options(user_cfg)
        for field in default_fields:
            if field in options and field != "base_url":
                cfg_path_dct[field] = get_abs_path(config.get(user_cfg, field))

    if "base_url" not in cfg_path_dct:
        cfg_path_dct["base_url"] = r"http://127.0.0.1:5000"

    return cfg_path_dct


def build_parser(rules_file: str) -> Tuple[dict, dict]:
    """
    Read predefined rules from configurations folder and export as a dictionary

    Args:
        rules_file: the path for the rules file

    Returns:
        dict contains the regular expressions as a dict

    """

    rules_df = pd.read_csv(
        rules_file,
        header=0,
        index_col=False,
        skip_blank_lines=True,
        na_filter=True,
        na_values=None,
    )
    rules_df.drop_duplicates(
        subset=["CLASS", "REGULAR_EXPRESSION"], keep="first", inplace=True
    )

    class_rules_dct = {}
    rules_class_dct = {}

    for i, r in rules_df.iterrows():
        if isinstance(r["CLASS"], str):
            rules_str = r["REGULAR_EXPRESSION"].strip('"')
            rules = re.compile(rules_str)

            rules_checker = True
            if isinstance(r["EXAMPLE"], str):
                if not rules.match(r["EXAMPLE"]):
                    rules_checker = False
                    logger.warning(
                        f'Rule example: "{r["EXAMPLE"]}" NOT fit with rule: "{rules_str}" -> skipped...'
                    )

            if rules_checker:
                rules_class_dct[rules] = r["CLASS"]
                if r["CLASS"] not in class_rules_dct:
                    class_rules_dct[r["CLASS"]] = [rules]
                else:
                    class_rules_dct[r["CLASS"]].append(rules)

                logger.debug(
                    f'Rule added: "{r["CLASS"]}" -> "{r["REMARK"]}" -> "{r["EXAMPLE"]}"'
                )

    return class_rules_dct, rules_class_dct


def get_cv_lst(cv_file: str) -> list:
    cv_df = pd.read_csv(cv_file)
    cv_lst = cv_df["CV"].tolist()

    return cv_lst


def build_mod_parser(cv_alias_info: Dict[str, List[str]]) -> dict:
    cv_patterns_dct = {}
    for cv in cv_alias_info:
        alias_lst = cv_alias_info[cv]
        for alia in alias_lst:
            cv_patterns_dct[alia] = re.compile(
                r"(\s*[;_+]\s*)?(?P<FRONT>\d\d?[xX]?)?(?P<MOD>{mod})(?P<END>\d\d?)?(?P<REPLACE>@[CHNOP]\d\d?)?".format(
                    mod=alia
                )
            )

    return cv_patterns_dct


def build_input_rules(folder: str) -> dict:

    input_rules = {}
    file_path_lst = load_folder(folder, file_type=".json")
    logger.debug(f"Fund JSON config files: \n {file_path_lst}")

    for f in file_path_lst:
        temp_rules = InputRules(f)
        idx_lst = [os.path.basename(f)] + temp_rules.source
        idx = "#".join(idx_lst)
        for c in temp_rules.rules:
            c_class_str = temp_rules.rules[c].get("CLASS", "")
            c_lmsd_classes = temp_rules.rules[c].get("LMSD_CLASSES", [])
            existed_c_info = input_rules.get(c_class_str, {})
            c_rgx = existed_c_info.get("SEARCH", re.compile(c_class_str))
            c_pattern = existed_c_info.get("MATCH", {})
            c_lmsd_classes.extend(existed_c_info.get("LMSD_CLASSES", []))
            c_lmsd_classes = natsorted(list(set(c_lmsd_classes)))
            c_info = temp_rules.rules[c]
            del c_info["CLASS"]
            del c_info["LMSD_CLASSES"]
            c_pattern[idx] = c_info

            input_rules[c_class_str] = {
                "LMSD_CLASSES": c_lmsd_classes,
                "SEARCH": c_rgx,
                "MATCH": c_pattern,
            }

    logger.debug(input_rules)

    return input_rules


def build_output_rules(folder: str) -> dict:

    output_rules = {}
    file_path_lst = load_folder(folder, file_type=".json")
    logger.debug(f"Fund JSON config files: \n {file_path_lst}")

    for f in file_path_lst:
        temp_rules = OutputRules(f)
        idx = f"{temp_rules.nomenclature}@{temp_rules.date}"
        output_rules[idx] = temp_rules.rules.get("LMSD_CLASSES")

    logger.debug(output_rules)

    return output_rules