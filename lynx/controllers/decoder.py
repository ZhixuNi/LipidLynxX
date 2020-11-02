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

from typing import Dict, List, Union

import regex as re
from natsort import natsorted

from lynx.controllers.formatter import Formatter

# from lynx.controllers.linker import get_lmsd_name, get_swiss_name
from lynx.models.alias import Alias
from lynx.models.defaults import default_input_rules
from lynx.utils.log import app_logger


class Decoder(object):
    def __init__(
        self,
        input_style: str = "",
        rules: dict = default_input_rules,
        logger=app_logger,
    ):
        self.input_style = input_style
        self.rules = rules
        self.formatter = Formatter()
        self.alias = Alias(logger=logger)
        self.logger = logger

    def check_segments(self, lipid_name: str, rule_class: str, rule: str):
        c = rule_class
        matched_info_dct = {}
        max_residues = 1
        is_this_class = False
        c_search_rgx = self.rules[c].get("SEARCH", None)
        if c_search_rgx.search(lipid_name):
            is_this_class = True
        else:
            if rule_class in [
                "RESIDUE",
                "RESIDUE_INFO",
                "RESIDUE_INFO_SUM",
                "RESIDUE_ALIAS",
                "ALIAS",
            ]:
                is_this_class = True
            else:
                pass
        c_match_rgx_dct = self.rules[c].get("MATCH", None)
        if is_this_class and isinstance(c_match_rgx_dct, dict):
            if rule:
                if rule in c_match_rgx_dct:
                    m_pattern = c_match_rgx_dct[rule]["MATCH"]
                    max_residues = c_match_rgx_dct[rule].get("MAX_RESIDUES", 1)
                    # m_groups = c_match_rgx_dct[rule]["GROUPS"]  # type: list
                    m_match = m_pattern.match(lipid_name)
                    if m_match:
                        # matched_dct = {}
                        matched_info_dct = m_match.capturesdict()
                else:
                    raise ValueError(f"Can not find rule: {rule} in configuration.")
            else:
                raise ValueError(f"Must provide a {rule} in configuration to search.")

        if any(
            [
                re.search(r"^residue[s]?(_alia[s]|_info)?$", c, re.IGNORECASE),
                re.search(r"mod", c, re.IGNORECASE),
            ]
        ):
            residues_count = 0
            pass
        else:
            residues_count = len(matched_info_dct.get("RESIDUE_INFO", []))
            residues_separator_count = len(
                matched_info_dct.get("RESIDUE_SEPARATOR", [])
            )
            if residues_count == residues_separator_count + 1 and residues_count > 0:
                pass
            else:
                # self.logger.debug(f"The numbers of residues and residue_separators do not fit. Skipped...")
                matched_info_dct = {}

        links = matched_info_dct.get("LINK", None)
        if links:
            if isinstance(links, str):
                if links.lower().startswith(("e", "o", "a")):
                    matched_info_dct["LINK"] = "O-"
                elif links.lower().startswith("p"):
                    matched_info_dct["LINK"] = "P-"
                elif links.lower().startswith("m"):
                    matched_info_dct["LINK"] = "m"
                elif links.lower().startswith("d"):
                    matched_info_dct["LINK"] = "d"
                elif links.lower().startswith("t"):
                    matched_info_dct["LINK"] = "t"
                elif links.lower().startswith("h"):
                    matched_info_dct["LINK"] = "h"
                else:
                    pass
            elif isinstance(links, list):
                std_links = []
                for link in links:
                    if link.lower().startswith(("e", "o", "a")):
                        std_links.append("O-")
                    elif link.lower().startswith("p"):
                        std_links.append("P-")
                    elif link.lower().startswith("m"):
                        std_links.append("m")
                    elif link.lower().startswith("d"):
                        std_links.append("d")
                    elif link.lower().startswith("t"):
                        std_links.append("t")
                    elif link.lower().startswith("h"):
                        std_links.append("t")
                matched_info_dct["LINK"] = std_links
            else:
                pass
        else:
            prefix = matched_info_dct.get("PREFIX", "")
            if prefix:
                if isinstance(prefix, str) and isinstance(links, str):
                    if re.match(r"(O|plasmanyl)[-]?$", prefix, re.IGNORECASE):
                        matched_info_dct["LINK"] = "O-"
                        matched_info_dct["PREFIX"] = ""
                        if residues_count > 0:
                            matched_info_dct["RESIDUE_INFO"] = (
                                "O-" + matched_info_dct["RESIDUE_INFO"]
                            )
                    elif re.match(r"(P|plasmenyl)[-]?$", prefix, re.IGNORECASE):
                        matched_info_dct["LINK"] = "P-"
                        matched_info_dct["PREFIX"] = ""
                        if residues_count > 0:
                            matched_info_dct["RESIDUE_INFO"] = (
                                "P-" + matched_info_dct["RESIDUE_INFO"]
                            )
                    elif re.match(r"((O|plasmanyl)[-]?)L?$", prefix, re.IGNORECASE):
                        matched_info_dct["LINK"] = "O-"
                        matched_info_dct["PREFIX"] = "L"
                        if residues_count > 0:
                            matched_info_dct["RESIDUE_INFO"] = (
                                "O-" + matched_info_dct["RESIDUE_INFO"]
                            )
                    elif re.match(r"((P|plasmenyl)[-]?)L?$", prefix, re.IGNORECASE):
                        matched_info_dct["LINK"] = "P-"
                        matched_info_dct["PREFIX"] = "L"
                        if residues_count > 0:
                            matched_info_dct["RESIDUE_INFO"] = (
                                "P-" + matched_info_dct["RESIDUE_INFO"]
                            )
                elif isinstance(prefix, list) and isinstance(links, list):
                    std_links = []
                    for px in prefix:
                        if re.match(r"(O|plasmanyl)[-]?$", px, re.IGNORECASE):
                            std_links.append("O-")
                            if residues_count > 0:
                                matched_info_dct["RESIDUE_INFO"][0] = (
                                    "O-" + matched_info_dct["RESIDUE_INFO"][0]
                                )
                                matched_info_dct["PREFIX"] = []
                        elif re.match(r"(P|plasmenyl)[-]?$", px, re.IGNORECASE):
                            std_links.append("P-")
                            if residues_count > 0:
                                matched_info_dct["RESIDUE_INFO"][0] = (
                                    "P-" + matched_info_dct["RESIDUE_INFO"][0]
                                )
                                matched_info_dct["PREFIX"] = []
                        elif re.match(r"((O|plasmanyl)[-]?)L?$", px, re.IGNORECASE):
                            std_links.append("O-")
                            if residues_count > 0:
                                matched_info_dct["RESIDUE_INFO"][0] = (
                                    "O-" + matched_info_dct["RESIDUE_INFO"][0]
                                )
                                matched_info_dct["PREFIX"] = ["L"]
                        elif re.match(r"((P|plasmenyl)[-]?)L?$", px, re.IGNORECASE):
                            std_links.append("P-")
                            if residues_count > 0:
                                matched_info_dct["RESIDUE_INFO"][0] = (
                                    "P-" + matched_info_dct["RESIDUE_INFO"][0]
                                )
                                matched_info_dct["PREFIX"] = ["L"]
                    if std_links:
                        matched_info_dct["LINK"] = std_links

        # Check parsed lenth
        rebuild_str_lst = []
        if matched_info_dct:
            for k in matched_info_dct:
                if re.match(
                    r".*((CLASS)|(FIX)|(SUM)|(COUNT)|(SEPARATOR)|(BRACKET)).*",
                    k,
                    re.IGNORECASE,
                ):
                    ks = matched_info_dct[k]
                    if ks:
                        rebuild_str_lst.extend(ks)
            rebuild_str = "".join(rebuild_str_lst)
            if len(rebuild_str) + 4 < len(lipid_name) and re.match(
                r".*[_|/].*", lipid_name, re.IGNORECASE
            ):  # if more than 4 char not parsed, reject matching
                matched_info_dct = {}
        return matched_info_dct

    def check_alias(self, alias: str, alias_type: str = "RESIDUE") -> str:
        defined_id = ""
        if alias_type.upper().startswith("RESIDUE"):
            lite_alias_info = self.alias.residue_alias
        elif alias_type.upper().startswith("LIPID"):
            lite_alias_info = self.alias.lipid_alias
        else:
            raise ValueError(
                f"Cannot load alias_type {alias_type} from defined_alias.json"
            )

        for alias_rgx in lite_alias_info:
            if re.search(alias_rgx, alias):
                defined_id = lite_alias_info[alias_rgx]
            else:
                pass

        # if not defined_id:
        #     self.logger.debug(
        #         f"Cannot decode alias: {alias} using alias_type: {alias_type}."
        #     )

        return defined_id

    def check_residues(
        self,
        rule: str,
        residues: list,
        sum_residues: str,
        alias=None,
        alias_rule: str = "LipidLynxX.json#LipidLynxX",
        max_residues: int = 1,
        separator_levels: dict = None,
        separator: str = "-|/",
        lmsd_classes: List[str] = None,
    ) -> dict:

        if alias is None:
            alias = []
        if separator_levels is None:
            separator_levels = {"B": "", "D": "_", "S": "/"}

        res_lst = residues
        res_sep_lst = re.findall(separator, sum_residues)
        if not res_sep_lst:
            res_sep_lst = [""]
        for s in res_sep_lst:
            res_lst = [res.strip(s) for res in res_lst]
        res_sep_levels = []
        for res_sep in res_sep_lst:
            for lv in separator_levels:
                if res_sep == separator_levels[lv]:
                    res_sep_levels.append(lv)
        if res_sep_levels:
            lv_min = natsorted(res_sep_levels)[0]
        else:
            lv_min = "B"

        out_res_dct = {}
        out_res_lst = []
        res_true_lst = []
        is_fa = False
        for c_lmsd in lmsd_classes:
            if c_lmsd.upper().startswith("FA"):
                is_fa = True
        if not is_fa and "0:0" in res_lst:
            res_true_lst = [res for res in res_lst if res != "0:0"]
        else:
            res_true_lst = res_lst

        if len(res_lst) <= max_residues or len(res_true_lst) <= max_residues:
            for res in res_lst:
                if res in alias:
                    if res != "":
                        def_res = self.check_alias(res, "RESIDUE_INFO")
                        if def_res:
                            res = def_res
                    matched_info_dct = self.check_segments(
                        res, "RESIDUE_INFO", rule=alias_rule
                    )
                else:
                    matched_info_dct = self.check_segments(
                        res, "RESIDUE_INFO", rule=rule
                    )
                # num_o_chk_lst = matched_info_dct.get("NUM_O", [""])
                # if isinstance(num_o_chk_lst, list) and re.match(r'\d?O|O\d?|\d', num_o_chk_lst[0]):
                #     matched_info_dct["MOD_TYPE"] = ["O"] + matched_info_dct.get("MOD_TYPE", [])
                matched_dct = self.formatter.format_residue(matched_info_dct)
                # self.logger.debug(f"{res} matched {rule}:\n {matched_dct}")
                out_res_lst.append(res)
                out_res_dct[res] = matched_dct

        if len(out_res_lst) > 1 and lv_min == "B":
            lv_min = "M"

        if lv_min == "M":
            no_res_lst = []
            spb_lst = []
            o_lst = []
            p_lst = []
            r_lst = []
            for r in out_res_lst:
                if r == "0:0":
                    no_res_lst.append(r)
                elif r.startswith("O-"):
                    o_lst.append(r)
                elif r.startswith("P-"):
                    p_lst.append(r)
                else:
                    r_lst.append(r)

            out_res_lst = (
                natsorted(no_res_lst)
                + natsorted(o_lst)
                + natsorted(p_lst)
                + natsorted(r_lst)
            )

        return {
            "residues_order": out_res_lst,
            "residues_info": out_res_dct,
            # "RESIDUES_SEPARATOR": res_sep_lst,
            "residues_separator_level": lv_min,
        }

    def extract_by_class_rule(
        self, lipid_name: str, c: str, lynx_rule_idx: str = "LipidLynxX.json#LipidLynxX"
    ) -> dict:
        c_lmsd_classes = self.rules[c].get("LMSD_CLASSES", None)

        res_sep = self.rules[c].get("RESIDUES_SEPARATOR", None)  # type: str
        sep_levels = self.rules[c].get("SEPARATOR_LEVELS", {})  # type: dict
        # c_rules = self.rules[c].get("MATCH", {})
        c_rules = self.rules[c].get("MATCH", {})

        matched_info_dct = {}
        for lr in c_rules:
            if re.search(r"Lynx", lr, re.IGNORECASE):
                lynx_rule_idx = lr
        if self.input_style in c_rules:
            use_c_rules = [self.input_style]
        else:
            use_c_rules = c_rules
        for r in use_c_rules:
            matched_dct = self.check_segments(lipid_name, c, r)
            sum_residues_lst = matched_dct.get("RESIDUE_INFO_SUM", [])
            obs_residues_lst = matched_dct.get("RESIDUE_INFO", [])
            alias_lst: list = matched_dct.get("ALIAS", [])
            c_max_res = c_rules[r].get("MAX_RESIDUES", 1)
            if sum_residues_lst and len(sum_residues_lst) == 1 and obs_residues_lst:
                residues_dct = self.check_residues(
                    r,
                    obs_residues_lst,
                    sum_residues_lst[0],
                    alias=alias_lst,
                    alias_rule=lynx_rule_idx,
                    max_residues=c_max_res,
                    separator_levels=sep_levels,
                    separator=res_sep,
                    lmsd_classes=c_lmsd_classes,
                )
                # set specific classes into
                for c_lmsd in c_lmsd_classes:
                    if c_lmsd.upper().startswith("FA"):
                        residues_dct["residues_separator_level"] = "S"
                    elif c_lmsd.upper().startswith("SP"):
                        residues_dct["residues_separator_level"] = "S"
                        sp_res_abbr_lst = residues_dct["residues_order"]
                        if len(sp_res_abbr_lst) == 2:
                            if (
                                ";" in sp_res_abbr_lst[1]
                                and ";" not in sp_res_abbr_lst[0]
                            ):
                                residues_dct["residues_order"] = [
                                    sp_res_abbr_lst[1],
                                    sp_res_abbr_lst[0],
                                ]
                            elif re.match(
                                r"[mdt].*", sp_res_abbr_lst[1]
                            ) and not re.match(r"[mdt].*", sp_res_abbr_lst[0]):
                                residues_dct["residues_order"] = [
                                    sp_res_abbr_lst[1],
                                    sp_res_abbr_lst[0],
                                ]
                    elif c_lmsd.upper().startswith("ST"):
                        residues_dct["residues_separator_level"] = "S"
                matched_info_dct[r] = {
                    "lmsd_classes": c_lmsd_classes,
                    "segments": matched_dct,
                    "residues": residues_dct,
                    # "RESIDUES_SEPARATOR": res_sep,
                    # "SEPARATOR_LEVELS": sep_levels,
                }
            elif sum_residues_lst and len(sum_residues_lst) > 1:
                self.logger.error(
                    f"More than two parts of SUM residues matched: {sum_residues_lst}"
                )
            else:
                pass  # nothing found. the rule is not used.
        return matched_info_dct

    def extract(self, lipid_name: str) -> Dict[str, Union[str, dict]]:

        """
        Main parser to read input abbreviations
        Args:
            lipid_name: input lipid abbreviation to be converted

        Returns:
            extracted_info_dct: parsed information stored as dict

        """

        extracted_info_dct = {}
        force_alias_check = False
        if len(lipid_name) < 16 and re.match(r"^([CDHSA])|(P[DLO]).*", lipid_name):
            force_alias_check = True

        for c in self.rules:
            if lipid_name:
                matched_info_dct = self.extract_by_class_rule(lipid_name, c)
                is_c_count_lst = False
                r_empty_key_lst = []
                for r in matched_info_dct:
                    r_info = matched_info_dct.get(r, {})
                    c_count_lst = r_info.get("segments", {}).get("C_COUNT", [])
                    if c_count_lst:
                        is_c_count_lst = True
                    else:
                        r_empty_key_lst.append(r)
                if r_empty_key_lst:
                    for r_k in r_empty_key_lst:
                        del matched_info_dct[r_k]
                if is_c_count_lst and matched_info_dct:
                    extracted_info_dct[c] = matched_info_dct
            else:
                self.logger.warning(
                    f"No lipid name is given. Please submit a lipid name."
                )

        obs_alias_lst = []
        if force_alias_check or not extracted_info_dct:

            def_alias = self.check_alias(lipid_name, "LIPID")
            if def_alias:
                obs_alias_lst.append(def_alias)
                for cx in self.rules:
                    alias_matched_info_dct = self.extract_by_class_rule(def_alias, cx)
                    is_x_c_count_lst = False
                    rx_empty_key_lst = []
                    for rx in alias_matched_info_dct:
                        rx_info = alias_matched_info_dct.get(rx, {})
                        alias_c_count_lst = rx_info.get("segments", {}).get(
                            "C_COUNT", []
                        )
                        if alias_c_count_lst:
                            is_x_c_count_lst = True
                        else:
                            rx_empty_key_lst.append(rx)
                    if rx_empty_key_lst:
                        for rx_k in rx_empty_key_lst:
                            del alias_matched_info_dct[rx_k]
                    if is_x_c_count_lst and alias_matched_info_dct:
                        extracted_info_dct[cx] = alias_matched_info_dct
            obs_alias_lst = list(set(obs_alias_lst))
        if obs_alias_lst:
            self.logger.debug(f"Using alias: {lipid_name} -> {obs_alias_lst}")

        return extracted_info_dct


if __name__ == "__main__":
    # LIPID MAPS
    # t_in = "GM3(d18:1/18:2(9Z,12Z))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5Z,8Z,11Z,14Z)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5,8,11,14)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/20:4(5,8,11,14)(7R-OH,12S-OH))"
    # t_in = "TG (P-18:1/18:2(9Z,12Z)/5S,15R-DiHETE)"

    # MS-DIAL
    # t_in = "TG(16:0/18:2/20:4<OH>)"
    # t_in = "TG(16:0/18:2/HETE)"
    # RefMet
    t_in = "PC(36:2e)"

    extractor = Decoder(rules=default_input_rules)
    t_out = extractor.extract(t_in)

    app_logger.info(t_out)
    app_logger.info("FIN")
