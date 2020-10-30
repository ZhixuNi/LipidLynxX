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
from typing import Dict, List

from natsort import natsorted

from lynx.controllers.decoder import Decoder
from lynx.models.defaults import (
    default_input_rules,
    default_output_rules,
    supported_levels,
)
from lynx.models.residue import merge_residues, Residue
from lynx.utils.log import app_logger
from lynx.utils.params_loader import load_output_rule


class Encoder(object):
    def __init__(
        self,
        style: str = "LipidLynxX",
        input_rules: dict = default_input_rules,
        output_rules: dict = default_output_rules,
        input_style: str = "",
        logger=app_logger,
    ):
        self.export_style = style
        self.input_rule = input_style
        self.output_rules = load_output_rule(output_rules, style)
        self.class_rules = self.output_rules.get("LMSD_CLASSES", {})
        self.mod_rules = self.output_rules.get("MOD", {}).get("MOD_INFO", {})
        self.sum_mod_rules = self.output_rules.get("MOD", {}).get("MOD_INFO_SUM", {})
        self.residue_rules = self.output_rules.get("RESIDUE", {}).get(
            "RESIDUE_INFO_SUM", {}
        )
        self.separator_levels = self.output_rules.get("SEPARATOR", {}).get(
            "SEPARATOR_LEVELS", {}
        )
        self.separators = self.output_rules.get("SEPARATOR", {})
        self.extractor = Decoder(
            input_style=self.input_rule, rules=input_rules, logger=logger
        )
        self.logger = logger

    def get_best_id(self, candidate: Dict[str, str]) -> str:
        c_lv_lst = natsorted(list(candidate.keys()))
        self.logger.debug(f"Export MAX level: {c_lv_lst[-1]}. ")
        c_max_str = candidate.get(c_lv_lst[-1], "")
        return c_max_str

    def get_best_id_series(self, candidates: List[dict]) -> [dict, str]:
        best_id_dct = {}
        best_id_score_dct = {}
        num_lv = 0
        max_str = ""
        max_sum_c_count = 0
        max_sum_sp_o_count = 0
        best_input_rule = ""
        is_modified = False
        for c_info_set in candidates:
            c_info = c_info_set.get("compiled_names", {})
            c_sum_c = c_info_set.get("sum_c", 0)
            c_sum_sp_o = c_info_set.get("sp_o_count", 0)
            c_in_rule = c_info_set.get("input_rule", 0)
            c_is_modified = c_info_set.get("is_modified", False)
            lmsd_classes = c_info_set.get("lmsd_classes", False)
            best_rule_score = 0
            if c_is_modified:
                is_modified = True
            c_info_max_c_lv = ""
            c_info_max_c = ""
            for c_lv in c_info:
                ci_lv_lst = natsorted(list(c_info[c_lv].keys()))
                if len(c_info[c_lv].get(ci_lv_lst[-1], "")) >= len(c_info_max_c_lv) and len(c_lv) >= len(c_info_max_c):
                    c_info_max_c_lv = c_lv
                    c_info_max_c = c_lv
            if c_info_max_c:
                c_info = {c_info_max_c: c_info[c_info_max_c]}
            for c in c_info:
                c_lv_lst = natsorted(list(c_info[c].keys()))
                c_num_lv = len(c_lv_lst)
                c_max_str = c_info[c].get(c_lv_lst[-1], "")
                if c_num_lv > num_lv:
                    max_str = c_max_str
                    num_lv = c_num_lv
                    max_sum_c_count = max(c_sum_c, max_sum_c_count)
                    max_sum_sp_o_count = max(c_sum_sp_o, max_sum_sp_o_count)
                    best_rule_score = 5
                else:
                    c_max_str = c_info[c].get(c_lv_lst[-1], "")
                    if len(c_max_str) > len(max_str):
                        max_str = c_max_str
                        num_lv = c_num_lv
                        max_sum_c_count = max(c_sum_c, max_sum_c_count)
                        max_sum_sp_o_count = max(c_sum_sp_o, max_sum_sp_o_count)
                        best_rule_score = 4
                    else:
                        if c_sum_c >= max_sum_c_count and c_sum_sp_o >= max_sum_sp_o_count:
                            max_str = c_max_str
                            num_lv = c_num_lv
                            max_sum_c_count = c_sum_c
                            max_sum_sp_o_count = c_sum_sp_o
                            best_rule_score = 3
                        elif c_sum_c >= max_sum_c_count and c_sum_sp_o < max_sum_sp_o_count:
                            max_str = c_max_str
                            num_lv = c_num_lv
                            max_sum_c_count = c_sum_c
                            best_rule_score = 2
                        elif c_sum_c < max_sum_c_count and c_sum_sp_o >= max_sum_sp_o_count:
                            max_sum_sp_o_count = c_sum_sp_o
                            best_rule_score = 1
                        else:
                            pass
                if best_rule_score:
                    best_id_score_dct[best_rule_score] = c_info[c]
            if re.match(r"BioPAN", self.export_style, re.IGNORECASE):
                if is_modified:
                    return {}, {}
            else:
                is_sp_class = False
                for lmsd in lmsd_classes:
                    if re.match(r'^SP.*$', lmsd, re.IGNORECASE):
                        is_sp_class = True
                if is_sp_class and max_sum_sp_o_count > 0:
                    best_rule_score = 1
        if best_id_score_dct:
            max_best_score = max(list(best_id_score_dct.keys()))
            best_id_dct = best_id_score_dct.get(max_best_score, {})
        # add levels for B0, D0, S0 lipids
        if best_id_dct and all(
            [re.match(r"^[BMS]0(.[12])?$", lv) for lv in best_id_dct]
        ):
            for add_lv in supported_levels:
                from_lv = f"{add_lv[0]}0{add_lv[2:]}"
                if from_lv in best_id_dct:
                    best_id_dct[add_lv] = best_id_dct.get(from_lv, "")
        # add levels for lipids with no DB
        updated_best_id_dct = {}
        if (
            best_id_dct
            and max_sum_c_count == 0
            and all([re.match(r"^[BMS][0-5]$", lv) for lv in best_id_dct])
        ):
            for from_lv in best_id_dct:
                updated_best_id_dct[from_lv] = best_id_dct.get(from_lv, "")
                for db_lv in [".1", ".2"]:
                    add_lv = f"{from_lv}{db_lv}"
                    if add_lv in supported_levels and add_lv not in best_id_dct:
                        updated_best_id_dct[add_lv] = best_id_dct.get(from_lv, "")
        if updated_best_id_dct and len(list(updated_best_id_dct.keys())) > len(
            list(best_id_dct.keys())
        ):
            best_id_dct = updated_best_id_dct

        if best_id_dct.get("M1") and not best_id_dct.get("B1"):
            best_id_dct["B1"] = best_id_dct["M1"]
        if best_id_dct.get("M2") and not best_id_dct.get("B2"):
            best_id_dct["B2"] = best_id_dct["M2"]

        return best_id_dct, best_input_rule

    # def check_rest(self, segment_text: str, segment_name: str, lmsd_class: str):
    #     patterns_dct = self.class_rules[lmsd_class].get(segment_name)
    #     out_seg_lst = []
    #     if segment_text and patterns_dct:
    #         for s_rgx in patterns_dct:
    #             self.logger.debug(
    #                 f"Test {segment_text} on {segment_name} of {lmsd_class} using {s_rgx}"
    #             )
    #             s_matched = s_rgx.match(segment_text)
    #             if s_matched:
    #                 defined_seg = patterns_dct[s_rgx]
    #                 if defined_seg == "EXCEPTIONS":
    #                     if lmsd_class in ["GP12"]:
    #                         out_seg_lst.append(segment_text)
    #                     if lmsd_class in ["SP05", "SP06"]:
    #                         out_seg_lst.append(segment_text)
    #                 else:
    #                     out_seg_lst.append(defined_seg)
    #
    #         out_seg_lst = list(filter(None, list(set(out_seg_lst))))
    #
    #     return self.get_best_candidate(out_seg_lst)

    def get_residues(self, residues: dict):
        residues_order = residues.get("residues_order", [])
        residues_sep_level = residues.get("residues_separator_level", "")
        residues_info = residues.get("residues_info", [])
        res_count = len(residues_order)
        # sum_residues_str = ""
        res_lv_id_dct = {}
        res_lv_dct = {}
        sum_lv_lst = []
        for res_abbr in residues_info:
            res_obj = Residue(residues_info[res_abbr])
            res_lv_id_dct[res_abbr] = res_obj.linked_ids
            res_lv_dct[res_abbr] = list(res_obj.linked_ids.keys())
            sum_lv_lst.extend(res_lv_dct[res_abbr])

        sum_lv_lst = natsorted(set(sum_lv_lst))
        for res in res_lv_id_dct:
            r_lv_dct = res_lv_id_dct[res]
            if list(r_lv_dct.keys()) == ["0"]:
                for lv in sum_lv_lst:
                    r_lv_dct[lv] = r_lv_dct["0"]
        res_lv_id_lst_dct = {}
        for sum_lv in sum_lv_lst:
            lv_id_lst = []
            for r in residues_order:
                r_lv_id = res_lv_id_dct.get(r, {}).get(sum_lv, None)
                if r_lv_id:
                    lv_id_lst.append(r_lv_id)
            if len(lv_id_lst) == res_count:
                res_lv_id_lst_dct[sum_lv] = lv_id_lst

        sum_res_id_lv_dct = {}
        sum_res_sep_lv_lst = []

        # set FA or class with one residue into level s
        if len(residues_order) == 1 and len(sum_lv_lst) > 0:
            if "0.1" in sum_lv_lst and residues_sep_level == "B":
                residues_sep_level = "S"

        if residues_sep_level == "S":
            sum_res_sep_lv_lst = ["B", "M", "S"]
        elif residues_sep_level == "M":
            sum_res_sep_lv_lst = ["B", "M"]
        elif residues_sep_level == "B":
            sum_res_sep_lv_lst = ["B"]

        for sep_lv in sum_res_sep_lv_lst:
            if sep_lv == "B":
                # prepare bulk level
                if len(residues_order) > 1:
                    merged_res_obj = merge_residues(residues_order, residues_info)
                    merged_res_linked_ids = merged_res_obj.linked_ids
                    merged_res_lv_lst = list(merged_res_obj.linked_ids.keys())
                elif len(residues_order) == 1:
                    merged_res_obj = Residue(residues_info.get(residues_order[0]))
                    merged_res_linked_ids = merged_res_obj.linked_ids
                    merged_res_lv_lst = list(merged_res_obj.linked_ids.keys())
                else:
                    merged_res_linked_ids = []
                    merged_res_lv_lst = []
                for merged_res_lv in merged_res_lv_lst:
                    sum_res_id_lv_dct[f"B{merged_res_lv}"] = merged_res_linked_ids[
                        merged_res_lv
                    ]
            else:
                for res_lv in res_lv_id_lst_dct:
                    sum_res_id_lv_dct[
                        f"{sep_lv}{res_lv}"
                    ] = f"{self.separator_levels[sep_lv]}".join(
                        res_lv_id_lst_dct[res_lv]
                    )

        return sum_res_id_lv_dct

    @staticmethod
    def check_head_seg(head_segments):
        has_segments = False
        head_seg = ""
        if len(head_segments) == 1 and head_segments[0]:
            head_seg = head_segments[0]
            has_segments = True
        elif len(head_segments) > 1:
            head_seg = "".join(list(set(head_segments)))
            if head_seg == "LO-":
                head_seg = "O-L"
            elif head_seg == "LP-":
                head_seg = "P-L"
            has_segments = True

        return has_segments, head_seg

    @staticmethod
    def check_biopan(lmsd_classes: list, c_prefix_lst: list, residues: dict):
        is_sp_class = False
        is_gl_class = False
        is_gp_class = False
        is_modified = False
        for c in lmsd_classes:
            if c.upper().startswith("GL"):
                is_gl_class = True
                break
            elif c.upper().startswith("GP"):
                is_gp_class = True
                break
            elif c.upper().startswith("SP"):
                is_sp_class = True
                break

        if is_gl_class or is_gp_class or is_sp_class:
            residues_order = residues.get("residues_order", [])
            residues_info = residues.get("residues_info", {})
            for res in residues_order:
                mod_level = (
                    residues_info.get(res, {})
                    .get("info", 0)
                    .get("mod_info_sum", {})
                    .get("level", 0)
                )
                mod_info = (
                    residues_info.get(res, {})
                    .get("info", 0)
                    .get("mod_info_sum", {})
                    .get("info", {})
                )
                if float(mod_level) > 0 or mod_info:
                    is_modified = True

            if is_modified:
                if len(residues_order) > 1:
                    residues = {}
            else:
                if is_gl_class or is_gp_class:
                    residues_separator_level = residues_info.get(
                        "residues_separator_level", "B"
                    )
                    for res in residues_order:
                        res_info = residues_info.get(res, {}).get("info", 0)
                        res_link = res_info.get("link", "")

                        if res_link == "O-":
                            c_prefix_lst.append("O-")
                            residues_info[res]["info"]["link"] = ""
                        elif res_link == "P-":
                            c_prefix_lst.append("O-")
                            residues_info[res]["info"]["link"] = ""
                            residues_info[res]["info"]["db_count"] = (
                                1 + residues_info[res]["info"]["db_count"]
                            )
                            residues_info[res]["info"]["db_info_sum"]["info"][
                                "0.01_DB"
                            ]["count"] = residues_info[res]["info"]["db_count"]
                            residues_info[res]["info"]["db_info_sum"]["info"][
                                "0.01_DB"
                            ]["site"] = []
                            residues_info[res]["info"]["db_info_sum"]["info"][
                                "0.01_DB"
                            ]["site_info"] = []
                    residues = {
                        "residues_order": residues_order,
                        "residues_info": residues_info,
                        "residues_separator_level": residues_separator_level,
                    }
                elif is_sp_class:
                    for res in residues_order:
                        res_info = residues_info.get(res, {}).get("info", 0)
                        res_c_count = res_info.get("c_count", 0)
                        res_sp_o_count = res_info.get("sp_o_count", 0)
                        residues_separator_level = residues_info.get(
                            "residues_separator_level", "B"
                        )
                        res_link = res_info.get("link", "")
                        if (
                            res_sp_o_count == 2
                            or res_link == "d"
                            or res.lower().startswith("d")
                        ):
                            res_db_count = res_info.get("db_count", 0)
                            if res_c_count == 18 and res_db_count in [0, 1]:
                                if len(residues_order) == 2:
                                    residues_order.remove(res)
                                    del residues_info[res]
                                residues = {
                                    "residues_order": residues_order,
                                    "residues_info": residues_info,
                                    "residues_separator_level": residues_separator_level,
                                }
                                if res_db_count == 0:
                                    c_prefix_lst.append("dh")
                            else:
                                residues = {}
                        elif (
                            res_sp_o_count == 0
                            and res_c_count < 27
                            and len(residues_order) == 1
                        ):
                            residues = {
                                "residues_order": residues_order,
                                "residues_info": residues_info,
                                "residues_separator_level": "B",
                            }
                        else:
                            residues = {}
        else:
            # residues_order = residues.get("residues_order", [])
            # residues_info = residues.get("residues_info", {})
            # is_modified = False
            # for res in residues_order:
            #     mod_level = (
            #         residues_info.get(res, {})
            #         .get("info", 0)
            #         .get("mod_info_sum", {})
            #         .get("level", 0)
            #     )
            #     mod_info = (
            #         residues_info.get(res, {})
            #         .get("info", 0)
            #         .get("mod_info_sum", {})
            #         .get("info", {})
            #     )
            #     if float(mod_level) > 0 or mod_info:
            #         is_modified = True
            #
            # if is_modified:
            #     residues = {}
            pass
        return c_prefix_lst, residues, is_modified

    def check_segments(self, parsed_info: dict):
        segments_dct = {}
        lmsd_classes = parsed_info.get("lmsd_classes", [])
        segments = parsed_info.get("segments", {})
        c_prefix_lst = segments.get("PREFIX", [])
        c_suffix_lst = segments.get("SUFFIX", [])
        residues = parsed_info.get("residues", {})
        if re.match(r"BioPAN", self.export_style, re.IGNORECASE):
            c_prefix_lst, residues, is_modified = self.check_biopan(
                lmsd_classes, c_prefix_lst, residues
            )
        c_has_prefix, c_prefix_seg = self.check_head_seg(c_prefix_lst)
        c_has_suffix, c_suffix_seg = self.check_head_seg(c_suffix_lst)
        if residues:
            sum_res_id_lv_dct = self.get_residues(residues)
            obs_c_seg_lst = segments.get("CLASS", [])
            c_seg = ""
            if obs_c_seg_lst and len(obs_c_seg_lst) == 1:
                obs_c_seg = obs_c_seg_lst[0]
                for c in lmsd_classes:
                    c_segments_dct = {}
                    c_orders = []
                    if c in self.class_rules:
                        c_orders = self.class_rules[c].get("ORDER", [])
                        c_identifier_dct = self.class_rules[c].get("CLASS", {})
                        for c_identifier in c_identifier_dct:
                            if c_identifier.match(obs_c_seg):
                                c_seg = c_identifier_dct.get(c_identifier, "")
                                for lv in sum_res_id_lv_dct:
                                    c_lv_segments = {
                                        "CLASS": c_seg,
                                        "RESIDUE_INFO_SUM": sum_res_id_lv_dct[lv],
                                    }
                                    if c_has_prefix:
                                        c_lv_segments["PREFIX"] = c_prefix_seg
                                    if c_has_suffix:
                                        c_lv_segments["SUFFIX"] = c_suffix_seg
                                    c_segments_dct[lv] = c_lv_segments
                            else:
                                pass
                    else:
                        pass
                    if c_segments_dct and c_orders:
                        segments_dct[c] = {"ORDER": c_orders, "INFO": c_segments_dct}
            else:
                self.logger.warning(f"No Class identified!")
        else:
            segments_dct = {}

        return segments_dct

    def compile_segments(self, segments: dict):
        comp_seg_dct = {}
        for c in segments:
            c_lv_id_dct = {}
            c_seg_dct = segments[c]
            c_seg_order = c_seg_dct.get("ORDER", [])
            c_comp_seg_dct = c_seg_dct.get("INFO", {})
            c_optional_seg = self.class_rules.get(c, {}).get("OPTIONAL", [])
            for lv in c_comp_seg_dct:
                lv_seg_info = c_comp_seg_dct[lv]
                lv_seg_lst = []
                for c_seg in c_seg_order:
                    if c_seg in lv_seg_info:
                        lv_seg_lst.append(lv_seg_info[c_seg])
                    elif c_seg in self.separators and c_seg != "SEPARATOR_LEVELS":
                        lv_seg_lst.append(re.sub(r"\\", "", self.separators[c_seg]))
                    else:
                        if c_seg not in c_optional_seg:
                            self.logger.debug(
                                f"Segments not found: {c_seg}, defined orders {c_seg_order}"
                            )
                        else:
                            pass
                c_lv_id_dct[lv] = "".join(lv_seg_lst)
            comp_seg_dct[c] = c_lv_id_dct

        return comp_seg_dct

    def export_all_levels(self, lipid_name: str) -> [dict, str]:

        extracted_info = self.extractor.extract(lipid_name)
        export_info = []
        if extracted_info:
            for p in extracted_info:
                p_info = extracted_info[p]
                self.logger.debug(p_info)
                for in_r in p_info:
                    r_info = p_info[in_r]  # type: dict
                    checked_seg_info = self.check_segments(r_info)
                    if checked_seg_info:
                        comp_dct = self.compile_segments(checked_seg_info)
                        res_info = r_info.get("residues", {}).get("residues_info", {})
                        sum_c = 0
                        sum_sp_o = 0
                        is_modified = False
                        for res in res_info:
                            r_info_dct = res_info[res].get("info", {})
                            sum_c += r_info_dct.get("c_count", 0)
                            sum_sp_o += r_info_dct.get("sp_o_count", 0)
                            mod_info_sum = r_info_dct.get("mod_info_sum", {})
                            mod_level = mod_info_sum.get("level", 0)
                            mod_info = mod_info_sum.get("info", {})
                            if float(mod_level) > 0 or mod_info:
                                is_modified = True
                        export_info.append(
                            {
                                "compiled_names": comp_dct,
                                "sum_c": sum_c,
                                "sp_o_count": sum_sp_o,
                                "input_rule": in_r,
                                "is_modified": is_modified,
                                "segments": r_info.get("segments", {}),
                                "lmsd_classes": r_info.get("lmsd_classes", {}),
                            }
                        )
                    else:
                        pass
            pre_best_export_dct, best_input_rule = self.get_best_id_series(export_info)
            # sort dict by keys
            best_export_dct = {
                k: pre_best_export_dct[k] for k in sorted(pre_best_export_dct)
            }
            self.logger.debug(f"Convert Lipid: {lipid_name} into:\n{best_export_dct}")
        else:
            best_export_dct = {}
            best_input_rule = ""

        return best_export_dct, best_input_rule

    def get_best_rule(self, lipid_name: str) -> str:
        pre_best_export_dct, best_input_rule = self.export_all_levels(lipid_name)
        return best_input_rule

    def convert(self, lipid_name: str, level: str = None) -> str:
        if level in supported_levels:
            best_id = self.export_level(lipid_name, level=level)
        else:
            all_lv_id_dct, best_in_rule = self.export_all_levels(lipid_name)
            best_id = ""
            if all_lv_id_dct:
                best_id = self.get_best_id(all_lv_id_dct)
            else:
                pass
        return best_id

    def export_level(
        self,
        lipid_name: str,
        level: str = "B1",
    ):

        lv_id = ""
        all_lv_id_dct, best_in_rule = self.export_all_levels(lipid_name)
        if level in supported_levels:
            if level in all_lv_id_dct:
                lv_id = all_lv_id_dct[level]
            elif level.upper() == "MAX":
                max_level = natsorted(list(all_lv_id_dct.keys()))[-1]
                self.logger.info(
                    f"Lipid: {lipid_name} cannot be converted into max level: {max_level}. "
                )
                lv_id = all_lv_id_dct[max_level]
            else:
                self.logger.warning(
                    f"Lipid: {lipid_name} cannot be converted into level: {level}. "
                    f"Can be converted into: {all_lv_id_dct}"
                )
        else:
            self.logger.warning(
                f"Level: {level} not supported. Supported levels: {supported_levels}"
            )

        return lv_id

    def export_levels(
        self,
        lipid_name: str,
        levels: list = None,
        import_rules: dict = default_input_rules,
    ) -> dict:

        if levels is None:
            levels = ["B0"]
        lv_id_dct = {}
        all_lv_id_dct, best_in_rule = self.export_all_levels(lipid_name)
        for level in levels:
            if level in supported_levels:
                if level in all_lv_id_dct:
                    lv_id_dct[level] = all_lv_id_dct[level]
                elif level.upper() == "MAX":
                    max_level = natsorted(list(all_lv_id_dct.keys()))[-1]
                    self.logger.info(
                        f"Lipid: {lipid_name} cannot be converted into max level: {max_level}. "
                    )
                    lv_id_dct[max_level] = all_lv_id_dct[max_level]
                    lv_id_dct[level] = all_lv_id_dct[max_level]
                else:
                    raise ValueError(
                        f"Lipid: {lipid_name} cannot be converted into level: {level}. "
                        f"Can be converted into: {all_lv_id_dct}"
                    )
            else:
                raise ValueError(
                    f"Level: {level} not supported. Supported levels: {supported_levels}"
                )

        return lv_id_dct


if __name__ == "__main__":
    t_in_lst = [
        # "GM3(d18:1/18:2(9Z,11Z)(12OH))",
        # "TG P-18:1_18:2(9Z,11Z)(12OH)_18:1(9)(11OH)",
        # "CL(1'-[18:1(9Z)/18:2(9Z,12Z)],3'-[18:2(9Z,12Z)/18:2(9Z,12Z)])",
        # "TG(16:0/18:2/PA)",
        # "PE O-p 32:1",
        # "PE O-a 36:2",
        # "PE O-18:1a/18:1",
        # "PE O-p 36:2",
        # "PE O-18:1p/18:1",
        # "PE-C16:0-C18:2",
        # "HETE",
        # "TG(16:0/18:2/18:2[2xDB,1xOH])",
        # "PC 16:0/18:2[9,12]",
        # "DG 16:2[9,12]_O-18:2_0:0",
        # "14,15-HxB3 (13R)",
        # "C22:5 CE",
        # "15-Keto-PGF2α",
        # "PGF2α",
        # "8-iso PGF2a III",
        # "palmitoleic acid",
        "FA 16:1(9Z)"
    ]
    lynx_gen = Encoder(style="ShorthandNotation")
    for t_in in t_in_lst:
        t1_out = lynx_gen.convert(t_in)
        app_logger.info(f"Input: {t_in} -> Best Output: {t1_out}")
        # t_lv = "B0"
        # t2_out = lynx_gen.export_level(
        #     t_in, level=t_lv, import_rules=default_input_rules
        # )
        # self.logger.info(f"Input: {t_in} -> Output @ Lv {t_lv}: {t2_out}")
        # t_lv_lst = ["B0"]
        # t3_out = lynx_gen.export_levels(
        #     t_in, levels=t_lv_lst, import_rules=default_input_rules
        # )
        # self.logger.info(f"Input: {t_in} -> Output @ Lv {t_lv_lst}: {t3_out}")
        t4_out = lynx_gen.export_all_levels(t_in)
        app_logger.info(f"Input: {t_in} -> Output @ all levels: {t4_out}")

    app_logger.info("fin")
