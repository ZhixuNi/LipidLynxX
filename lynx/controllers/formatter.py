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
import os

from natsort import natsorted

from lynx.models.cv import CV
from lynx.models.defaults import default_cv_file, elem_nominal_info
from lynx.utils.log import app_logger


# upper case keys are used in dicts from regex matching of input or formatted for output
# lower case keys are LipidLynxX internal keys for processing


class Formatter(object):
    def __init__(self, cv_file: str = default_cv_file, logger=app_logger):
        if os.path.isfile(cv_file):
            pass
        else:
            cv_file = default_cv_file
        self.alias2cv = CV(cv_file).info
        self.raw_cv = CV(cv_file).raw_cv
        self.logger = logger

    @staticmethod
    def to_mass_shift(elements: dict) -> int:

        delta = 0
        for elem in elements:
            if elem in elem_nominal_info:
                delta += elem_nominal_info[elem] * elements[elem]
        return delta

    def format_site_info(self, info: str) -> dict:
        info = re.sub(r"[(\[{)\]}]", "", info)
        info_lst = info.split(",")

        site_lst = []
        site_info_lst = []
        for site_seg in info_lst:
            if re.match(r"^\d{1,2}[EZRS]$", site_seg):
                site_lst.append(site_seg[:-1])
                site_info_lst.append(site_seg)
            elif re.match(r"^\d{1,2}$", site_seg):
                site_lst.append(site_seg)
                site_info_lst.append(site_seg)
            else:
                if info:
                    self.logger.debug(f"Can NOT decode site: {site_seg} in {info}")
        chk_site_info_lst = [
            si for si in site_info_lst if re.match(r"^\d{1,2}[EZRS]$", si)
        ]
        if chk_site_info_lst:
            pass
        else:
            site_info_lst = []
        site_info_dct = {"site": site_lst, "site_info": site_info_lst}

        return site_info_dct

    def format_sp_o(self, info: dict):
        o_count = 0
        o_str = ""
        o_lv = 0
        o_site_lst = []
        o_site_info_lst = []
        o_info_lst = info.get("SP_O_COUNT", [])
        o_info_sum_lst = info.get("SP_O_INFO_SUM", [])
        if o_info_lst:
            o_str = str(o_info_lst[0]).strip(" ").upper()
            if o_str:
                if o_str in ["O", "OH"]:
                    o_count = 1
                else:
                    o_count_str = o_str.strip("OH")
                    o_count_str = o_count_str.strip("O")
                    try:
                        o_count = int(o_count_str)
                    except ValueError:
                        pass
            else:
                pass
            if o_count > 0:
                o_site_str_lst = info.get("SP_O_SITE", [])
                if o_site_str_lst:
                    o_site_str = o_site_str_lst[0]
                else:
                    o_site_str = ""
                o_site_info_dct = self.format_site_info(o_site_str)
                o_site_lst = natsorted(o_site_info_dct.get("site", []))
                o_site_info_lst = natsorted(o_site_info_dct.get("site_info", []))
                if re.match(r".*OH.*", o_str):
                    if o_site_info_lst:
                        o_lv = 0.2
                    else:
                        if o_site_lst:
                            o_lv = 0.1
                        else:
                            o_lv = 0.0
                elif re.match(r"^\d{0,2}O\d{0,2}$", o_str):
                    o_lv = 2
        info_sum_o_count = 0
        if len(o_info_sum_lst) == 1:
            info_sum = o_info_sum_lst[0]
            if isinstance(info_sum, str) and len(info_sum) < 4:
                info_sum_o_count_str = info_sum.strip(";")
                info_sum_o_count_str = info_sum_o_count_str.strip("OH")
                info_sum_o_count_str = info_sum_o_count_str.strip("O")
                try:
                    info_sum_o_count += int(info_sum_o_count_str)
                except ValueError:
                    pass
        if 0 < info_sum_o_count < 3 and o_count in [0, 1]:
            o_count = info_sum_o_count
            o_lv = 0
            o_site_lst = []
            o_site_info_lst = []

        if o_lv > 0:
            sp_o_cv = "OH"
        else:
            sp_o_cv = "O"
        o_info = {
            "count": o_count,
            "cv": sp_o_cv,
            "level": o_lv,
            "order": 0.02,
            "site": o_site_lst,
            "site_info": o_site_info_lst,
        }

        return o_info

    def format_db(self, info: dict) -> dict:
        db_lv = 0
        db_site_lst = []
        db_site_info_lst = []
        db_info_lst = info.get("DB_COUNT", ["0"])
        if db_info_lst:
            db_count = int(db_info_lst[0])
            if db_count > 0:
                db_site_str_lst = info.get("DB_SITE", [])
                if db_site_str_lst:
                    db_site_str = db_site_str_lst[0]
                else:
                    db_site_str = ""
                db_site_info_dct = self.format_site_info(db_site_str)
                db_site_lst = natsorted(db_site_info_dct.get("site", []))
                db_site_info_lst = natsorted(db_site_info_dct.get("site_info", []))
                if db_site_info_lst:
                    db_lv = 0.2
                else:
                    if db_site_lst:
                        db_lv = 0.1
                    else:
                        db_lv = 0
        else:
            db_count = 0
        db_info = {
            "count": db_count,
            "cv": "",
            "level": db_lv,
            "order": 0.01,
            "site": db_site_lst,
            "site_info": db_site_info_lst,
        }

        return db_info

    @staticmethod
    def format_link(info: dict) -> str:
        link = ""
        link_lst = info.get("LINK", [""])
        if link_lst:
            link = link_lst[0]
            link = link.strip(" ")
            if re.match(r"^O[-]?$", link, re.IGNORECASE) or link == "e":
                link = "O-"
            elif re.match(r"^P[-]?$", link, re.IGNORECASE) or link == "p":
                link = "P-"
            else:
                pass
        else:
            pass

        return link

    def format_mod(self, info: dict) -> dict:
        formatted_mod_lst = []
        raw_mod_type_lst = info.get("MOD_TYPE", [])
        unique_raw_mod_type_lst = raw_mod_type_lst
        # unique_raw_mod_type_lst = list(set(raw_mod_type_lst))
        mod_type_lst = []
        if unique_raw_mod_type_lst and len(unique_raw_mod_type_lst) > 1:
            if unique_raw_mod_type_lst[0] not in ["DB", ""]:
                mod_type_lst.append(unique_raw_mod_type_lst[0])
                for idx in range(1, len(mod_type_lst) + 1):
                    if unique_raw_mod_type_lst[idx] not in ["DB", ""]:
                        mod_type_lst.append(unique_raw_mod_type_lst[idx])
            else:
                mod_type_lst = unique_raw_mod_type_lst
        else:
            mod_type_lst = unique_raw_mod_type_lst
        mod_count_lst = info.get("MOD_COUNT", [])
        mod_site_seg_sum_lst = info.get("MOD_SITE", [""])
        mod_site_lst = []
        mod_site_info_lst = []
        for mod_site_seg in mod_site_seg_sum_lst:
            mod_site_info_dct = self.format_site_info(mod_site_seg)
            mod_site_lst.append(mod_site_info_dct.get("site", []))
            mod_site_info_lst.append(mod_site_info_dct.get("site_info", []))
        if mod_type_lst == [""]:
            mod_count_lst = [""]
        formatted_mod_type_lst = []
        mod_lv_dct = {}
        mass_shift_dct = {}
        for mod_type in mod_type_lst:
            for cv_alias in self.alias2cv:
                alia_rgx = self.alias2cv[cv_alias].get("MATCH", None)
                matched_cv = self.alias2cv[cv_alias].get("CV", None)
                if alia_rgx and matched_cv is not None:
                    alia_match = alia_rgx.match(mod_type)
                    if alia_match:
                        self.logger.debug(
                            f"mod_type: {mod_type} identified as {matched_cv}"
                        )
                        # if matched_cv not in formatted_mod_type_lst:
                        #     formatted_mod_type_lst.append(matched_cv)
                        formatted_mod_type_lst.append(matched_cv)
                        if matched_cv == "Delta":
                            mod_lv_dct[matched_cv] = self.raw_cv["Delta"].get(
                                "LEVEL", 0
                            )
                            mass_shift_dct[matched_cv] = int(mod_type)
                            break
                        else:
                            mod_lv_dct[matched_cv] = max(
                                mod_lv_dct.get(matched_cv, 0),
                                self.alias2cv[cv_alias].get("LEVEL", 0),
                            )
        mod_type_count = len(formatted_mod_type_lst)
        if (
            mod_type_lst
            and formatted_mod_type_lst
            and len(mod_type_lst) == mod_type_count
        ):
            info["MOD_TYPE"] = formatted_mod_type_lst

        if len(mod_count_lst) == mod_type_count:
            pass
        elif len(mod_count_lst) > mod_type_count:
            mod_count_lst = mod_count_lst[:mod_type_count]
        elif len(mod_count_lst) < mod_type_count:
            mod_count_lst = mod_count_lst + [""] * (mod_type_count - len(mod_count_lst))

        if mod_type_count == len(mod_site_lst) == len(mod_site_info_lst):
            formatted_mod_lst = zip(
                mod_count_lst, formatted_mod_type_lst, mod_site_lst, mod_site_info_lst
            )
        elif mod_type_count == len(mod_site_info_lst):
            mod_site_lst = [""] * mod_type_count
            formatted_mod_lst = zip(
                mod_count_lst, formatted_mod_type_lst, mod_site_lst, mod_site_info_lst
            )
        elif mod_type_count == len(mod_site_lst):
            mod_site_info_lst = [""] * mod_type_count
            formatted_mod_lst = zip(
                mod_count_lst, formatted_mod_type_lst, mod_site_lst, mod_site_info_lst
            )
        elif mod_type_count > 0 and not mod_site_lst and not mod_site_info_lst:
            mod_site_lst = [""] * mod_type_count
            mod_site_info_lst = [""] * mod_type_count
            formatted_mod_lst = zip(
                mod_count_lst, formatted_mod_type_lst, mod_site_lst, mod_site_info_lst
            )
        else:
            if 0 < len(mod_site_lst) < mod_type_count:
                self.logger.debug(
                    f"mod_site_lst: {mod_site_lst} | formatted_mod_type_lst: {formatted_mod_type_lst}"
                )
                formatted_mod_lst = []
            # elif 0 < mod_type_count < len(mod_site_lst):
            #     if list(set(formatted_mod_type_lst)) == ["DB"]:
            #         formatted_mod_type_lst = ["DB"] * len(mod_site_lst)
            #         if not mod_site_info_lst:
            #             mod_site_info_lst = [""] * len(mod_site_lst)
            #         else:
            #             mod_site_info_lst = mod_site_info_lst + [""] * (
            #                 len(mod_site_lst) - len(mod_site_lst)
            #             )
            #         formatted_mod_lst = zip(
            #             mod_count_lst,
            #             formatted_mod_type_lst,
            #             mod_site_lst,
            #             mod_site_info_lst,
            #         )
            else:
                formatted_mod_lst = []

        mod_info_dct = {}
        if formatted_mod_lst:
            formatted_mod_lst = list(formatted_mod_lst)
            for mod_tp in formatted_mod_lst:
                delta_mod_count = mod_tp[0]
                if delta_mod_count and isinstance(delta_mod_count, str):
                    if delta_mod_count == "+":
                        delta_mod_count = 1
                    elif delta_mod_count == "-":
                        delta_mod_count = -1
                    else:
                        try:
                            delta_mod_count = int(delta_mod_count)
                        except (ValueError, TypeError):
                            delta_mod_count = 1
                else:
                    delta_mod_count = 1
                mod_type = mod_tp[1]
                mod_order = self.raw_cv[mod_type].get("ORDER", 0)
                existed_mod_count = mod_info_dct.get(f"{mod_order}_{mod_type}", {}).get(
                    "mod_count", 0
                )
                existed_mod_site_lst = mod_info_dct.get(
                    f"{mod_order}_{mod_type}", {}
                ).get("mod_site", [])
                existed_mod_site_info_lst = mod_info_dct.get(
                    f"{mod_order}_{mod_type}", {}
                ).get("mod_site_info", [])
                # e.g. mod_tp: ('', 'DB', '9', 'Z')
                existed_mod_site_lst.extend(mod_tp[2]),
                existed_mod_site_info_lst.extend(mod_tp[3])
                mod_level = 0
                mod_count = existed_mod_count + delta_mod_count

                mod_level = mod_lv_dct.get(mod_type, 0)
                if mod_tp[2] and mod_tp[3]:
                    mod_level += 2
                elif mod_tp[2] and not mod_tp[3]:
                    mod_level += 1
                else:
                    pass

                updated_mod_info = {
                    "count": mod_count,
                    "cv": mod_type,
                    "level": mod_level,
                    "order": mod_order,
                    "site": natsorted(existed_mod_site_lst),
                    "site_info": natsorted(existed_mod_site_info_lst),
                }
                verbose = {}
                if mod_type in self.raw_cv:
                    verbose["elements"] = self.raw_cv[mod_type].get("ELEMENTS", {})
                    if (
                        mod_type not in ["Delta", "DB"]
                        and mod_type not in mass_shift_dct
                    ):
                        verbose["mass_shift"] = self.to_mass_shift(verbose["elements"])
                    elif mod_type == "Delta" and mod_type in mass_shift_dct:
                        verbose["mass_shift"] = mass_shift_dct.get("Delta", 0)
                    else:
                        verbose["mass_shift"] = 0
                else:
                    raise ValueError(f"Unsupported modification type: {mod_type}")
                updated_mod_info["verbose"] = verbose
                mod_info_dct[f"{mod_order}_{mod_type}"] = updated_mod_info
        else:
            pass
        mod_seg_levels_lst = []
        if mod_info_dct:
            for mod_seg in mod_info_dct:
                mod_seg_levels_lst.append(mod_info_dct[mod_seg].get("level", 0))
        if mod_seg_levels_lst:
            max_mod_level = max(mod_seg_levels_lst)
        else:
            max_mod_level = 0
        # mod_info_dct["MOD_LEVEL"] = max_mod_level

        return {"level": max_mod_level, "info": mod_info_dct}

    def format_residue(self, info: dict) -> dict:

        link = self.format_link(info)
        c_count = int(info.get("C_COUNT", ["0"])[0])
        db_info = self.format_db(info)
        db_info_sum = {"level": db_info.get("level", 0), "info": {"0.01_DB": db_info}}
        sp_o_info = self.format_sp_o(info)
        sp_o_info_sum = {
            "level": sp_o_info.get("level", 0),
            "info": {"0.02_SP_O": sp_o_info},
        }
        mod_info_sum = self.format_mod(info)

        res_lv = mod_info_sum.get("level") + max(
            db_info_sum.get("level"), sp_o_info_sum.get("level")
        )

        residue_info_dct = {
            "link": link,
            "c_count": c_count,
            "db_count": db_info.get("count", 0),
            "db_info_sum": db_info_sum,
            "sp_o_count": sp_o_info.get("count", 0),
            "sp_o_info_sum": sp_o_info_sum,
            "mod_info_sum": mod_info_sum,
        }

        residue_dct = {"level": res_lv, "info": residue_info_dct}

        return residue_dct

    def format(self, info: dict) -> dict:
        formatted_info = {"residues": self.format_residue(info)}

        return formatted_info


if __name__ == "__main__":
    pass
