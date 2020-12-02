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

import regex as re


def import_biopan(lipid):

    biopan_name = "UNPROCESSED"
    lipid_class = lipid.lipid_class
    lv = lipid.lv
    residues = lipid.residues
    is_modified = lipid.is_modified

    if re.match(r"[MS].*", lv.upper()):
        if lipid_class in ["Cer", "SM"]:
            is_spb = False
            if len(residues) == 1:
                residue = residues[0]
                num_c = residue.num_c
                num_db = residue.num_db
                num_o = residue.num_o

                if num_c == 18 and num_db == 1 and num_o == 2:
                    is_spb = True
                    return biopan_name
                else:
                    lipid.residues.append(residue("SPB"))
    else:
        if lipid_class in ["Cer", "SM"]:
            has_spb = False
            for residue in residues:
                num_c = residue.num_c
                num_db = residue.num_db
                num_o = residue.num_o
                if num_c < 30 and num_o < 2:
                    lipid.residues.append(residue("SPB"))

    return lipid


def export_biopan(lipid):

    biopan_name = "UNPROCESSED"
    lipid_class = lipid.lipid_class
    lv = lipid.lv
    residues = lipid.residues
    is_modified = lipid.is_modified

    if is_modified:
        pass
    else:
        if re.match(r"[MS].*", lv.upper()):
            if lipid_class in ["Cer", "SM"]:
                has_spb = False
                fa_obj_lst = []
                for residue in residues:
                    num_c = residue.num_c
                    num_db = residue.num_db
                    num_o = residue.num_o

                    if num_c == 18 and num_db == 1 and num_o == 2:
                        has_spb = True
                    else:
                        fa_obj_lst.append(residue)
                if has_spb and len(fa_obj_lst) == 1:
                    biopan_name = lipid
        else:
            if lipid_class in ["Cer", "SM"]:
                has_spb = False
                for residue in residues:
                    num_c = residue.num_c
                    num_db = residue.num_db
                    num_o = residue.num_o
                    if num_c < 30 and num_o < 2:
                        biopan_name = lipid

    return biopan_name
