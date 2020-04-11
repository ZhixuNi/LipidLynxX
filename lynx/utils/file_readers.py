# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2020  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import json
import os
from io import BytesIO
from typing import List, Union

import pandas as pd
from werkzeug.datastructures import FileStorage
from werkzeug.utils import secure_filename

from lynx.utils.log import logger


def get_abs_path(file_path: str) -> str:
    """
    check and get absolute file path from given file
    Args:
        file_path: The relative path to a file

    Returns:
        abs_path: the absolute path of input file

    """

    abs_path = ""

    if os.path.isdir(file_path):
        abs_path = os.path.abspath(file_path)

    elif os.path.isfile(file_path):
        abs_path = os.path.abspath(file_path)
    else:
        in_file_lst = [
            r"{path}".format(path=file_path),
            r"./{path}".format(path=file_path),
            r"../{path}".format(path=file_path),
            r"../../{path}".format(path=file_path),
            r"../../../{path}".format(path=file_path),
        ]
        for f in in_file_lst:
            if os.path.isdir(f):
                abs_path = os.path.abspath(f)
                break
            elif os.path.isfile(f):
                abs_path = os.path.abspath(f)
                break

    if not abs_path:
        raise FileNotFoundError(f"Can not find file: {file_path}")

    return abs_path


def get_table(file: Union[str, FileStorage]) -> dict:
    if isinstance(file, str):
        try:
            abs_path = get_abs_path(file)
            file = abs_path
        except FileNotFoundError:
            raise FileNotFoundError
    elif isinstance(file, FileStorage):
        abs_path = secure_filename(file.filename)
    else:
        raise FileNotFoundError

    if abs_path.lower().endswith(".csv"):
        df = pd.read_csv(file)
    elif abs_path.lower().endswith(".tsv"):
        df = pd.read_csv(file, sep="\t")
    elif abs_path.lower().endswith(".xlsx") or abs_path.lower().endswith(".xls"):
        df = pd.read_excel(file)
    else:
        df = pd.DataFrame()

    if not df.empty:
        dct = df.to_dict(orient="list")
    else:
        dct = {}

    return dct


def get_json(file: str) -> dict:
    file = get_abs_path(file)
    if file.lower().endswith(".json"):
        with open(file) as file_obj:
            js_obj = json.load(file_obj)
            return js_obj
    else:
        raise IOError(f"Input file: {file} is not json file")


def load_folder(folder: str, file_type: str = "") -> List[str]:
    """
     Load all files under given folder, optional with selected file suffix
     Args:
         folder: path of the folder.
         file_type: type of the file, default value is "" for no file type filter

     Returns:
         file_abs_path_lst: the list of files under given folder in absolute path

     """
    abs_path = get_abs_path(folder)
    file_lst = os.listdir(abs_path)
    file_abs_path_lst = [os.path.join(abs_path, x) for x in file_lst]
    if file_type:
        file_abs_path_lst = [
            f for f in file_abs_path_lst if f.lower().endswith(file_type.lower())
        ]
    file_abs_path_lst = [abs_f for abs_f in file_abs_path_lst if os.path.isfile(abs_f)]
    logger.debug(
        f"Fund {file_type} files:\nunder folder: {folder}\nfiles:\n {file_abs_path_lst}"
    )

    return file_abs_path_lst


def create_output(data: dict) -> BytesIO:
    excel_io = None
    converted_df = pd.DataFrame()
    not_converted_df = pd.DataFrame()
    if data:
        not_converted_dct = {}
        df_lst = []
        for k in data:
            if isinstance(data[k], dict):
                k_pairs = data[k].get("converted", [])
                k_not_converted = data[k].get("skipped", [])
                if k_pairs and isinstance(k, str):
                    df_lst.append(pd.DataFrame(k_pairs, columns=[k, f"{k}_converted"]))

                if k_not_converted:
                    not_converted_dct[f"{k}_skipped"] = k_not_converted
            elif isinstance(data[k], list) and k == "converted":
                k_pairs = data.get("converted", [])
                if k_pairs:
                    df_lst.append(
                        pd.DataFrame(k_pairs, columns=["input", f"converted"])
                    )
            elif isinstance(data[k], list) and k == "skipped":
                k_not_converted = data.get("skipped", [])
                if k_not_converted:
                    not_converted_dct[f"skipped"] = k_not_converted

        if df_lst:
            converted_df = pd.concat(df_lst, axis=1)

        if not_converted_dct:
            not_converted_df = pd.DataFrame.from_dict(
                not_converted_dct, orient="index"
            ).T

        if not converted_df.empty:
            excel_io = BytesIO()
            output_writer = pd.ExcelWriter(
                excel_io, engine="openpyxl"
            )  # write to BytesIO instead of file path
            converted_df.to_excel(output_writer, sheet_name="converted", index=False)
            if not not_converted_df.empty:
                not_converted_df.to_excel(
                    output_writer, sheet_name="skipped", index=False
                )
            output_writer.save()
            excel_io.seek(0)

    return excel_io


def create_equalizer_output(sum_data: dict) -> BytesIO:

    if sum_data:
        table_io = BytesIO()
        table_writer = pd.ExcelWriter(
            table_io, engine="openpyxl"
        )  # write to BytesIO instead of file path
        for lv in sum_data:
            data = sum_data[lv]
            for k in data:
                if k.lower().startswith("match"):
                    matched_dct = data[k]
                    if matched_dct:
                        out_matched_df = pd.DataFrame.from_dict(
                            matched_dct, orient="index"
                        )
                        out_matched_df.index.names = [f"ID@Lv_{lv}"]
                        out_matched_df.sort_index().to_excel(
                            table_writer, sheet_name=f"matched_{lv}"
                        )
                elif k.lower().startswith("equalized"):
                    equalized_dct = data[k]
                    if equalized_dct:
                        pd.DataFrame.from_dict(
                            equalized_dct, orient="index"
                        ).sort_index().to_excel(table_writer, sheet_name=f"unmatched")
                elif k.lower().startswith("skipped"):
                    skipped_dct = data[k]
                    if skipped_dct:
                        pd.DataFrame.from_dict(
                            skipped_dct, orient="index"
                        ).T.sort_index().to_excel(
                            table_writer, sheet_name="skipped", index=False
                        )
                else:
                    pass

        table_writer.save()
        table_io.seek(0)
    else:
        table_io = None
    return table_io


def save_table(df: pd.DataFrame, file_name: str) -> (bool, str):
    is_output = False
    abs_output_path = None
    if not df.empty:
        df.to_excel(file_name)
        is_output = True
        abs_output_path = get_abs_path(file_name)

    return is_output, abs_output_path
