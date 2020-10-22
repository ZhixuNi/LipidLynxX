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

from io import BytesIO
import os
from pathlib import Path
from typing import Union

import pandas as pd
import typer


def cli_get_table(file: Union[Path, str]):
    table_header_lst = []
    if isinstance(file, Path) and file.is_file():
        in_file_name_low_str = file.name.lower()
    elif isinstance(file, str) and os.path.isfile(file):
        in_file_name_low_str = file.lower()
    else:
        typer.echo(f"[IO Error] Can not find file: {file}")
        raise typer.Exit(code=1)
    if in_file_name_low_str:
        if in_file_name_low_str.endswith("xlsx"):
            table_df = pd.read_excel(file)
            table_dct = table_df.to_dict(orient="list")
            table_header_lst = table_df.columns.values.tolist()
        elif in_file_name_low_str.endswith("csv"):
            try:
                table_df = pd.read_csv(file)
                table_dct = table_df.to_dict(orient="list")
                tab_table_df = pd.read_csv(file, sep="\t")
                tab_table_dct = tab_table_df.to_dict(orient="list")
                if len(list(tab_table_dct.keys())) > len(list(table_dct)):
                    print(f"{in_file_name_low_str} is identified as Tab separated csv.")
                    table_dct = tab_table_dct
                    table_df = tab_table_df
                else:
                    pass
            except pd.errors.ParserError:
                print(f"{in_file_name_low_str} is Tab separated csv.")
                table_df = pd.read_csv(file, sep="\t")
                table_dct = table_df.to_dict(orient="list")
            table_header_lst = table_df.columns.values.tolist()

        elif in_file_name_low_str.endswith("tsv"):
            table_df = pd.read_csv(file, sep="\t")
            table_dct = table_df.to_dict(orient="list")
            table_header_lst = table_df.columns.values.tolist()
        # elif in_file_name_low_str.endswith("txt"):
        #     with open(file) as f_obj:
        #         table_dct = {"input": f_obj.readlines()}
        else:
            typer.echo("File type not supported")
            raise typer.Exit(code=1)
    else:
        typer.echo(f"[IO Error] Can not find file.")
        raise typer.Exit(code=1)

    return table_dct, table_header_lst


def cli_save_output(output_info: Union[str, BytesIO], output_file: Path):
    if isinstance(output_info, str):
        typer.echo(
            typer.style(
                f"Save output as: {output_file.as_posix()}", fg=typer.colors.CYAN
            )
        )
        raise typer.Exit(code=0)
    else:
        typer.echo(f"Failed to generate output: {output_file.as_posix()}")
        raise typer.Exit(code=1)


if __name__ == '__main__':

    f_csv = r"../../doc/sample_data/input/LipidLynxX_test.csv"
    f_t_csv = r"../../doc/sample_data/input/LipidLynxX_test_tab.csv"
    f_tsv = r"../../doc/sample_data/input/LipidLynxX_test.tsv"
    f_xlsx = r"../../doc/sample_data/input/LipidLynxX_test.xlsx"
    for f in [f_csv, f_t_csv, f_tsv, f_xlsx]:
        t = cli_get_table(f)
        print(t)

    t2 = cli_get_table(f_t_csv)
    print(t2)
