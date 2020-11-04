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

import pandas as pd

# Get rhea-directions.tsv
# https://www.rhea-db.org/help/reaction-side-direction
# ftp://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv
rhea_tsv = r'../../test/test_input/rhea-directions.tsv'
rhea_json = r'../../test/test_output/rhea-directions.json'
rhea_col = ["RHEA_ID_MASTER", "RHEA_ID_LR", "RHEA_ID_RL", "RHEA_ID_BI"]

rhea_df = pd.read_csv(rhea_tsv, sep="\t")

raw_rhea_dct = {}
rhea_dct = {}

for col in rhea_col:
    tmp_df = pd.DataFrame(data=rhea_df)
    tmp_df[f"{col}_index"] = tmp_df[col]
    tmp_df.set_index(f"{col}_index", drop=True, inplace=True)
    print(col)
    print(tmp_df.head())
    raw_rhea_dct.update(tmp_df.to_dict(orient="index"))
    del tmp_df

rhea_ids = sorted(list(set(list(raw_rhea_dct.keys()))))

for rhea_id in rhea_ids:
    rhea_dct[rhea_id] = raw_rhea_dct.get(rhea_id)

with open(rhea_json, "w+") as js_obj:
    json.dump(rhea_dct, js_obj)
print(f"Rhea JSON saved as: {rhea_json}")
