# ユーザーが実行するメインの解析スクリプト


#!/bin/bash

# run.py (または単なる実行指示メモ)
#
# 解析プログラム (tes.py) を実行する例

# 1. 最小構成の実行例: data/test ディレクトリを解析し、結果をPDFで出力 (デフォルト)
echo "--- Running Basic Analysis ---"
python3 ../tes.py test --data_column Normal

# 2. 詳細設定の実行例: data/LOR3A_TEST ディレクトリを解析し、HTMLで出力
echo ""
echo "--- Running Detailed Analysis ---"
python3 ../tes.py LOR3A_TEST