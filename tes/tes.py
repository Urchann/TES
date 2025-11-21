#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
 _                               
| |_   ___  ___     _ __   _   _ 
| __| / _ \/ __|   | '_ \ | | | |
| |_ |  __/\__ \ _ | |_) || |_| |
 \__| \___||___/(_)| .__/  \__, |
                   |_|     |___/ 
"""
# ==============================================================================================
#  File:        tes/tes.py
#  Project:     Transition Edge Sensor (TES) Microcalorimeter
#  Module:      Main CLI for TES Analysis
#
#  Description:
#      This script serves as the central entry point for the TES analysis framework.
#      It utilizes `argparse` to parse command-line arguments and dispatches 
#      execution to specific analysis modules.
#
#  Requirements:
#      - Python 3.10+
#      - scripts.thickness_cmd (Internal module)
#
#  Author:      Ryota Fukuda (mailto:25la018c@rikkyo.ac.jp)
#  Maintainer:  Yamada Laboratory, Department of Physics, College of Science, Rikkyo University
#  Url:         https://github.com/Urchann/TES/tree/main/tes/tes.py
#  Created:     2025-10-09
#  Updated:     2025-12-31
#  Version:     1.1.0
#
#  ChangeLog:
#      v1.0.0 (2025-10-09)
#         - Built the main CLI framework using `argparse`.
#           Implemented the routing system to `scripts.thickness_cmd`.
#      V1.0.1 (2025-11-18)
#         - Implemented board configuration that allows switching options 
#           via comment toggles depending on the measurement target.
#      v1.1.0 (2025-11-20)
#         - Added `--exclude_files` argument to the thickness subcommand
#           to allow skipping specific TXT files during analysis.
#         - Added `--mode` argument to support 'step' (staircase) waveform analysis.
#         - Added `--margin` and `--boundary` for threshold detection.
#
#  Notes:
#      - Execute this script from the project root directory.
#
# ==============================================================================================


import os
import sys
import argparse
from scripts.thickness_cmd import run_thickness


def main():
    """
    解析ツールのメインエントリポイント。
    コマンドライン引数をパースし、設定されたパラメータに基づいて解析を実行する。

    Parameters:
        dir_input (str): 入力TXTファイルを含むディレクトリ (必須)。
        --exclude_files (list): 解析から除外するファイル名のリスト。
        --output_dir (str): 結果の保存先ディレクトリ。
        --data_column (str): 使用するデータ列 (Normal, Raw, etc.)。
        --window_size (int): 平滑化のための移動平均ウィンドウサイズ。
        --margin
        --boundary
        --threshold (float): 境界判定や異常値除去のための閾値パラメータ。
        --mode (str): 解析モード ('square' or 'step')。

    Returns:
    -------
    None
        戻り値はありませんが、実行結果として以下の処理が行われます:
        - 指定された出力ディレクトリへのグラフ・レポートの保存。
        - 処理経過の標準出力への表示。
    
    Notes
    -----
    解析対象の基板設定（pos_group, board_img）は現在この関数内にハードコードされています。
    対象を変更する場合は、コード内の該当箇所のコメントアウトを切り替えてください。
    """
    parser = argparse.ArgumentParser(
        description="Analysis framework",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_thickness = subparsers.add_parser("thickness")
    parser_thickness.add_argument("dir_input",
                                  help="Input directory containing the TXT files (e.g., 'test').\n" \
                                       "This will be automatically resolved relative to the 'data/' directory.")
    parser_thickness.add_argument("--exclude_files",
                                  nargs='*',
                                  default=None,
                                  help="Space-separated list of TXT files to exclude from analysis (e.g., 'file1.txt file2.txt').")
    parser_thickness.add_argument("--output_dir",
                                  default=None,
                                  help="Directory to save plots and results table. (Default: 'input_dir_name'_results)")
    parser_thickness.add_argument("--data_column",
                                  choices=['Raw', 'RawLevel', 'Normal', 'Rough', 'Wavi'],
                                  default='Normal',
                                  help="Column name to use for measurement data. (Default: 'Normal')")
    parser_thickness.add_argument("--output_format",
                                  choices=['pdf', 'html'],
                                  default=None,
                                  help="Output format for the results table: 'pdf' (LaTeX) or 'html'. (Default: None)")
    parser_thickness.add_argument("--window_size",
                                  type=int,
                                  default=50,
                                  help="Window size for data smoothing (moving average). (Default: 50)")
    parser_thickness.add_argument("--margin",
                                  type=float,
                                  default=50,
                                  help="Safety buffer margin to ensure measurement accuracy. (Default: 50)")
    parser_thickness.add_argument("--boundary",
                                  type=float,
                                  default=0.1,
                                  help="Threshold for determining rising and falling boundaries. (Default: 0.1)")
    parser_thickness.add_argument("--threshold",
                                  type=float,
                                  default=0.5,
                                  help="Relative threshold for extreme values. (Default: 0.5)")
    parser_thickness.add_argument("--mode",
                                  choices=['square', 'step'],
                                  default='square',
                                  help="Roughness data waveform. (Default: 'square')")

    args = parser.parse_args()

    if args.command == "thickness":
        # ---------------------------------------------------------
        # 基板設定（測定対象に合わせてコメントアウトを切り替える）
        # ---------------------------------------------------------
        # Config 1: LOR3A Board
        # pos_group = ['Ea1-3-5', 'Ea2-3-4']
        # board_img = './data/roughness/LOR3A_BoardSchematic.png'

        # Config 2: AZECI3012 Board
        pos_group = ['W01-W08-W15', 'W07-W08-W09']
        board_img = './data/roughness/AZECI3012_BoardSchematic.png'

        # Config 3: RD250617 Board
        # pos_group = ['NL-NR-CL-CR-SL-SR', 'WL-WR-CL-CR-EL-ER']
        # board_img = './data/roughness/RD250617_BoardSchematic.png'

        run_thickness(args, board_img, pos_group)
    else:
        parser.print_help()


if __name__ == "__main__":
    script_dir  = os.path.dirname(os.path.abspath(__file__))
    current_dir = os.getcwd()

    if script_dir != current_dir:
        print(f"\n\033[91mFATAL ERROR\033[0m: Execution location is incorrect.")
        print(f"Please run the script from the project root directory: '{script_dir}'")
        print(f"Current directory: '{current_dir}'")
        sys.exit(1)
    
    main()
