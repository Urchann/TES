#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
 _    _      _        _                                                        _                  
| |_ | |__  (_)  ___ | | __ _ __    ___  ___  ___          ___  _ __ ___    __| |    _ __   _   _ 
| __|| '_ \ | | / __|| |/ /| '_ \  / _ \/ __|/ __|        / __|| '_ ` _ \  / _` |   | '_ \ | | | |
| |_ | | | || || (__ |   < | | | ||  __/\__ \\__ \       | (__ | | | | | || (_| | _ | |_) || |_| |
 \__||_| |_||_| \___||_|\_\|_| |_| \___||___/|___/ _____  \___||_| |_| |_| \__,_|(_)| .__/  \__, |
                                                  |_____|                           |_|     |___/ 
"""
# ==============================================================================================
#  File:        tes/scripts/thickness_cmd.py
#  Project:     Transition Edge Sensor (TES) Microcalorimeter
#  Module:      Thickness Analysis Command Execution
#
#  Description:
#      This module is intended to drive `tes/libs/core/thickness.py`.
#      It automatically corrects the naming of analysis files
#      according to predefined rules during processing.
#
#  Requirements:
#      - Python 3.10+
#      - libs.core.thickness (Internal module)
#      - libs.utils (Internal module)
#
#  Author:      Ryota Fukuda (mailto:25la018c@rikkyo.ac.jp)
#  Maintainer:  Yamada Laboratory, Department of Physics, College of Science, Rikkyo University
#  Url:         https://github.com/Urchann/TES/tree/main/tes/libs/utils.py
#  Created:     2025-10-09
#  Updated:     2025-12-31
#  Version:     1.1.2
#
#  ChangeLog:
#      v1.0.0 (2025-10-09)
#         - Build the execution system for `tes/libs/core/thickness.py`.
#      v1.1.0 (2025-10-21)
#         - Adapted to `tes/libs/utils.py`.
#      v1.1.1 (2025-11-20)
#         - Added exclude_files argument to skip specific files from analysis
#           for fix_thickness_filename function.
#      v1.1.2 (2025-11-20)
#         - Added floor and mode argument to instances of the Thickness class
#           for analyzing step-function shaped roughness data.
#
#  Notes:
#      - For detailed information, please see thickness.py and utils.py files.
#
# ==============================================================================================


import os
import sys
import time
import itertools
import threading
import contextlib
import io
from libs.core.thickness import Thickness
from libs.utils import fix_thickness_filename


_stop_animation = threading.Event()


def _animate_status(method):
    """
    バックグラウンドでアニメーションステータスをコンソールに表示し続ける関数。
    三点リーダーのアニメーションを付与し、ユーザにプログラムが動作中であることを伝える。

    Parameters:
        method (str): 実行中のメソッド名。
    
    Returns:
        None        : 返り値なし。標準出力にアニメーションを出力する。
    """
    chars   = ['   ', '.  ', '.. ', '...'] 
    spinner = itertools.cycle(chars)

    real_stdout = sys.__stdout__
    while not _stop_animation.is_set():
        status_text = f" Running {method} {next(spinner)}"
        real_stdout.write(f"\r{status_text}   ")
        real_stdout.flush()
        time.sleep(0.3) 
    
    real_stdout.write("\r" + " " * 50 + "\r")
    real_stdout.flush()

def _run_method(method, target_func=None, *args):
    """
    時間が可変の解析処理をシミュレートし、アニメーションステータスを表示する関数。

    Parameters:
        method (str)      : 実行中のメソッド名。
        target_func (func): 実行したい実際の関数。
        *args             : target_funcに渡す引数。
    
    Returns:
        None          : 返り値なし。
    """
    global _stop_animation
    _stop_animation.clear()
    buf = io.StringIO()
    
    animation_thread = threading.Thread(target=_animate_status, args=(method,), daemon=True)
    animation_thread.start()

    result = None
    try:
        if target_func:
            with contextlib.redirect_stdout(buf):
                result = target_func(*args)
        else:
            result = None
    except Exception as e:
        _stop_animation.set()
        animation_thread.join()
        raise e
    
    _stop_animation.set()
    animation_thread.join()
    real_stdout = sys.__stdout__
    sys.stdout.write(f"\r Running {method} --> \033[1mCompleted!\033[0m\n")
    real_stdout.flush()
    
    logs = buf.getvalue()
    if logs:
        real_stdout.write(logs + "\n")
        real_stdout.flush()

    return result


def run_analysis(thickness, board_img, pos_group):
    """
    解析結果のプロットおよびファイル出力処理を順序立てて実行するパイプライン関数。
    各処理の実行をシミュレートし、アニメーション表示を制御する。

    Parameters:
        thickness (Thickness): 厚み解析メソッドを持つThicknessクラスのインスタンス。
        board_img (str)      : 厚み解析の結果で、相関を見る測定位置の概略図のパス。
        pos_group (list)     : 厚み解析の結果で、相関を見る測定位置のリスト。
    
    Returns:
        None                 : 返り値なし。
    """
    pos = [item.split('-') for item in pos_group]
    def _plot_all_results():
        thickness.plot_roughtness()
        thickness.plot_test_correlation()
        for i in range(len(pos)):
            thickness.plot_pos_correlation(pos[i])
    _run_method("Plot Results", _plot_all_results)
    
    _run_method("Output Results", thickness.export_results, board_img, *pos_group)


def run_thickness(args, board_img, pos_group):
    """
    厚み解析のメイン処理を実行する制御関数。
    ファイル名のチェック、Thicknessクラスの初期化、解析の計算、および結果の出力パイプラインの実行を順序立てて管理する。

    Parameters:
        args (Namespace): コマンドライン引数（入力/出力ディレクトリ、設定値など）を含むオブジェクト。
        board_img (str) : 厚み解析の結果で、相関を見る測定位置の概略図のパス。
        pos_group (list): 厚み解析の結果で、相関を見る測定位置のリスト。
    
    Returns:
        None: 返り値なし。処理の成否は標準出力にメッセージとして出力される。
    """
    if fix_thickness_filename(args.dir_input, args.exclude_files, os.getcwd()):
        print(" Please rerun after correcting filenames.")
        sys.exit(0)

    try:
        thickness = Thickness(
            dir_input     = args.dir_input,
            dir_output    = args.output_dir,
            dat_exclude   = args.exclude_files,
            dat_column    = args.data_column,
            output_format = args.output_format,
            smooth_window = args.window_size,
            margin        = args.margin,
            boundary      = args.boundary,
            threshold     = args.threshold,
            mode          = args.mode
        )
        
        thickness.parse_filename()
        thickness.open_file()
        thickness.calculate_thickness()

        os.makedirs(thickness.dir_output, exist_ok=True)
        print(" --- Summarizing Results ---")
        print(f" Output directory created/verified: {thickness.dir_output}")
        
        run_analysis(thickness, board_img, pos_group)
        print("\033[32m Analysis Success! All plots and results generated.\033[0m")

    except FileNotFoundError as e:
        print(f'\n \033[91mFATAL ERROR\033[0m: {e}')
        print(" Please ensure the input directory is correct relative to the 'data' folder.")
    except Exception as e:
        print(f"\n \033[91mAn unexpected error occurred: {e}\033[0m")
