#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
        _    _  _                       
 _   _ | |_ (_)| | ___     _ __   _   _ 
| | | || __|| || |/ __|   | '_ \ | | | |
| |_| || |_ | || |\__ \ _ | |_) || |_| |
 \__,_| \__||_||_||___/(_)| .__/  \__, |
                          |_|     |___/ 
"""
# ==============================================================================================
#  File:        tes/libs/utils.py
#  Project:     Transition Edge Sensor (TES) Microcalorimeter
#  Module:      Utility functions for file naming correction
#
#  Description:
#      This module checks and corrects roughness measurement filenames.
#      It ensures wafer name, date, position, direction, and numbering
#      follow the project naming rules.
#
#  Requirements:
#      - Python 3.10+
#      - No external dependencies (Standard library only)
#
#  Author:      Ryota Fukuda (mailto:25la018c@rikkyo.ac.jp)
#  Maintainer:  Yamada Laboratory, Department of Physics, College of Science, Rikkyo University
#  Url:         https://github.com/Urchann/TES/tree/main/tes/libs/utils.py
#  Created:     2025-10-21
#  Updated:     2025-12-31
#  Version:     1.0.2
#
#  ChangeLog:
#      v1.0.0 (2025-10-21)
#         - Establishment of automatic file names correction logic
#           for thickness analysis data files.
#      v1.0.1 (2025-11-18)
#         - Fixed to prevent creation of input directory if it does not exist,
#           and cause the process to exit with an error.
#      v1.0.2 (2025-11-20)
#         - Added exclude_files argument to skip specific files from analysis
#           for fix_thickness_filename function.
#
#  Notes:
#      - To avoid false positives in automatic corrections,
#        be sure to update the logic whenever you change the processing rules.
#      - As a renaming program, it is not yet versatile enough to handle file names.
#
# ==============================================================================================


import os
import glob
import re
import sys
import shutil
from datetime import datetime


RE_DATE   = re.compile(r'^\d{6}$|^\d{8}$')
RE_TEST   = re.compile(r'^test(\d+)$')
RE_NUMBER = re.compile(r'^\d+$')


def _correct_date(part):
    """
    ファイル名の要素から測定年月日を抽出し、8桁に変換して返す関数。
    測定日要素の形式は6桁または8桁が有効 (例: 2025年12月1日 --> 251201 or 20251201)。
    6桁の日付に対しては、前に '20' を付加して8桁に変換する。

    Parameters:
        part (str)     : ファイル名の要素。
    
    Returns:
        - (str or None): 8桁に変換された測定日要素。
                         測定日要素が6桁または8桁以外の場合は None を返す。
    """
    if not RE_DATE.match(part):
        return None
    if len(part) == 8:
        return part
    return '20' + part


def _correct_numbering(part_key, exist_numbers):
    """
    指定されたファイル名の要素に対して、重複しない一意の連番を見つけて返す関数。
    使用済みの連番を更新し、次回以降の重複チェックに利用する。

    Parameters:
        part_key (str)      : 連番を割り当てる基準となるファイル名の要素。
        exist_numbers (dict): ファイル名の要素に対する使用済みの連番の集合。
    
    Returns:
        unique_number (int) : 重複していない最小の連番。
    """
    used_numbers = exist_numbers.get(part_key, set())
    unique_number = 1
    while unique_number in used_numbers:
        unique_number += 1
    used_numbers.add(unique_number)
    exist_numbers[part_key] = used_numbers
    return unique_number


def fix_thickness_filename(input_dir, exclude_files, root_dir):
    """
    tes/libs/core/thickness_measurement.py の解析対象のファイル名を校正する関数。
    命名規則は <基板名>_<test番号>_<測定日>_<その他基板情報>_<測定方向>_<測定位置>_<番号>.txt で、
    各要素に対して以下の処理を行う:
    - <基板名>
        > 抽出文字列: 'ファイルが格納されているディレクトリ名 (アンダースコアはハイフンに変更)'
        > 抽出失敗時: 上記を取得して <基板名> とする。
    - <test番号>
        > 抽出文字列: 'test(数字)'
        > 抽出失敗時: test1 から連番でファイル名が重複しないように <test番号> とする。
    - <測定日>
        > 抽出文字列: '6桁または8桁の数字 (数字は8桁に修正して統一)'
        > 抽出失敗時: プログラム実行時の年月日を8桁で付ける。
    - <測定方向>
        > 抽出文字列: 'LtoR または RtoL'
        > 抽出失敗時: エラーを出力してプログラムを強制終了する。
    - <測定位置>
        > 抽出文字列: *抽出なし*
        > 抽出失敗時: 基板名、test番号、測定日、測定方向、番号を除く要素のうち、最後の文字列を <測定位置> とする。
    - <番号>
        > 抽出文字列: 'ファイル名の末端要素の数字'
        > 抽出失敗時: 1 から連番でファイル名が重複しないように <番号> とする。
    <基板名> と <その他基板情報> については、アンダースコアが含まれることを想定して、アンダースコアはハイフンに改める。

    Parameters:
        input_dir (str)     : 全解析ファイルが格納されているディレクトリ名。
        root_dir (str)      : プログラム実行時のカレントディレクトリ名。
        exclude_files (list): 解析から除外するファイル名（ベースネーム）のリスト。
    
    Returns:
        -                   : 全解析ファイル名を命名ルールに基づき改名更新。
    """
    input_path = os.path.join(root_dir, 'data', 'roughness', input_dir)
    if not os.path.isdir(input_path):
        return False
    
    print("\n --- File Naming Check and Correction ---")
    staging_dir = os.path.join(input_path, '_staging_rename_')

    if exclude_files is None:
        exclude_files = []
    all_files = sorted(glob.glob(os.path.join(input_path, '*.txt')))
    file_list = []
    for file in all_files:
        basename = os.path.basename(file)
        if basename not in exclude_files:
            file_list.append(file)

    wafer_name    = input_dir.replace('_', '-')
    dir_name      = input_dir.split('_')
    current_date  = datetime.now().strftime('%Y%m%d')
    exist_numbers = {}
    rename_list   = []

    if os.path.exists(staging_dir):
        print(f"    \033[93mWARNING:\033[0m Found existing staging directory. Removing: {os.path.basename(staging_dir)}")
        shutil.rmtree(staging_dir)
    os.makedirs(staging_dir)
    print(f"       \033[94mINFO:\033[0m Created staging directory: {os.path.basename(staging_dir)}")

    for file in file_list:
        basename = os.path.basename(file)
        filename = os.path.splitext(basename)[0]
        parts    = filename.split('_')

        part_wafer_name  = None
        part_test_number = None
        part_date        = None
        part_direction   = None
        part_position    = None
        part_numbering   = None
        part_other_info  = []

        for idx, p in enumerate(parts):
            if p.replace('_', '-') == wafer_name:
                part_wafer_name = wafer_name
                continue

            test_number = RE_TEST.match(p)
            if test_number:
                part_test_number = f"test{test_number.group(1)}"
                continue

            date = _correct_date(p)
            if date:
                part_date = date
                continue

            if p in ['LtoR', 'RtoL']:
                part_direction = p
                continue

            if idx == len(parts) - 1 and RE_NUMBER.match(p):
                part_numbering = int(p)
                continue

            if p in dir_name and idx != 0:
                continue

            wafer_name_element = wafer_name.split('-')
            wafer_name_duplicate = False
            for element in wafer_name_element:
                if element == p:
                    wafer_name_duplicate = True
                    break
            
            if wafer_name_duplicate:
                continue

            part_other_info.append(p)

        if part_wafer_name is None:
            print(f"    \033[93mWARNING:\033[0m Wafer name in '{filename}' is incorrect.")
            part_wafer_name = wafer_name

        if part_test_number is None:
            print(f"    \033[93mWARNING:\033[0m Test number in '{filename}' is missing.")
            part_test_number = "test1"

        if part_date is None:
            print(f"    \033[93mWARNING:\033[0m Date information in '{filename}' is missing.")
            part_date = current_date

        if part_direction is None:
            print(f"      \033[91mERROR:\033[0m Measurement direction (LtoR or RtoL) in '{basename}' is missing.")
            sys.exit(0)

        if part_numbering is None:
            print(f"    \033[93mWARNING:\033[0m File numberig in '{filename}' is missing.")
            key = f"{part_wafer_name}_{part_test_number}_{part_date}_{'-'.join(part_other_info)}_{part_direction}_{part_position}"
            part_numbering = _correct_numbering(key, exist_numbers)

        if len(part_other_info) > 0:
            part_position   = part_other_info[-1]
            part_other_info = part_other_info[:-1]
        else:
            print(f"    \033[91mERROR:\033[0m Measurement position in '{basename}' is missing.")
            continue

        re_part_other_info = '-'.join(p.replace('_', '-') for p in part_other_info)
        if re_part_other_info:
            new_filename = f"{part_wafer_name}_{part_test_number}_{part_date}_{re_part_other_info}_{part_direction}_{part_position}_{part_numbering}.txt"
        else:
            new_filename = f"{part_wafer_name}_{part_test_number}_{part_date}_{part_direction}_{part_position}_{part_numbering}.txt"
        new_filepath = os.path.join(staging_dir, new_filename)

        if new_filename != basename:
            print(f"     \033[96mRENAME:\033[0m '{basename}' -> '{new_filename}'")
            try:
                shutil.copy2(file, new_filepath) 
                
                rename_list.append({
                    'old_path': file,
                    'old_name': basename,
                    'new_path': new_filepath,
                    'new_name': new_filename
                })
                
            except Exception as e:
                print(f"    \033[91mFATAL ERROR:\033[0m Failed to copy file to staging: {e}")
                shutil.rmtree(staging_dir)
                sys.exit(1)
        else:
            shutil.copy2(file, os.path.join(staging_dir, basename))


    rename_count = len(rename_list)
    
    if rename_count > 0:
        print('\n --- ATTENTION: File Naming Correction Proposal ---')
        print(f' {rename_count} filename modifications have been suggested.')
        
        user_input = input(' Apply all the proposed changes to the original directory? (y/n): ').strip().lower()
        if user_input == 'y':
            print('\n Applying autofix to the original directory...')
            
            try:
                for item in rename_list:
                    if os.path.exists(item['old_path']):
                        os.remove(item['old_path'])
                    final_path = os.path.join(input_path, item['new_name'])
                    shutil.move(item['new_path'], final_path)
                    
                for staged_file in os.listdir(staging_dir):
                    source = os.path.join(staging_dir, staged_file)
                    destination = os.path.join(input_path, staged_file)
                    if not any(item['new_name'] == staged_file for item in rename_list):
                        shutil.move(source, destination)

                shutil.rmtree(staging_dir)
                print(f' \033[92mSUCCESS:\033[0m {rename_count} file(s) automatically corrected and applied.')
                print(f"       \033[94mINFO:\033[0m Staging directory deleted.")
                return True
                
            except Exception as e:
                print(f"  \033[91mFATAL ERROR:\033[0m An error occurred during application. Cleanup failed: {e}")
                return False
        
        else:
            print(" \033[91mCANCELED:\033[0m Changes rejected by user. Original files remain unchanged.")
            shutil.rmtree(staging_dir)
            print(f"       \033[94mINFO:\033[0m Staging directory deleted.")
            return False

    else:
        shutil.rmtree(staging_dir)
        print(f"       \033[94mINFO:\033[0m Staging directory deleted.")
        print('No naming issues requiring automatic correction found.')
        return False
