#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
 _    _      _        _                                            
| |_ | |__  (_)  ___ | | __ _ __    ___  ___  ___     _ __   _   _ 
| __|| '_ \ | | / __|| |/ /| '_ \  / _ \/ __|/ __|   | '_ \ | | | |
| |_ | | | || || (__ |   < | | | ||  __/\__ \\__ \ _ | |_) || |_| |
 \__||_| |_||_| \___||_|\_\|_| |_| \___||___/|___/(_)| .__/  \__, |
                                                     |_|     |___/ 
"""
# ==============================================================================================
#  File:        tes/libs/core/thickness.py
#  Project:     Transition Edge Sensor (TES) Microcalorimeter
#  Module:      Analysis functions for thickness analysis
#
#  Description:
#      This module analyzes and calculates thickness from roughness measurement data.
#      The roughness data waveform corresponds to a square wave (__|‾‾|__)
#      and a step function (__|‾‾ or ‾‾|__).
#
#  Requirements:
#      - Python 3.10+
#      - No external dependencies (Standard library only)
#
#  Author:      Ryota Fukuda (mailto:25la018c@rikkyo.ac.jp)
#  Maintainer:  Yamada Laboratory, Department of Physics, College of Science, Rikkyo University
#  Url:         https://github.com/Urchann/TES/tree/main/tes/libs/core/thickness.py
#  Created:     2025-10-09
#  Updated:     2025-12-31
#  Version:     2.1.0
#
#  ChangeLog:
#      v1.0.0 (2025-10-09)
#         - Establishment of thickness analysis logic for square wave roughness datasets.
#      v2.0.0 (2025-11-20)
#         - Addition of thickness analysis logic for step function roughness datasets.
#      v2.1.0 (2025-11-22)
#         - Added support for table overflow and page breaks in PDF output.
#
#  Notes:
#      - Analysis is executed in the 'tes' directory by running tes.py.
#      Execution command:
#          `python tes.py thickness <directory_name_containing_roughness_datasets>`
#      - The following options are supported:
#          > Excluding files from analysis.
#          > Specifying the output directory for results.
#          > Specifying the columns of the dataset to be used.
#          > Specifying the smoothing window for the roughness data.
#          > Specifying the waveform variation discrimination threshold.
#      - For detailed information, please see the help message:
#          `python tes.py thickness -h`
#      - If you analyze your own roughness data, change __maintainer__.
#
# ==============================================================================================


__author__     = 'Ryota Fukuda'
__maintainer__ = 'Ryota Fukuda'
__version__    = '2.1.0'

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'
import pandas as pd
import os
import glob
import sys
import re
import shutil
import textwrap
from scipy.ndimage import uniform_filter1d
from collections import defaultdict


class Thickness:
    def __init__(self, dir_input, dir_output, dat_exclude, dat_column, output_format, smooth_window, margin, boundary, threshold, mode):
        """
        Thicknessクラスのコンストラクタ。
        インスタンスを初期化し、解析設定およびファイルリストを作成する。

        Parameters:
            dir_input (str)    : 解析対象のデータ（*.txt）が含まれるディレクトリ名。
            dir_output (str)   : 解析結果を出力するディレクトリ名（デフォルトは'<dir_input>_result'）。
            dat_exclude (list) : 解析から除外するファイル名のリスト（デフォルトはなし）。
            dat_column (str)   : データファイル内で解析対象とする列名（デフォルトは'Normal'列）。
            output_format (str): 結果の出力形式（'html'(デフォルト)または'pdf'）。
            smooth_window (int): 移動平均スムージングのウィンドウサイズ（デフォルト値は50）。
            margin (int)       : ステップモードにおける立ち上がり・立ち下がり検出のマージン（デフォルト値は50）。
            boundary (float)   : ステップモードにおける勾配の境界検出の閾値係数（デフォルト値は0.1）。
            threshold (float)  : 矩形波モードにおける立ち上がり・立ち下がり検出の閾値係数（デフォルト値は0.5）。
            mode (str)         : 解析モード（'square'(デフォルト)または'step'）。

        Returns:
            None
        """
        self.dir_root = os.getcwd()
        self.dir_path = os.path.join(self.dir_root, 'data', 'roughness')
        self.dir_input_name = os.path.join(self.dir_path, dir_input)

        if not os.path.isdir(self.dir_input_name):
            raise FileNotFoundError(
                f"The specified input directory cannot be found: '{dir_input}'"
                f" -> Reference PATH: '{self.dir_input_name}'"
            )
        
        dat_all_list = sorted(glob.glob(os.path.join(self.dir_input_name, '*.txt')))
        dat_exclude  = set(dat_exclude) if dat_exclude else set()
        self.dat_input_list = [f for f in dat_all_list if os.path.basename(f) not in dat_exclude]

        if dir_output is None:
            dir_output_name = os.path.basename(self.dir_input_name) + '_result'
        else:
            dir_output_name = dir_output
        self.dir_output = os.path.join(self.dir_path, dir_output_name)

        print("==================================================")
        print(f" Module Version {__version__}")
        print(f" Reading Files in '{os.path.basename(self.dir_input_name)}' Directory.")
        print("==================================================")
        print(f" List of TXT Files to be analyzed ({len(self.dat_input_list)} files)")
        for i, datname in enumerate(self.dat_input_list):
            print(f' {i+1} : {os.path.basename(datname)}')
        print("--------------------------------------------------")

        self.dat_column    = dat_column
        self.output_format = output_format
        self.smooth_window = smooth_window
        self.margin        = margin
        self.threshold     = threshold
        self.boundary      = boundary
        self.mode          = mode

        self.dat_input_group  = defaultdict(lambda: defaultdict(list))
        self.dat_group        = defaultdict(lambda: defaultdict(dict))
        self.dat_smooth_group = defaultdict(lambda: defaultdict(dict))
        self.results          = defaultdict(lambda: defaultdict(dict))


    def _parse_filename_key(self, filename):
        """
        ファイル名を解析し、テスト番号と測定位置を抽出する関数。

        Parameters:
            filename (str): 解析対象のファイル名（拡張子込み）。

        Returns:
            part_test (str): テスト番号（test<数字>の型）。抽出できない場合はNone。
            part_pos (str) : 測定位置（数字のみ以外）。抽出できない場合はNone。
        """
        base_name = os.path.splitext(filename)[0]
        parts = base_name.split('_')
        if len(parts) < 2:
            return None, None
        
        part_test = None
        for part in parts:
            if re.match(r'^test\d+$', part):
                part_test = part
                break
        
        part_pos = parts[-2]
        if part_pos.isdigit():
            part_pos = None
        
        return part_test, part_pos


    def parse_filename(self):
        """
        入力ディレクトリ内のファイルについて、テスト番号と測定位置ごとにグループ化してメンバ変数に格納する関数。

        Parameters:
            None

        Returns:
            None
        """
        if self.dat_input_list:
            print(" --- File Grouping by Test Number & Position ---")
        else:
            sys.exit(" \033[31mNo files selected for analysis.\033[0m")
        
        for dat_name in self.dat_input_list:
            file_name = os.path.basename(dat_name)
            part_test, part_pos = self._parse_filename_key(file_name)

            
            if part_test is None:
                print(f" \033[33mWarning: No test number found in '{file_name}'. Substituting 'test1'.\033[0m")
                part_test = 'test1'
            elif part_pos is None:
                sys.exit(f" \033[31mFile name of '{file_name}' is too short. Check file naming rules. -> Reference PATH: {dat_name}\033[0m")
            
            self.dat_input_group[part_test][part_pos].append(dat_name)
    

    def open_file(self):
        """
        グループ化されたリストに基づき、実際のデータファイルを読み込み、解析用データとして格納する関数。

        Parameters:
            None

        Returns:
            None
        """
        if not self.dat_input_group:
            sys.exit(" \033[31mFATAL ERROR: No files were successfully grouped for analysis. Check file names.\033[0m")
        
        for test_group, pos_group in self.dat_input_group.items():
            print(f"    [Test Number: {test_group}]")
            for pos, files in pos_group.items():
                print(f"       > Pos. {pos} ({len(files)} Files)")
                for file in files:
                    file_name = os.path.basename(file)
                    header = 0
                    try:
                        with open(file, 'r') as f:
                            for line in f:
                                header += 1
                                if self.dat_column in line:
                                    break
                    except Exception as e:
                        print(f"          Warning: Column '{self.dat_column}' not found in header. Skipping. Error: {e}")
                        continue

                    try:
                        df = pd.read_csv(file, sep=r'\s+', skiprows=header-1, skipinitialspace=True)
                        self.dat_group[test_group][pos][file] = df[self.dat_column]
                        print(f"          - {file_name} (Data Shape: {df.shape})")
                    except Exception as e:
                        print(f"           Error: Failed to load data for '{file_name}'. Error: {e}")
                        continue
        print("--------------------------------------------------")
    

    def _find_extrema_indices(self, dat_smooth):
        """
        取得データをスムージングしたものから、立ち上がりと立ち下がりの主要な変化点インデックスを検出する関数。
        squareモードでは、隣接要素の差分のthresholdを超える点を取り、
        stepモードでは、隣接要素の中心差分の最大値の点を取る。

        Parameters:
            dat_smooth (ndarray): 取得データをスムージング済みの1次元データ配列。

        Returns:
            rise_idx (list)     : 立ち上がりが検出されたインデックスのリスト。
            fall_idx (list)     : 立ち下がりが検出されたインデックスのリスト。
        """
        dat_diff = np.diff(dat_smooth)
        dat_grad = np.gradient(dat_smooth)
        rise_idx = []
        fall_idx = []
        
        if self.mode == 'square':
            max_diff = np.max(dat_diff)
            rising   = False
            rise_threshold = max_diff * self.threshold
            for i in range(1, len(dat_diff)):
                if dat_diff[i] > rise_threshold and not rising:
                    rise_idx.append(i)
                    rising = True
                elif dat_diff[i] < rise_threshold and rising:
                    rising = False
            
            min_diff = np.min(dat_diff)
            falling  = False
            fall_threshold = min_diff * self.threshold
            for i in range(1, len(dat_diff)):
                if dat_diff[i] < fall_threshold and not falling:
                    fall_idx.append(i)
                    falling = True
                elif dat_diff[i] > fall_threshold and falling:
                    falling = False
        
        if self.mode == 'step':
            rise_grad_idx = np.argmax(dat_grad)
            rise_grad_val = dat_grad[rise_grad_idx]
            fall_grad_idx = np.argmin(dat_grad)
            fall_grad_val = dat_grad[fall_grad_idx]
            if np.abs(rise_grad_val) > np.abs(fall_grad_val):
                rise_idx.append(int(rise_grad_idx))
            else:
                fall_idx.append(int(fall_grad_idx))
                
        return rise_idx, fall_idx
    

    def _find_extrema_boundary(self, dat_smooth, rise_idx, fall_idx):
        """
        変化点周辺を探索し、立ち上がり・立ち下がりの開始点と終了点を決定する関数。
        データの勾配や値の増減を利用して境界を特定する。

        Parameters:
            dat_smooth (ndarray): 取得データをスムージング済みの1次元データ配列。
            rise_idx (list)     : 立ち上がり変化点のインデックスのリスト。
            fall_idx (list)     : 立ち下がり変化点のインデックスのリスト。

        Returns:
            pre_rise_idx (list) : 立ち上がり開始点のリスト。
            post_rise_idx (list): 立ち上がり終了点のリスト。
            pre_fall_idx (list) : 立ち下がり開始点のリスト。
            post_fall_idx (list): 立ち下がり終了点のリスト。
        """
        N = len(dat_smooth)
        dat_grad = np.gradient(dat_smooth)
        grad_abs = np.abs(dat_grad)
        pre_rise_idx  = []
        post_rise_idx = []
        pre_fall_idx  = []
        post_fall_idx = []
        
        if self.mode == 'square':
            for idx in rise_idx:
                pre_idx = idx
                for i in range(idx, 0, -1):
                    if dat_smooth[i] <= dat_smooth[i-1]:
                        pre_idx = i
                        break
                pre_rise_idx.append(pre_idx)
                
                post_idx = idx
                for i in range(idx, N-1, 1):
                    if dat_smooth[i] >= dat_smooth[i+1]:
                        post_idx = i
                        break
                post_rise_idx.append(post_idx)
            
            for idx in fall_idx:
                pre_idx = idx
                for i in range(idx, 0, -1):
                    if dat_smooth[i] >= dat_smooth[i-1]:
                        pre_idx = i
                        break
                pre_fall_idx.append(pre_idx)
                
                post_idx = idx
                for i in range(idx, N-1, 1):
                    if dat_smooth[i] <= dat_smooth[i+1]:
                        post_idx = i
                        break
                post_fall_idx.append(post_idx)

        if self.mode == 'step':
            for idx in rise_idx:
                max_grad_val  = grad_abs[idx]
                threshold_val = max_grad_val * self.boundary
                
                start = idx
                while start > 0:
                    if grad_abs[start] < threshold_val:
                        break
                    start -= 1
                pre_rise_idx.append(int(max(0, start - self.margin)))

                end = idx
                while end < N - 1:
                    if grad_abs[end] < threshold_val:
                        break
                    end += 1
                post_rise_idx.append(int(min(N - 1, end + self.margin)))
            
            for idx in fall_idx:
                max_grad_val  = grad_abs[idx]
                threshold_val = max_grad_val * self.boundary
                
                start = idx
                while start > 0:
                    if grad_abs[start] < threshold_val:
                        break
                    start -= 1
                pre_fall_idx.append(int(max(0, start - self.margin)))
                
                end = idx
                while end < N - 1:
                    if grad_abs[end] < threshold_val:
                        break
                    end += 1
                post_fall_idx.append(int(min(N - 1, end + self.margin)))
                
        return pre_rise_idx, post_rise_idx, pre_fall_idx, post_fall_idx


    def calculate_thickness(self):
        """
        読み込んだ全データに対して厚み計算を実行する関数。
        取得データのスムージングから変化点検出、境界特定を行い、
        立ち上がり・立ち下がりの段差（厚み）およびその誤差（標準偏差等）を算出・集計する。

        Parameters:
            None

        Returns:
            None
        """
        if self.dat_input_list:
            print(" --- Run Thickness Calculation ---")
        else:
            sys.exit(" \033[31mFATAL ERROR: No data loaded for analysis.\033[0m")

        for test_group, pos_group in self.dat_group.items():
            print(f" {test_group} Analysis")
            for pos, dat_dict in pos_group.items():
                print(f"    Pos. {pos}")
                pos_rise_height = []
                pos_fall_height = []
                
                for i, (file_path, dat) in enumerate(dat_dict.items()):
                    print(f'       {i+1}. ', end='')
                    filename = os.path.basename(file_path)
                    
                    if len(dat) < 3:
                        print(f"\033[33mSkipping '{filename}': Data too short.\033[0m")
                        continue

                    dat_smooth = uniform_filter1d(dat.values, size=self.smooth_window)
                    self.dat_smooth_group[test_group][pos][file_path] = dat_smooth
                    
                    rise_idx, fall_idx = self._find_extrema_indices(dat_smooth)
                    if not rise_idx and not fall_idx:
                        print(f"\033[33mSkipping '{filename}': No clear rise/fall events found.\033[0m")
                        continue

                    pre_rise_idx, post_rise_idx, pre_fall_idx, post_fall_idx = self._find_extrema_boundary(dat_smooth, rise_idx, fall_idx)
                    avg_rise_height     = np.nan
                    avg_fall_height     = np.nan
                    avg_rise_height_err = np.nan
                    avg_fall_height_err = np.nan

                    if pre_rise_idx and post_rise_idx:
                        rise_height     = np.abs(dat_smooth[post_rise_idx] - dat_smooth[pre_rise_idx])
                        avg_rise_height = np.mean(rise_height)
                        
                        if len(rise_height) < 2:
                            err_base_idx = pre_rise_idx[0] if pre_rise_idx[0] > 0 else 1
                            avg_rise_height_err = np.std(dat[0:err_base_idx], ddof=1)
                        else:
                            avg_rise_height_err = np.std(rise_height, ddof=1)
                        
                        pos_rise_height.append(avg_rise_height)

                    if pre_fall_idx and post_fall_idx:
                        fall_height = np.abs(dat_smooth[pre_fall_idx] - dat_smooth[post_fall_idx])
                        avg_fall_height = np.mean(fall_height)

                        if len(fall_height) < 2:
                            if pre_rise_idx:
                                err_base_idx = pre_rise_idx[0]
                            else:
                                err_base_idx = pre_fall_idx[0]
                            
                            err_base_idx = err_base_idx if err_base_idx > 0 else 1
                            avg_fall_height_err = np.std(dat[0:err_base_idx], ddof=1)
                        else:
                            avg_fall_height_err = np.std(fall_height, ddof=1)
                        
                        pos_fall_height.append(avg_fall_height)
                    
                    self.results[test_group][pos][file_path] = {
                        'rise_thickness'      : avg_rise_height,
                        'fall_thickness'      : avg_fall_height,
                        'rise_thickness_error': avg_rise_height_err,
                        'fall_thickness_error': avg_fall_height_err,
                        'pre_rise_idx'        : pre_rise_idx,
                        'post_rise_idx'       : post_rise_idx,
                        'pre_fall_idx'        : pre_fall_idx,
                        'post_fall_idx'       : post_fall_idx
                    }
                    
                    if pre_rise_idx and post_rise_idx and pre_fall_idx and post_fall_idx:
                        rise_height = np.abs(dat_smooth[post_rise_idx] - dat_smooth[pre_rise_idx])
                        fall_height = np.abs(dat_smooth[pre_fall_idx] - dat_smooth[post_fall_idx])
                        avg_rise_height = np.mean(rise_height)
                        avg_fall_height = np.mean(fall_height)

                        if len(rise_height) < 2 or len(fall_height) < 2:
                            avg_rise_height_err = np.std(dat[0:pre_rise_idx[0]], ddof=1) / np.sqrt(pre_rise_idx[0])
                            avg_fall_height_err = np.std(dat[0:pre_rise_idx[0]], ddof=1) / np.sqrt(pre_rise_idx[0])
                        else:
                            avg_rise_height_err = np.std(rise_height, ddof=1) / np.sqrt(len(rise_height))
                            avg_fall_height_err = np.std(fall_height, ddof=1) / np.sqrt(len(fall_height))
                        
                        pos_rise_height.append(avg_rise_height)
                        pos_fall_height.append(avg_fall_height)
                        
                    self.results[test_group][pos][file_path] = {
                        'rise_thickness'      : avg_rise_height,
                        'fall_thickness'      : avg_fall_height,
                        'rise_thickness_error': avg_rise_height_err,
                        'fall_thickness_error': avg_fall_height_err,
                        'pre_rise_idx'        : pre_rise_idx,
                        'post_rise_idx'       : post_rise_idx,
                        'pre_fall_idx'        : pre_fall_idx,
                        'post_fall_idx'       : post_fall_idx
                    }
                    
                    r_str = f"{avg_rise_height:.2f} ± {avg_rise_height_err:.2f}" if not np.isnan(avg_rise_height) else "---"
                    f_str = f"{avg_fall_height:.2f} ± {avg_fall_height_err:.2f}" if not np.isnan(avg_fall_height) else "---"
                    print(f"Calculated '{filename}': Rise Height = {r_str} Å, Fall Height = {f_str} Å")

                if pos_rise_height or pos_fall_height:
                    pos_avg_rise = np.mean(pos_rise_height) if pos_rise_height else np.nan
                    pos_avg_fall = np.mean(pos_fall_height) if pos_fall_height else np.nan
                    
                    if len(pos_rise_height) > 1:
                        pos_err_rise = np.std(pos_rise_height, ddof=1) / np.sqrt(len(pos_rise_height))
                    elif len(pos_rise_height) == 1:
                        pos_err_rise = avg_rise_height_err 
                    else:
                        pos_err_rise = np.nan

                    if len(pos_fall_height) > 1:
                        pos_err_fall = np.std(pos_fall_height, ddof=1) / np.sqrt(len(pos_fall_height))
                    elif len(pos_fall_height) == 1:
                        pos_err_fall = avg_fall_height_err
                    else:
                        pos_err_fall = np.nan

                    self.results[test_group][pos]['pos_thickness'] = {
                        'avg_rise_thickness': pos_avg_rise,
                        'avg_fall_thickness': pos_avg_fall,
                        'avg_rise_thickness_error': pos_err_rise,
                        'avg_fall_thickness_error': pos_err_fall
                    }
                    
                    r_res = f"{pos_avg_rise:.2f}" if not np.isnan(pos_avg_rise) else "---"
                    f_res = f"{pos_avg_fall:.2f}" if not np.isnan(pos_avg_fall) else "---"
                    print(f" -> Average Thickness: \033[30;47mRise Height = {r_res} Å, Fall Height = {f_res} Å\033[0m")
                else:
                    print(" \033[31mNothing is Calculated in this Position.\033[0m")
                    continue
        print('--------------------------------------------------')


    def plot_roughtness(self):
        """
        各測定データの波形プロットを作成し、画像ファイルとして保存する関数。
        元データ、スムージングデータ、および検出された立ち上がり・立ち下がり位置を描画する。

        Parameters:
            None

        Returns:
            None
        """
        for test_group, pos_group in self.dat_group.items():
            for pos, dat_dict in pos_group.items():
                for file_path, dat in dat_dict.items():
                    filename = os.path.basename(file_path)
                    output_file = os.path.splitext(filename)[0] + ".png" 
                    output_path = os.path.join(self.dir_output, output_file)

                    dat_smooth = self.dat_smooth_group[test_group][pos].get(file_path)
                    result = self.results[test_group][pos].get(file_path)

                    if dat_smooth is None or result is None:
                        continue

                    pre_rise_idx  = result.get('pre_rise_idx') or []
                    post_rise_idx = result.get('post_rise_idx') or []
                    pre_fall_idx  = result.get('pre_fall_idx') or []
                    post_fall_idx = result.get('post_fall_idx') or [] 
                    rise_height   = result.get('rise_thickness', np.nan)
                    fall_height   = result.get('fall_thickness', np.nan)

                    plt.figure(figsize=(12, 6))
                    plt.plot(dat.index, dat, label='Original Data', alpha=0.7, lw=0.8)
                    plt.plot(dat.index, dat_smooth, label=f"Smoothed Data (Window: {self.smooth_window})", c='orange', lw=1.5)

                    if pre_rise_idx and post_rise_idx:
                        plt.plot(pre_rise_idx, dat_smooth[pre_rise_idx], 'r^', ms=8)
                        plt.plot(post_rise_idx, dat_smooth[post_rise_idx], 'rv', ms=8)

                        if not np.isnan(rise_height):
                            plt.text(0.01, 0.98,
                                    f"Rise Thickness: {rise_height:.2f} Å",
                                    transform=plt.gca().transAxes,
                                    ha='left', va='top',
                                    color='darkgreen',
                                    fontsize=12,
                                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='lightgray'))
                            
                    if pre_fall_idx and post_fall_idx:
                        plt.plot(pre_fall_idx, dat_smooth[pre_fall_idx], 'bv', ms=8)
                        plt.plot(post_fall_idx, dat_smooth[post_fall_idx], 'b^', ms=8)

                        text_y = 0.92 if (pre_rise_idx and post_rise_idx) else 0.98
                        if not np.isnan(fall_height):
                            plt.text(0.01, text_y,
                                    f"Fall Thickness: {fall_height:.2f} Å",
                                    transform=plt.gca().transAxes,
                                    ha='left', va='top',
                                    color='darkred',
                                    fontsize=12,
                                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='lightgray'))

                    plt.title(f"Roughness Analysis: {filename}")
                    plt.xlabel("Measurement Index")
                    plt.ylabel("Roughness (Å)")
                    plt.legend(loc=1)
                    plt.grid(True, linestyle=':', alpha=0.6)
                    plt.tight_layout()
                    plt.savefig(output_path, dpi=300) 
                    plt.close()

                    self.results[test_group][pos][file_path].update({'plot_roughtness': output_file})

        
    def plot_test_correlation(self):
        """
        テスト番号後とに、各測定点の厚み（立ち上がり・立ち下がり）の相関プロットを作成する関数。
        横軸にファイル識別子、縦軸に計算された厚みを表示する。

        Parameters:
            None

        Returns:
            None
        """
        for test_group, pos_group in self.results.items():
            points = []
            rise_height = []
            fall_height = []
            rise_height_err = []
            fall_height_err = []
            
            for _, result_dict in pos_group.items():
                for file_path, results in result_dict.items():
                    if 'rise_thickness' not in results:
                        continue
                    
                    filename = os.path.splitext(file_path)[0]
                    point_label = '_'.join(filename.split('_')[-3:])
                    points.append(point_label)

                    rise_height.append(results['rise_thickness'])
                    fall_height.append(results['fall_thickness'])
                    rise_height_err.append(results.get('rise_thickness_error', 0))
                    fall_height_err.append(results.get('fall_thickness_error', 0))
                    
            output_file = f'thickness_correlation_{test_group}.png'
            output_path = os.path.join(self.dir_output, output_file)

            self.test_correlation_map = self.test_correlation_map if hasattr(self, 'test_correlation_map') else {}
            self.test_correlation_map[test_group] = output_file
            
            plt.figure(figsize=(8, 6))
            plt.errorbar(points, rise_height, yerr=rise_height_err, color='red',  capsize=3, markersize=5, marker='o', markeredgecolor='red',  markerfacecolor='red',  ls='None', label='Rise Thickness (Å)')
            plt.errorbar(points, fall_height, yerr=fall_height_err, color='blue',  capsize=3, markersize=5, marker='o', markeredgecolor='blue',  markerfacecolor='blue',  ls='None', label='Fall Thickness (Å)')
            plt.title(f"Thickness Correlation for {test_group} measurement points")
            plt.xlabel("Measurement Points")
            plt.ylabel("Thickness (Å)")
            plt.xticks(rotation=90)
            plt.legend()
            plt.grid(True, linestyle=':', alpha=0.6)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
    

    def plot_pos_correlation(self, pos_list):
        """
        測定位置ごとの平均厚みの相関プロットを作成する関数。
        指定された位置リストに基づき、位置間の厚みの変化を可視化する。

        Parameters:
            pos_list (list): プロット対象とする測定位置キーのリスト。Noneの場合は全位置。

        Returns:
            None
        """
        for test_group, pos_group in self.results.items():
            points = []
            rise_height = []
            fall_height = []
            rise_height_err = []
            fall_height_err = []
            
            iterator = pos_list if pos_list is not None else pos_group.keys()
            for pos in iterator:
                if pos not in pos_group:
                    continue
                
                results = pos_group[pos]
                if 'pos_thickness' in results:
                    pos_results = results['pos_thickness']
                    points.append(pos)
                    rise_height.append(pos_results['avg_rise_thickness'])
                    fall_height.append(pos_results['avg_fall_thickness'])
                    rise_height_err.append(pos_results.get('avg_rise_thickness_error', 0))
                    fall_height_err.append(pos_results.get('avg_fall_thickness_error', 0))

            if self.mode == 'step':
                pos_str = '_'.join(points)
                sorted_indices = np.arange(len(points))
            if self.mode == 'square':
                pos_str = '_'.join(sorted(points))
                sorted_indices = np.argsort([re.sub(r'[^0-9]', '', p) for p in points if re.sub(r'[^0-9]', '', p)])
            output_file = f"thickness_correlation_{test_group}_{pos_str}.png"
            output_path = os.path.join(self.dir_output, output_file)

            points = np.array(points)[sorted_indices]
            rise_height = np.array(rise_height)[sorted_indices]
            fall_height = np.array(fall_height)[sorted_indices]
            rise_height_err = np.array(rise_height_err)[sorted_indices]
            fall_height_err = np.array(fall_height_err)[sorted_indices]
            
            self.pos_correlation_map = self.pos_correlation_map if hasattr(self, 'pos_correlation_map') else {}
            self.pos_correlation_map[f"{test_group}_{pos_str}"] = output_file
            
            x_indices = np.arange(len(points))
            plt.figure(figsize=(10, 6))
            plt.errorbar(x_indices, rise_height, yerr=rise_height_err, color='red', capsize=3, markersize=5, marker='o', markeredgecolor='red', markerfacecolor='red', label='Rise Thickness (Å)')
            plt.errorbar(x_indices, fall_height, yerr=fall_height_err, color='blue', capsize=3, markersize=5, marker='o', markeredgecolor='blue', markerfacecolor='blue', label='Fall Thickness (Å)')
            plt.xticks(x_indices, points)
            plt.title(f"Thickness correlation for {test_group} measurement points")
            plt.xlabel("Measurement Position")
            plt.ylabel("Thickness (Å)")
            plt.legend()
            plt.grid(True, linestyle=':', alpha=0.6)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
        

    def export_results(self, board_img, *groups):
        """
        全解析結果、プロット画像、集計テーブルをまとめてHTML形式ないしはLaTeX (PDF) 形式で出力する関数。
        PDF形式での出力を望む際は、.texファイルになるため、別途コンパイルが必要。
        
        Parameters:
            board_img (str): 厚み解析の相関結果における、基板の概略図（測定位置図）の画像ファイルパス。
            *groups (list) : 厚み解析の相関結果における、比較・集計対象とする測定位置グループ（可変長引数）。
        """
        results    = []
        ea_summary = {}
        for test_group, pos_group in self.results.items():
            for pos, dat_dict in pos_group.items():
                pos_summary = dat_dict.get('pos_thickness')
                if pos_summary:
                    ea_summary[f"{test_group}_{pos}"] = {
                        'ea_ave_rise': pos_summary['avg_rise_thickness'],
                        'ea_err_rise': pos_summary['avg_rise_thickness_error'],
                        'ea_ave_fall': pos_summary['avg_fall_thickness'],
                        'ea_err_fall': pos_summary['avg_fall_thickness_error']
                    }
                
                for file_path, result in dat_dict.items():
                    if 'rise_thickness' in result:
                        results.append({
                            'InputFile'      : os.path.basename(file_path),
                            'RiseHeight'     : result.get('rise_thickness'),
                            'RiseHeightError': result.get('rise_thickness_error'),
                            'FallHeight'     : result.get('fall_thickness'),
                            'FallHeightError': result.get('fall_thickness_error'),
                            'EachPlotResult' : result.get('plot_roughtness', '')
                        })

        test_correlations = self.test_correlation_map if hasattr(self, 'test_correlation_map') else {}
        ea_correlations   = self.pos_correlation_map if hasattr(self, 'pos_correlation_map') else {}
        test_nums = sorted(list(set([k.split('_')[0] for k in ea_summary.keys() if 'test' in k])))

        def _safe_sort_key(s):
            match = re.search(r'\d+', s)
            if match:
                return (0, int(match.group()))
            else:
                return (1, s)
        ea_nums = sorted({k.split('_')[1] for k in ea_summary.keys()}, key=_safe_sort_key)
        if not results:
            print("\n\033[31mError: No valid data found for export.\033[0m")
            return
        
        pos_group_tuple = [tuple(item.split('-')) for item in groups]
        desc_summary = textwrap.indent(
            textwrap.dedent("""\
                <p>　理研CRの表面粗さ計（触針型）による測定データから，試料の厚みを計算する。
                   <br>　粗さデータセットの波形は矩形波上になっており，本解析では，その矩形波の凸の上昇部分を「立ち上がり」，下降部分を「立ち下がり」と呼ぶ。"""\
                   """この立ち上がりと立ち下がりの波高値をそれぞれ計算することにより，試料の厚みを評価する。
                   <br>　波高値の計算方法は，まず素データをスムージングした後，隣り合うデータ点の差を取り，その最大値と最小値をそれぞれ立ち上がり点および立ち下がり点と決定する。"""\
                   """次に，立ち上がり点からデータ始点，立ち下がり点からデータ終点までについてスキャンし，それぞれデータ点が増加する点を立ち上がり始点と立ち下がり始点とする。"""\
                   """さらに，立ち上がり点から立ち下がり点，立ち下がり点から立ち上がり点までについてスキャンし，それぞれデータ点が減少する点を立ち上がり終点と立ち下がり終点とする。"""\
                   """最後に，立ち上がりおよび立ち下がりの始点と終点の差を計算することで，試料の厚みをそれぞれ評価する。
                </p>"""),
        "   "
        )
        desc_plot = textwrap.indent(
            textwrap.dedent("""\
                <p>　全測定データを各々プロットする。
                   <br>　横軸は測定位置のインデックスで，縦軸が膜厚値である。青色のラインは素データ，橙色のラインはuniform_filter1dを用いてスムージングしたデータを表す。"""\
                   """また，赤色の上矢印マーカーは立ち上がり始めの点，赤色の下矢印マーカーは立ち上がり終わりの点，青色の下矢印マーカーは立ち下がりは始めの点，青色の上矢印マーカーは立ち下がり終わりの点を表す。
                </p>"""),
        "      "
        )
        desc_test = textwrap.indent(
            textwrap.dedent("""\
                <p>　基板ごとに，立ち上がりと立ち下がりによる膜厚値の全計算結果をまとめてプロットする。
                   <br>　横軸は測定点で，縦軸が膜厚値である。また，赤色の誤差付きマーカーは立ち上がりの膜厚値，青色の誤差付きマーカーは立ち下がりの膜厚値を表す。
                   <br>　誤差は，各データについて，Measurement Indexにおける0から立ち上がり始めまでのデータ点から，不偏標準偏差で計算する。
                </p>"""),
        "      "
        )
        desc_ea = textwrap.indent(
            textwrap.dedent(f"""\
                <p>　基板ごとに，測定位置による相関をプロットする。以下，上の図は基板の測定位置（{pos_group_tuple[0][0][0]}番号）の詳細。"""\
                   f"""相関を見る位置は，({', '.join(pos_group_tuple[0])}) と ({', '.join(pos_group_tuple[1])}) の2つとし，以下の下の図にまとめる。
                   <br>　誤差は，各データについて，Measurement Indexにおける0から立ち上がり始めまでのデータ点から，不偏標準偏差で計算する。
                </p>"""),
        "         "
        )
        desc_table = textwrap.indent(
            textwrap.dedent("""\
                <p>　全測定データにおける膜厚値の表にまとめる。
                   <br>　1列目から順に，データファイル名，立ち上がり始めの値，立ち上がり終わりの値，立ち上がりの膜厚値，立ち下がり始めの値，立ち下がり終わりの値，立ち下がりの膜厚値である。
                </p>"""),
        "      "
        )
        desc_ave_table = textwrap.indent(
            textwrap.dedent("""\
                <p>　各基板の位置ごとに，立ち上がり時間と立ち下がり時間による膜厚値の結果（平均値）を表にまとめる。
                   <br>　誤差は，平均値の不偏標準偏差で計算する。
                </p>"""),
        "         "
        )
        description = [desc_summary, desc_plot, desc_test, desc_ea, desc_table, desc_ave_table]
        
        if os.path.exists(board_img):
            shutil.copy(board_img, self.dir_output)
            print(f"\n BoardSchematic.png has been copied to {self.dir_output}.")
        else:
            print(f"\n \033[93mWARNING:\033[0m BoardSchematic.png not found. Create HTML or PDF without board schematic.")
    
        if self.output_format == 'html':
            output_path = os.path.join(self.dir_output, "thickness_results.html")
            output_dir  = os.path.basename(self.dir_output)
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write("<!DOCTYPE html>\n")
                f.write("<html lang='ja'>\n")
                f.write("<head>\n")
                f.write("   <meta charset='UTF-8'>\n")
                f.write("   <meta name='viewport' content='width=device-width, initial-scale=1.0'>\n")
                f.write("   <title>Roughness Analysis Results</title>\n")
                f.write("   <style>\n")
                f.write("      body { font-family: 'Times New Roman', Times, serif; margin: 20px; line-height: 1.6; }\n")
                f.write("      h1 { font-size: 2rem; color: #333; }\n")
                f.write("      h2 { font-size: 1.5rem; color: #333; }\n")
                f.write("      h3 { font-size: 1rem; color: #333; }\n")
                f.write("      table { width: 100%; border-collapse: collapse; margin-top: 20px; }\n")
                f.write("      th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n")
                f.write("      th { background-color: #f2f2f2; text-align: center; }\n")
                f.write("      th:first-child { text-align: left; }\n")
                f.write("      td { text-align: right; }\n")
                f.write("      td:first-child { text-align: left; }\n")
                f.write("      .plot-section { margin-top: 40px; border-top: 1px solid #ccc; padding-top: 20px; }\n")
                f.write("      .plot-section img { max-width: 50%; height: auto; display: block; margin: auto; border: 1px solid #eee; box-shadow: 2px 2px 5px rgba(0,0,0,0.1); }\n")
                f.write("      .plot-section a { display: inline-block; }\n")
                f.write("      .plot-section .plot-container { text-align: center; }\n")
                f.write("      .correlation-grid-container { overflow-x: auto; }\n")
                f.write("      .correlation-grid { display: grid; grid-template-columns: repeat(2, 1fr); gap: 20px; justify-items: center; align-items: center; }\n")
                f.write("      .correlation-grid img { max-width: 80%; height: auto; }\n")
                f.write("      .summary-box { border: 1px solid #ccc; padding: 15px; margin: 15px 0; background-color: #f9f9f9; border-radius: 5px; }\n")
                f.write("      .summary-item { margin-bottom: 10px; }\n")
                f.write("      .summary-table-container { overflow-x: auto; }\n")
                f.write("      .summary-table thead th { background-color: #e6e6e6; border: 1px solid #b0b0b0; text-align: center; }\n")
                f.write("      .summary-table tbody td { border: 1px solid #b0b0b0; text-align: center; white-space: nowrap; }\n")
                f.write("      .footer { text-align: center; margin-top: 50px; font-size: 0.8em; color: #777; }\n")
                f.write("      hr { border: none; height: 1px; background-color: #ccc; margin: 20px 0; }\n")
                f.write("   </style>\n")
                f.write("   <script src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>\n")
                f.write("</head>\n")
                f.write("<body>\n")
                f.write(f"   <h1>Roughness Analysis Results for {'_'.join(output_dir.split('_')[:-1])}</h1>\n")
                f.write(f"   <p><strong>Author:</strong> {__maintainer__}</p>\n")
                f.write(description[0] + "\n")
                
                f.write("   <div class='plot-section'>\n")
                f.write("      <h2>Thickness Results</h2>\n")
                f.write("      <div class='summary-box'>\n")
                f.write("         <h3>Correlation Plot</h3>\n")
                f.write(description[3] + "\n")
                f.write("         <div class='plot-container'>\n")
                f.write(f"            <a href='./{os.path.basename(board_img)}' target='_blank'><img src='./{os.path.basename(board_img)}' alt='Measurement position on the board'></a>\n")
                f.write("         </div>\n")
                f.write("         <div class='correlation-grid'>\n")
                def _sort_ea_key(key):
                    group  = {}
                    match  = re.search(r'test(\d+)', key)
                    number = int(match.group(1)) if match else float('inf')
                    for i in range(len(groups)):
                        group[groups[i]] = i
                    suffix = key.split('_')[-1]
                    order  = group.get(suffix, float('inf'))
                    return (number, order)
                sorted_ea_keys = sorted(ea_correlations.keys(), key=_sort_ea_key)
                for group_key in sorted_ea_keys:
                    plot_file = ea_correlations[group_key]
                    if plot_file:
                        f.write("            <div class='plot-container'>\n")
                        f.write(f"               <a href='./{plot_file}' target='_blank'><img src='./{plot_file}' alt='Correlation of thickness for {group_key}'></a>\n")
                        f.write(f"               <p>{group_key.split('_')[0]+": "+"-".join(group_key.split('_')[1:])}</p>\n")
                        f.write("            </div>\n")
                f.write("         </div>\n")

                f.write("         <h3>Average Thickness Table</h3>\n")
                f.write(description[5] + "\n")
                f.write("         <div class='summary-table-container'>\n")
                f.write("            <table class='summary-table' style='width: 100%; border-collapse: collapse; margin-top: 20px;'>\n")
                f.write("               <thead>\n")
                f.write("                  <tr>\n")
                f.write("                     <th rowspan='2'></th>\n")
                for test_name in test_nums:
                    f.write(f"                     <th colspan='{len(ea_nums)}'>{test_name}</th>\n")
                f.write("                  </tr>\n")
                f.write("                  <tr>\n")
                for _ in test_nums:
                    for ea_name in ea_nums:
                        f.write(f"                     <th>{ea_name}</th>\n")
                f.write("                  </tr>\n")
                f.write("               </thead>\n")
                f.write("               <tbody>\n")
                f.write("                  <tr>\n")
                f.write("                     <th style='background-color: #e6e6e6; border: 1px solid #b0b0b0;'>Rise (Å)</th>\n")
                for test_name in test_nums:
                    for ea_name in ea_nums:
                        key = f"{test_name}_{ea_name}"
                        summary_data = ea_summary.get(key)
                        if summary_data:
                            rise_val = f"{summary_data['ea_ave_rise']:.2f} \\(\\pm\\) {summary_data['ea_err_rise']:.2f}"
                        else:
                            rise_val = '---'
                        f.write(f"                     <td>{rise_val}</td>\n")
                f.write("                  </tr>\n")
                f.write("                  <tr>\n")
                f.write("                     <th style='background-color: #e6e6e6; border: 1px solid #b0b0b0;'>Fall (Å)</th>\n")
                for test_name in test_nums:
                    for ea_name in ea_nums:
                        key = f"{test_name}_{ea_name}"
                        summary_data = ea_summary.get(key)
                        if summary_data:
                            fall_val = f"{summary_data['ea_ave_fall']:.2f} \\(\\pm\\) {summary_data['ea_err_fall']:.2f}"
                        else:
                            fall_val = '---'
                        f.write(f"                     <td>{fall_val}</td>\n")
                f.write("                  </tr>\n")
                f.write("               </tbody>\n")
                f.write("            </table>\n")
                f.write("         </div>\n")
                f.write("      </div>\n")
                f.write("   </div>\n")
                f.write("   <br>\n")

                f.write("   <div class='plot-section'>\n")
                f.write("      <h2>Thickness Results Table</h2>\n")
                f.write(description[4] + "\n")
                f.write("      <table>\n")
                f.write("         <thead>\n")
                f.write("            <tr><th>File Name</th><th>Rise Height (Å)</th><th>Fall Height (Å)</th></tr>\n")
                f.write("         </thead>\n")
                f.write("         <tbody>\n")
                for r in results:
                    rise_height       = r['RiseHeight'] if r['RiseHeight'] is not None else '---'
                    fall_height       = r['FallHeight'] if r['FallHeight'] is not None else '---'
                    rise_height_error = r['RiseHeightError'] if r['RiseHeightError'] is not None else '---'
                    fall_height_error = r['FallHeightError'] if r['FallHeightError'] is not None else '---'
                    f.write(f"            <tr><td>{r['InputFile']}</td>"
                            f"<td>{rise_height:.2f} \\(\\pm\\) {rise_height_error:.2f}</td>"
                            f"<td>{fall_height:.2f} \\(\\pm\\) {fall_height_error:.2f}</td>"
                            f"</tr>\n")
                f.write("         </tbody>\n")
                f.write("      </table>\n")
                f.write("   </div>\n")

                f.write("   <br><br><br>\n")
                f.write("   <div class='plot-section'>\n")
                f.write("      <h2>Correlation Results</h2>\n")
                f.write(description[2] + "\n")
                for i, (group_name, correlation_file) in enumerate(test_correlations.items()):
                    if i != 0:
                        f.write("      <hr>\n")
                    f.write(f"      <h3>{group_name} Correlation</h3>\n")
                    f.write("      <div class='plot-container'>\n")
                    f.write(f"         <a href='./{correlation_file}' target='_blank'><img src='./{correlation_file}' alt='Correlation of thickness for {group_name}'></a>\n")
                    f.write("      </div>\n")
                f.write("   </div>\n")
                
                f.write("   <div class='plot-section'>\n")
                f.write("      <h2>Each Plot</h2>\n")
                f.write(description[1] + "\n")
                for i, r in enumerate(results):
                    plot_png_path = r['EachPlotResult']
                    if i != 0:
                        f.write("      <hr>\n")
                    f.write(f"      <h3>{r['InputFile']}</h3>\n")
                    f.write("      <div class='plot-container'>\n")
                    f.write(f"         <a href='./{plot_png_path}' target='_blank'><img src='./{plot_png_path}' alt='Roughness Plot for {r['InputFile']}'></a>\n")
                    f.write("      </div>\n")
                f.write("   </div>\n")
                
                f.write("</body>\n")
                f.write("</html>\n")
            print(f" HTML saved to: {output_path}")

        elif self.output_format == 'pdf':
            output_path = os.path.join(self.dir_output, "thickness_results.tex")
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(r"\RequirePackage{plautopatch}" + "\n")
                f.write(r"\RequirePackage[l2tabu,orthodox]{nag}" + "\n")
                f.write(r"\documentclass[platex,dvipdfmx,10pt,twoside,a4paper,jis2004]{jsarticle}" + "\n")
                f.write(r"\usepackage[top=4cm,bottom=3cm,left=2cm,right=2cm]{geometry}" + "\n")
                f.write("\n")
                f.write(r"\usepackage[deluxe]{otf}" + "\n")
                f.write(r"\usepackage[T1]{fontenc}" + "\n")
                f.write(r"\usepackage{lmodern}" + "\n")
                f.write(r"\usepackage{textcomp}" + "\n")
                f.write(r"\usepackage[geometry,electronic,weather,clock,alpine,misc]{ifsym}" + "\n")
                f.write(r"\usepackage{xcolor}" + "\n")
                f.write(r"\renewcommand{\sfdefault}{cmr}" + "\n")
                f.write("\n")
                f.write(r"\usepackage{ascmac}" + "\n")
                f.write(r"\usepackage{fancybox}" + "\n")
                f.write(r"\usepackage{tcolorbox}" + "\n")
                f.write(r"\usepackage{ulem}" + "\n")
                f.write(r"\usepackage{pxrubrica}" + "\n")
                f.write(r"\usepackage{seqsplit}" + "\n")
                f.write(r"\usepackage{enumitem}" + "\n")
                f.write("\n")
                f.write(r"\usepackage{amssymb,amsmath,amsfonts}" + "\n")
                f.write(r"\usepackage{physics}" + "\n")
                f.write(r"\usepackage[version=4]{mhchem}" + "\n")
                f.write(r"\usepackage{bm}" + "\n")
                f.write("\n")
                f.write(r"\usepackage{graphicx}" + "\n")
                f.write(r"\usepackage{float}" + "\n")
                f.write(r"\usepackage{booktabs}" + "\n")
                f.write(r"\usepackage{longtable}" + "\n")
                f.write(r"\usepackage{tabularx}" + "\n")
                f.write(r"\usepackage{xltabular}" + "\n")
                f.write(r"\usepackage{colortbl}" + "\n")
                f.write(r"\usepackage{multicol}" + "\n")
                f.write(r"\usepackage{multirow}" + "\n")
                f.write(r"\usepackage{dcolumn}" + "\n")
                f.write(r"\usepackage{caption}" + "\n")
                f.write(r"\captionsetup[figure]{labelformat=empty}" + "\n")
                f.write(r"\captionsetup[table]{labelformat=empty}" + "\n")
                f.write("\n")
                f.write(r"\usepackage[colorlinks=true,linkcolor=blue,urlcolor=cyan,citecolor=red]{hyperref}" + "\n")
                f.write("\n")
                f.write(r"\usepackage{titlesec}" + "\n")
                f.write(r"\titleformat{\section}{\LARGE\bfseries}{\normalfont\thesection}{1em}{}" + "\n")
                f.write(r"\titleformat{\subsection}{\Large\bfseries}{\normalfont\thesection}{1em}{}" + "\n")
                f.write("\n")

                f.write(r"\title{Roughness Analysis Results}" + "\n")
                f.write(r"\author{}" + "\n")
                f.write(r"\vspace{-5zw}" + "\n")
                f.write(r"\date{\today}" + "\n")
                f.write("\n")
                f.write(r"\begin{document}" + "\n")
                f.write(r"\maketitle" + "\n")
                text_summary = desc_summary.replace("   <p>", "").replace("</p>", "").replace("      <br>　", r"\par"+"\n").strip()
                f.write(text_summary + "\n")
                f.write(r"\clearpage" + "\n")
                f.write("\n")

                f.write(r"\section*{Thickness Results}" + "\n")
                f.write(r"\subsection*{Correlation Plots}" + "\n")
                text_ea = desc_ea.replace("         <p>", "").replace("</p>", "").replace("            <br>　", r"\par"+"\n").strip()
                f.write(text_ea + "\n")
                if os.path.exists(os.path.join(self.dir_output, os.path.basename(board_img))):
                    f.write(r"\begin{figure}[H]" + "\n")
                    f.write(r"    \centering" + "\n")
                    f.write(r"    \includegraphics[width=0.5\textwidth]{" + os.path.basename(board_img) + "}" + "\n")
                    f.write(r"    \label{fig:BoardSchematic}" + "\n")
                    f.write(r"\end{figure}" + "\n")
                ea_test = defaultdict(list)
                for group_key, file_path in ea_correlations.items():
                    match = re.search(r'(test\d+)', group_key)
                    if match:
                        ea_test[match.group(1)].append(file_path)
                sorted_test_nums = sorted(ea_test.keys(), key=lambda x: int(x.replace('test','')))
                
                for test_num in sorted_test_nums:
                    f.write(r"\begin{figure}[H]" + "\n")
                    f.write(r"    \centering" + "\n")
                    group_data = ea_test[test_num]
                    for i, plot_file in enumerate(group_data):
                        f.write(r"    \begin{minipage}[t]{0.48\textwidth}" + "\n")
                        f.write(r"        \centering" + "\n")
                        f.write(r"        \includegraphics[width=\linewidth]{" + plot_file + "}" + "\n")
                        f.write(r"    \end{minipage}" + "\n")
                        if i % 2 == 0 and i < len(group_data) - 1:
                            f.write(r"    \hfill" + "\n")
                    f.write(r"    \label{fig:" + test_num + "_pos_correlations}" + "\n")
                    f.write(r"\end{figure}" + "\n")
                f.write(r"\clearpage" + "\n")
                f.write("\n")

                f.write(r"\subsection*{Average Thickness Table}" + "\n")
                text_ave_table = desc_ave_table.replace("   <p>", "").replace("</p>", "").replace("            <br>　", r"\par"+"\n").strip()
                f.write(text_ave_table + "\n")
                f.write(r"\begin{table}[H]" + "\n")
                f.write(r"    \resizebox{\textwidth}{!}{" + "\n")
                f.write(r"        \begin{tabular}{cc" + r"D{+}{\,\pm\,}{6.6}" * len(ea_nums) + r"}" + "\n") 
                f.write(r"            \toprule" + "\n")
                f.write(r"            \multicolumn{2}{c}{} & " + r" & ".join([f"\\multicolumn{{1}}{{c}}{{{x}}}" for x in ea_nums]) + r" \\" + "\n")
                f.write(r"            \midrule" + "\n")
                for i, test_name in enumerate(test_nums):
                    if i > 0:
                        f.write(r"            \addlinespace" + "\n")
                    f.write(r"            \multirow{2}{*}{" + test_name + r"} & Rise ($\rm\AA$) & ")
                    rise_val = []
                    for ea_name in ea_nums:
                        key = f"{test_name}_{ea_name}"
                        summary_data = ea_summary.get(key)
                        if summary_data:
                            rise_val.append(f"{summary_data['ea_ave_rise']:.2f}+{summary_data['ea_err_rise']:.2f}")
                        else:
                            rise_val.append(r'\multicolumn{1}{c}{---}')
                    f.write(' & '.join(rise_val) + r" \\" + "\n")
                    f.write(r"        & Fall ($\rm\AA$) & ")
                    fall_val = []
                    for ea_name in ea_nums:
                        key = f"{test_name}_{ea_name}"
                        summary_data = ea_summary.get(key)
                        if summary_data:
                            fall_val.append(f"{summary_data['ea_ave_fall']:.2f}+{summary_data['ea_err_fall']:.2f}")
                        else:
                            fall_val.append(r'\multicolumn{1}{c}{---}')
                    f.write(' & '.join(fall_val) + r" \\" + "\n")
                f.write(r"            \bottomrule" + "\n")
                f.write(r"        \end{tabular}" + "\n")
                f.write(r"    }" + "\n")
                f.write(r"\end{table}" + "\n")
                f.write(r"\clearpage" + "\n")
                f.write("\n")

                f.write(r"\section*{Thickness Results Table}" + "\n")
                text_table = desc_table.replace("   <p>", "").replace("</p>", "").replace("         <br>　", r"\par"+"\n").strip()
                f.write(text_table + "\n")
                f.write(r"\begin{xltabular}{\textwidth}{>{\raggedright\arraybackslash}X" + r"D{+}{\,\pm\,}{6.6}"*2 + r"}" + "\n")
                f.write(r"    \toprule" + "\n")
                f.write(r"    \multicolumn{1}{c}{File Name} & \multicolumn{1}{c}{Rise Height ($\rm\AA$)} & \multicolumn{1}{c}{Fall Height ($\rm\AA$)} \\" + "\n")
                f.write(r"    \midrule" + "\n")
                for r in results:
                    file_name         = r['InputFile'].replace('_', r'\_')
                    file_label        = os.path.splitext(r['InputFile'])[0].replace('_', '')
                    rise_height       = r['RiseHeight'] if r['RiseHeight'] is not None else '---'
                    fall_height       = r['FallHeight'] if r['FallHeight'] is not None else '---'
                    rise_height_error = r['RiseHeightError'] if r['RiseHeightError'] is not None else '---'
                    fall_height_error = r['FallHeightError'] if r['FallHeightError'] is not None else '---'
                    f.write(f"            \\hyperref[fig:{file_label}]{{{file_name}}} & {rise_height:.2f}+{rise_height_error:.2f} & {fall_height:.2f}+{fall_height_error:.2f} \\\\" + "\n")
                f.write(r"    \bottomrule" + "\n")
                f.write(r"\end{xltabular}" + "\n")
                f.write(r"\clearpage" + "\n")
                f.write("\n")

                f.write(r"\section*{Correlation Results}" + "\n")
                text_test = desc_test.replace("      <p>", "").replace("</p>", "").replace("         <br>　", r"\par"+"\n").strip()
                f.write(text_test + "\n")
                for group_name, correlation_file in test_correlations.items():
                    f.write(r"\begin{figure}[H]" + "\n")
                    f.write(r"    \centering" + "\n")
                    f.write(r"    \includegraphics[width=0.8\textwidth]{" + correlation_file + "}" + "\n")
                    f.write(r"    \label{fig:Correlation" + group_name.replace(' ', '') + "}" + "\n")
                    f.write(r"\end{figure}" + "\n")
                f.write(r"\clearpage" + "\n")
                f.write("\n")

                f.write(r"\section*{Each Plot}" + "\n")
                text_plot = desc_plot.replace("      <p>", "").replace("</p>", "").replace("         <br>　", r"\par"+"\n").replace("_", r"\_").strip()
                f.write(text_plot + "\n")
                for r in results:
                    plot_png_path = r['EachPlotResult']
                    if not plot_png_path: continue
                    f.write(r"\begin{figure}[H]" + "\n")
                    f.write(r"    \centering" + "\n")
                    f.write(r"    \includegraphics[width=\textwidth]{" + plot_png_path + "}" + "\n")
                    f.write(r"    \label{fig:" + os.path.splitext(r['InputFile'])[0].replace('_', '') + "}" + "\n")
                    f.write(r"\end{figure}" + "\n")
                f.write("\n")
                f.write(r"\end{document}" + "\n")
            print(f"\n LaTeX saved to: {output_path}")
            print(" Please run the LaTeX file separately.")

        else:
            print(f"\n \033[31mError: Output format '{self.output_format}' is incorrect or not specified. If specified, only 'html' or 'pdf' is allowed.\033[0m")
