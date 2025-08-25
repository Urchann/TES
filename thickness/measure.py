#!/usr/bin/env python

# -*- coding: utf-8 -*-

r"""
 __  __                                                       
|  \/  |  ___   __ _  ___  _   _  _ __   ___     _ __   _   _ 
| |\/| | / _ \ / _` |/ __|| | | || '__| / _ \   | '_ \ | | | |
| |  | ||  __/| (_| |\__ \| |_| || |   |  __/ _ | |_) || |_| |
|_|  |_| \___| \__,_||___/ \__,_||_|    \___|(_)| .__/  \__, |
                                                |_|     |___/ 
"""

# メモ
__version__ = '1.0' # 2025.07.04
__version__ = '1.1' # 2025.07.05 --> plot_correlation関数の追加と、細かな修正を行いました。
__version__ = '1.2' # 2025.07.07 --> plot_correlation関数において各基板でのプロットを可能にし、オプション引数で指定ファイルを除外した処理の構築をしました。
__version__ = '1.3' # 2025.08.17 --> plot_correlation関数をplot_test_correlation関数に命名変更、cal_ea_thickness関数とplot_ea_correlation関数の追加など、Ea番号ごとの処理を可能にし、その他大幅な修正を加えました。
__author__  = 'Ryota Fukuda' # mailto:25la018c@rikkyo.ac.jp

# インポート
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'
import pandas as pd
import argparse
import os
import subprocess
import shutil
from scipy.ndimage import uniform_filter1d
import re

# 解析内容の説明
desc_summary = """   <p>　理研CRの表面粗さ計（触針型）による測定データから，試料の厚みを計算する。
      <br>　粗さデータセットの波形は矩形波上になっており，本解析では，その矩形波の凸の上昇部分を「立ち上がり」，下降部分を「立ち下がり」と呼ぶ。"""\
      """この立ち上がりと立ち下がりの波高値をそれぞれ計算することにより，試料の厚みを評価する。
      <br>　波高値の計算方法は，まず素データをスムージングした後，隣り合うデータ点の差を取り，その最大値と最小値をそれぞれ立ち上がり点および立ち下がり点と決定する。"""\
      """次に，立ち上がり点からデータ始点，立ち下がり点からデータ終点までについてスキャンし，それぞれデータ点が増加する点を立ち上がり始点と立ち下がり始点とする。"""\
      """さらに，立ち上がり点から立ち下がり点，立ち下がり点から立ち上がり点までについてスキャンし，それぞれデータ点が減少する点を立ち上がり終点と立ち下がり終点とする。"""\
      """最後に，立ち上がりおよび立ち下がりの始点と終点の差を計算することで，試料の厚みをそれぞれ評価する。
   </p>"""

# ----------------------------------------------------
# 立ち上がりおよび立ち下がりを検出する関数
# ----------------------------------------------------
def get_indices(data, window_size=50):
    """
    理研CRの表面粗さ計（触針型）による測定データから、試料の厚みを計算する。
    粗さデータセットから、立ち上がりの始めと終わり、および立ち下がりの始めと終わりの検出をして、それぞれの高さ（厚み）を計算する。
    
    Args:
        data (pd.Series)  : 得られた測定データ
        window_size (int) : 平滑化（移動平均）を行う際の平均をとる範囲の長さ

    Returns:
        detects (dict): 
            'smoothed_data' : 測定データを平滑化したもの
            'pre_rise_idx'  : 立ち上がり始めのインデックス
            'post_rise_idx' : 立ち上がり終わりのインデックス
            'pre_fall_idx'  : 立ち下がり始めのインデックス
            'post_fall_idx' : 立ち下がり終わりのインデックス
            'pre_rise_val'  : 立ち上がり始めの値
            'post_rise_val' : 立ち上がり終わりの値
            'pre_fall_val'  : 立ち下がり始めの値
            'post_fall_val' : 立ち下がり終わりの値
            'rise_height'   : 立ち上がりの高さ
            'fall_height'   : 立ち下がりの高さ
    """
    # 返り値の用意
    detects = {
        'smoothed_data':None,
        'pre_rise_idx': None,
        'post_rise_idx': None,
        'pre_fall_idx': None,
        'post_fall_idx': None,
        'pre_rise_val': None,
        'post_rise_val': None,
        'pre_fall_val': None,
        'post_fall_val': None,
        'rise_height': None,
        'fall_height': None
    }

    if len(data) < 3:
        print("Data is too short.")
        return detects
    
    # 測定データのスムージング
    detects['smoothed_data'] = uniform_filter1d(data, size=window_size)

    # 立ち上がりおよび立ち下がりの探索
    diff = np.diff(detects['smoothed_data'])
    mid_rise_idx = np.argmax(diff)
    mid_fall_idx = np.argmin(diff)

    # 立ち上がり始めの探索
    for i in range(mid_rise_idx, 0, -1):
        if detects['smoothed_data'][i] <= detects['smoothed_data'][i-1]:
            detects['pre_rise_idx'] = i
            detects['pre_rise_val'] = detects['smoothed_data'][i]
            break

    # 立ち上がり終わりの探索
    for i in range(mid_rise_idx, mid_fall_idx):
        if detects['smoothed_data'][i] >= detects['smoothed_data'][i+1]:
            detects['post_rise_idx'] = i
            detects['post_rise_val'] = detects['smoothed_data'][i]
            break

    # 立ち下がり始めの探索
    for i in range(mid_fall_idx-1, mid_rise_idx, -1):
        if detects['smoothed_data'][i] >= detects['smoothed_data'][i-1]:
            detects['pre_fall_idx'] = i
            detects['pre_fall_val'] = detects['smoothed_data'][i]
            break

    # 立ち下がり終わりの探索
    for i in range(mid_fall_idx+1, len(data)-1):
        if detects['smoothed_data'][i] <= detects['smoothed_data'][i+1]:
            detects['post_fall_idx'] = i
            detects['post_fall_val'] = detects['smoothed_data'][i]
            break
    
    # 立ち上がりおよび立ち下がりの高さ（厚み）を計算
    detects['rise_height'] = detects['post_rise_val'] - detects['pre_rise_val']
    detects['fall_height'] = detects['pre_fall_val'] - detects['post_fall_val']

    return detects

# ----------------------------------------------------
# 基板の位置ごとの厚みを計算する関数
# ----------------------------------------------------
def cal_ea_thickness(results_group):
    """
    各データセット (基板) の位置ごとの立ち上がりおよび立ち下がりの高さ (厚み) の結果について、平均値の計算をする。
    膜厚値の誤差はMeasurement Index (0:pre_rise_idx) の不偏標準偏差を用い、厚みの平均値は各データセットの厚みの不偏標準偏差を用いる。

    Args:
        results_ea_group (list) : 任意のEa番号 (例: Ea1) に属するresultsリスト

    Returns:
        ave_rise (float)  : 立ち上がり (厚み) の平均値
        ave_fall (float)  : 立ち下がり (厚み) の平均値
        err_rise (float)  : 立ち上がり (厚み) の平均値の誤差
        err_fall (float)  : 立ち下がり (厚み) の平均値の誤差
    """
    # 値と誤差の取得・計算
    point = []
    rise  = []
    fall  = []
    yerr  = []
    for r in results_group:
        match = re.search(r'(LtoR|RtoL)_(Ea\d+)_(\d+)', r['InputFile'])
        if match:
            point.append(f"{match.group(2)}_{match.group(1)}_{match.group(3)}")
        else:
            point.append('null')
        rise.append(r['RiseHeight'])
        fall.append(r['FallHeight'])
        
        # 誤差の計算
        yerr.append(np.std(r['DataList'][0:r['PreRiseIndex']], ddof=1) / np.sqrt(r['PreRiseIndex']))
        
    # 立ち上がりおよび立ち下がり（厚み）の平均値とその誤差を計算
    ave_rise = np.mean(rise)
    ave_fall = np.mean(fall)
    if len(rise) < 2 or len(fall) < 2:
        err_rise = yerr[0]
        err_fall = yerr[0]
    else:
        err_rise = np.std(rise, ddof=1) / np.sqrt(len(rise))
        err_fall = np.std(fall, ddof=1) / np.sqrt(len(fall))
    return ave_rise, ave_fall, err_rise, err_fall

# ----------------------------------------------------
# 各データセットの結果をプロットする関数
# ----------------------------------------------------
def plot_results(file, data, smoothed_data, output_dir, window_size, pre_rise_idx, post_rise_idx, pre_fall_idx, post_fall_idx, pre_rise_val, post_rise_val, pre_fall_val, post_fall_val, rise_height, fall_height):
    """
    測定データをプロットし、立ち上がりの始めと終わり、および立ち下がりの始めと終わりの検出結果、さらに立ち上がりと立ち下がりそれぞれについて高さ（厚み）の結果を表示する。
    各データセットに対するプロット結果は、 (データセットが格納されているディレクトリの名前)_results という名前のディレクトリに、PNGファイルで格納する。

    Args:
        file (str)            : 得られた測定データのパス
        data (pd.Series)      : 得られた測定データ
        smoothed_data (np.array) : 測定データを平滑化したもの
        output_dir (str)      : プロット結果の出力先ディレクトリの名前
        window_size (int)     : 平滑化（移動平均）を行う際の平均をとる範囲の長さ
        pre_rise_idx (int)    : 立ち上がり始めのインデックス
        post_rise_idx (int)   : 立ち上がり終わりのインデックス
        pre_fall_idx (int)    : 立ち下がり始めのインデックス
        post_fall_idx (int)   : 立ち下がり終わりのインデックス
        pre_rise_val (float)  : 立ち上がり始めの値
        post_rise_val (float) : 立ち上がり終わりの値
        pre_fall_val (float)  : 立ち下がり始めの値
        post_fall_val (float) : 立ち下がり終わりの値
        rise_height (float)   : 立ち上がりの高さ（厚み）
        fall_height (float)   : 立ち下がりの高さ（厚み）

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        description (str) : 図と計算方法の説明
    """
    description = rf"""      <p>　全測定データを各々プロットする。
         <br>　横軸は測定位置のインデックスで，縦軸が膜厚値である。青色のラインは素データ，橙色のラインはuniform_filter1dを用いてスムージングしたデータを表す。"""\
         """また，緑色のマーカーは立ち上がり始めの点，赤色のマーカーは立ち上がり終わりの点，シアン色のマーカーは立ち下がりは始めの点，マゼンタ色のマーカーは立ち下がり終わりの点を表す。
      </p>"""
    
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.basename(file)
    output_file = os.path.splitext(filename)[0] + ".png" 
    output_path = os.path.join(output_dir, output_file)

    # 測定データとそれを平滑化したものをプロット
    plt.figure(figsize=(12, 6))
    plt.plot(data.index, data, label='Original Data', alpha=0.7, linewidth=0.8)
    plt.plot(data.index, smoothed_data, label=f'Smoothed Data (Window: {window_size})', color='orange', linewidth=1.5)

    # 立ち上がり始めと終わり、および立ち下がりの始めと終わりの点をプロット
    if pre_rise_idx is not None:
        plt.plot(pre_rise_idx, smoothed_data[pre_rise_idx], 'o', color='green', markersize=8, label=f'Pre-Rise ({pre_rise_val:.2f} Å)')
    if post_rise_idx is not None:
        plt.plot(post_rise_idx, smoothed_data[post_rise_idx], 'o', color='red', markersize=8, label=f'Post-Rise ({post_rise_val:.2f} Å)')
    if pre_fall_idx is not None:
        plt.plot(pre_fall_idx, smoothed_data[pre_fall_idx], 'o', color='cyan', markersize=8, label=f'Pre-Fall ({pre_fall_val:.2f} Å)')
    if post_fall_idx is not None:
        plt.plot(post_fall_idx, smoothed_data[post_fall_idx], 'o', color='magenta', markersize=8, label=f'Post-Fall ({post_fall_val:.2f} Å)')

    plt.axhline(pre_rise_val, color='black', linestyle='--')
    plt.axhline(post_rise_val, color='black', linestyle='--')
    plt.axhline(pre_fall_val, color='gray', linestyle='-.')
    plt.axhline(post_fall_val, color='gray', linestyle='-.')

    # 厚みの計算結果を表示
    plt.text(0.05, 0.9,
             f'Rise Thickness: {rise_height:.2f} (Å)',
             transform=plt.gca().transAxes,
             ha='left', va='top',
             color='darkgreen',
             fontsize=12,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='lightgray')
    )
    plt.text(0.05, 0.8,
             f'Fall Thickness: {fall_height:.2f} (Å)',
             transform=plt.gca().transAxes,
             ha='left', va='top',
             color='darkred',
             fontsize=12,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='lightgray')
    )

    # プロットの設定
    plt.title(f'Roughness Analysis for {filename}')
    plt.xlabel('Measurement Index')
    plt.ylabel('Roughness (Å)')
    plt.legend(loc=1)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Plot saved to: {output_path}")
    return output_file, description

# ----------------------------------------------------
# 基板ごとの相関をプロットする関数
# ----------------------------------------------------
def plot_test_correlation(output_dir, results_group, group_name):
    """
    各データセット (基板) の立ち上がりおよび立ち下がりの高さ (厚み) の結果について、平均値の計算および相関関係のプロットをする。
    計算方法は、cal_ea_thickness 関数と同様。
    プロット結果は (データセットが格納されているディレクトリの名前)_results という名前のディレクトリに、PNGファイルで格納する。

    Args:
        output_dir (str)     : プロット結果の出力先ディレクトリの名前
        results_group (list) : 任意のグループ (例: test1) に属するresultsリスト
        group_name (str)     : 相関をプロットする基板の名前 (例: 'test1')

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        ave_rise (float)  : 立ち上がり (厚み) の平均値
        ave_fall (float)  : 立ち下がり (厚み) の平均値
        err_rise (float)  : 立ち上がり (厚み) の平均値の誤差
        err_fall (float)  : 立ち下がり (厚み) の平均値の誤差
        description (str) : 図と計算方法の説明
    """
    description = rf"""      <p>　基板ごとに，立ち上がりと立ち下がりによる膜厚値の全計算結果をまとめてプロットする。
         <br>　横軸は測定点で，縦軸が膜厚値である。また，赤色の誤差付きマーカーは立ち上がりの膜厚値，青色の誤差付きマーカーは立ち下がりの膜厚値を表す。
         <br>　誤差は，各データについて，Measurement Indexにおける0から立ち上がり始めまでのデータ点から，不偏標準偏差で計算する。
      </p>"""
    
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    output_file = f"thickness_correlation_{group_name}.png"
    output_path = os.path.join(output_dir, output_file)

    # 値の取得
    point = []
    rise  = []
    fall  = []
    yerr  = []
    for r in results_group:
        match = re.search(r'(LtoR|RtoL)_(Ea\d+)_(\d+)', r['InputFile'])
        if match:
            point.append(f"{match.group(1)}_{match.group(2)}_{match.group(3)}")
        else:
            point.append('null')
        rise.append(r['RiseHeight'])
        fall.append(r['FallHeight'])
        
        # 誤差の計算
        yerr.append(np.std(r['DataList'][0:r['PreRiseIndex']], ddof=1) / np.sqrt(r['PreRiseIndex']))
    
    # 立ち上がりおよび立ち下がり（厚み）の平均値とその誤差を計算
    ave_rise = sum(rise) / len(rise)
    ave_fall = sum(fall) / len(fall)
    err_rise = np.std(rise, ddof=1) / np.sqrt(len(rise))
    err_fall = np.std(fall, ddof=1) / np.sqrt(len(fall))

    # プロット
    plt.figure(figsize=(8, 6))
    plt.errorbar(point, rise, label='Rise Thickness (Å)', yerr=yerr, color='red',  capsize=3, markersize=5, marker='o', markeredgecolor='red',  markerfacecolor='red',  ls='None')
    plt.errorbar(point, fall, label='Fall Thickness (Å)', yerr=yerr, color='blue', capsize=3, markersize=5, marker='o', markeredgecolor='blue', markerfacecolor='blue', ls='None')
    plt.title(f'Thickness correlation for {group_name} measurement points')
    plt.xlabel('Measurement Points')
    plt.ylabel('Thickness (Å)')
    plt.xticks(rotation=90)
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    return output_file, ave_rise, ave_fall, err_rise, err_fall, description

# ----------------------------------------------------
# 基板の位置ごとの相関をプロットする関数
# ----------------------------------------------------
def plot_ea_correlation(output_dir, results_group, ea_list, group_name, group_suffix):
    """
    各データセット (基板) の立ち上がりおよび立ち下がりの高さ (厚み) の結果について、相関関係のプロットをする。
    計算方法は、cal_ea_thickness関数と同様。相関プロットは、基板の斜めの位置関係にあるEa番号の相関を見る。
    プロット結果は (データセットが格納されているディレクトリの名前)_results という名前のディレクトリに、PNGファイルで格納する。

    Args:
        output_dir (str)     : プロット結果の出力先ディレクトリの名前
        results_group (list) : 任意のグループ (例: test1) に属するresultsリスト
        ea_list (list)       : プロットに含めるEa番号のリスト (例: ['Ea1', 'Ea3', 'Ea5'])
        group_name (str)     : プロットする基板位置の名前 (例: 'test1_Ea1')
        group_suffix (str)   : プロットに含めるEa番号グループ (例: 'Ea1-3-5')

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        description (str) : 図と計算方法の説明
    """
    description = rf"""         <p>　基板ごとに，測定位置による相関をプロットする。以下，上の図は基板の測定位置（Ea番号）の詳細。"""\
            """相関を見る位置は，(Ea1, Ea3, Ea5) と (Ea2, Ea3, Ea4) の2つとし，以下の下の図にまとめる。
            <br>　誤差は，各データについて，Measurement Indexにおける0から立ち上がり始めまでのデータ点から，不偏標準偏差で計算する。
         </p>"""
    
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    output_file = f"thickness_correlation_{group_name}_{group_suffix}.png"
    output_path = os.path.join(output_dir, output_file)

    # プロットに含めるEaを抽出
    sample = []
    for r in results_group:
        match = re.search(r'(Ea\d+)', r['InputFile'])
        if match and match.group(1) in ea_list:
            sample.append(r)

    # データをEa番号でソート
    def sort_key(data_dict):
        match = re.search(r'Ea(\d+)', data_dict['InputFile'])
        return int(match.group(1)) if match else float('inf')
    filter = sorted(sample, key=sort_key)

    # 値の取得
    point = []
    rise  = []
    fall  = []
    yerr  = []
    for f in filter:
        match = re.search(r'(LtoR|RtoL)_(Ea\d+)_(\d+)', f['InputFile'])
        if match:
            point.append(f"{match.group(2)}")
        else:
            point.append('null')
        rise.append(f['RiseHeight'])
        fall.append(f['FallHeight'])
        
        # 誤差の計算
        yerr.append(np.std(f['DataList'][0:f['PreRiseIndex']], ddof=1) / np.sqrt(f['PreRiseIndex']))
    
    # Ea番号ごとに平均値を計算
    ea_data = {}
    for f in filter:
        match = re.search(r'(Ea\d+)', f['InputFile'])
        if match:
            ea_num = match.group(1)
            if ea_num not in ea_data:
                ea_data[ea_num] = []
            ea_data[ea_num].append(f)

    ave_point = sorted(ea_data.keys(), key=lambda x: int(re.search(r'\d+', x).group()))
    ave_rise  = []
    ave_fall  = []
    ave_err_rise = []
    ave_err_fall = []
    for ea_num in ave_point:
        ea_results = ea_data[ea_num]
        ea_ave_rise, ea_ave_fall, ea_err_rise, ea_err_fall = cal_ea_thickness(ea_results)
        ave_rise.append(ea_ave_rise)
        ave_fall.append(ea_ave_fall)
        ave_err_rise.append(ea_err_rise)
        ave_err_fall.append(ea_err_fall)

    # プロット
    plt.figure(figsize=(8, 6))
    plt.errorbar(point, rise, label='Rise Thickness (Å)', yerr=yerr, color='red',  capsize=3, markersize=5, marker='o', markeredgecolor='red',  markerfacecolor='red',  ls='None')
    plt.errorbar(point, fall, label='Fall Thickness (Å)', yerr=yerr, color='blue', capsize=3, markersize=5, marker='o', markeredgecolor='blue', markerfacecolor='blue', ls='None')
    plt.errorbar(ave_point, ave_rise, label='Rise Average Thickness (Å)', yerr=ave_err_rise, color='green', capsize=3, marker='o', markersize=5, linestyle='-', linewidth=2)
    plt.errorbar(ave_point, ave_fall, label='Fall Average Thickness (Å)', yerr=ave_err_fall, color='orange', capsize=3, marker='s', markersize=5, linestyle='-', linewidth=2)
    plt.title(f'Thickness correlation for {group_name} measurement points')
    plt.xlabel('Measurement Points')
    plt.ylabel('Thickness (Å)')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    return output_file, description

# ----------------------------------------------------
# 解析結果をまとめる関数
# ----------------------------------------------------
def export_results(results, test_correlations, ea_correlations, ea_summary, output_dir, output_format, description):
    """
    解析結果をHTML形式ないしはLaTeX (PDF) 形式でまとめる。
    解析ファイル名、立ち上がり始め値、立ち上がり終わり値、立ち上がり高さ（厚み）、立ち下がり始め値、立ち下がり終わり値、立ち下がり高さ（厚み）を順に表にまとめて、
    プロット結果を順に表示する。

    Args:
        results (dict)           : 解析ファイル名、データセットリスト、立ち上がり始めインデックス、立ち上がり終わりインデックス、立ち下がり始めインデックス、立ち下がり終わりインデックス、
                                   立ち上がり始め値、立ち上がり終わり値、立ち下がり始め値、立ち下がり終わり値、立ち上がり高さ（厚み）、立ち下がり高さ（厚み）
        test_correlations (dict) : 各基板の相関プロット結果を格納したディクショナリ
        ea_correlations (dict)   : 各基板位置の相関プロット結果を格納したディクショナリ
        ea_summary (dict)        : 各基板位置の厚みの計算結果を格納したディクショナリ
        output_dir (str)         : 解析結果の出力先ディレクトリ
        output_format (str)      : 出力フォーマット（HTML形式またはLaTeX (PDF) 形式）
    """
    # 基板概略図をカレントディレクトリから出力先ディレクトリにコピー
    img_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BoardSchematic.png')
    if os.path.exists(img_path):
        shutil.copy(img_path, output_dir)
        print(f"\nBoardSchematic.png has been copied to {output_dir}.")
    else:
        print(f"\nBoardSchematic.png not found. Create HTML without board schematic.")
    
    # HTML出力
    if output_format == 'html':
        output_path = os.path.join(output_dir, "thickness_results.html")
        output_dir  = os.path.basename(os.path.dirname(output_path))
        with open(output_path, 'w') as f:
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
            f.write(f"   <p><strong>Author:</strong> {__author__}</p>\n")
            f.write(description[0] + "\n")
            
            f.write("   <div class='plot-section'>\n")
            f.write("      <h2>Thickness Results</h2>\n")
            f.write("      <div class='summary-box'>\n")
            f.write("         <h3>Correlation Plot</h3>\n")
            f.write(description[3] + "\n")
            f.write("         <div class='plot-container'>\n")
            f.write(f"            <a href='./BoardSchematic.png' target='_blank'><img src='./BoardSchematic.png' alt='Measurement position on the board'></a>\n")
            f.write("         </div>\n")
            f.write("         <div class='correlation-grid'>\n")
            def sort_ea_key(key):
                match = re.search(r'test(\d+)', key)
                number = int(match.group(1)) if match else float('inf')
                group  = {'Ea1-3-5': 0, 'Ea2-3-4': 1}
                suffix = key.split('_')[-1]
                order  = group.get(suffix, float('inf'))
                return (number, order)
            sorted_ea_keys = sorted(ea_correlations.keys(), key=sort_ea_key)
            for group_key in sorted_ea_keys:
                plot_file = ea_correlations[group_key]
                if plot_file:
                    f.write("            <div class='plot-container'>\n")
                    f.write(f"               <a href='./{plot_file}' target='_blank'><img src='./{plot_file}' alt='Correlation of thickness for {group_key}'></a>\n")
                    f.write(f"               <p>{group_key.replace('_', ' ')}</p>\n")
                    f.write("            </div>\n")
            f.write("         </div>\n")

            f.write("         <h3>Average Thickness Table</h3>\n")
            f.write("         <p>　各基板の位置ごとに，立ち上がり時間と立ち下がり時間による膜厚値の結果（平均値）を表にまとめる。<br>　誤差は，平均値の不偏標準偏差で計算する。</p>\n")
            f.write("         <div class='summary-table-container'>\n")
            f.write("            <table class='summary-table' style='width: 100%; border-collapse: collapse; margin-top: 20px;'>\n")
            f.write("               <thead>\n")
            f.write("                  <tr>\n")
            f.write("                     <th rowspan='2'></th>\n")
            test_nums = sorted(list(set([k.split('_')[0] for k in ea_summary.keys() if 'test' in k])))
            ea_nums   = sorted(list(set([k.split('_')[1] for k in ea_summary.keys() if 'Ea' in k])), key=lambda x: int(re.search(r'Ea(\d+)', x).group(1)))
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
            f.write("                     <th style='background-color: #e6e6e6; border: 1px solid #b0b0b0;'>Rise</th>\n")
            for test_name in test_nums:
                for ea_name in ea_nums:
                    key = f"{test_name}_{ea_name}"
                    summary_data = ea_summary.get(key)
                    rise_val = f"{summary_data['ea_ave_rise']:.2f} \\(\\pm\\) {summary_data['ea_err_rise']:.2f}"
                    f.write(f"                     <td>{rise_val}</td>\n")
            f.write("                  </tr>\n")
            f.write("                  <tr>\n")
            f.write("                     <th style='background-color: #e6e6e6; border: 1px solid #b0b0b0;'>Fall</th>\n")
            for test_name in test_nums:
                for ea_name in ea_nums:
                    key = f"{test_name}_{ea_name}"
                    summary_data = ea_summary.get(key)
                    fall_val = f"{summary_data['ea_ave_fall']:.2f} \\(\\pm\\) {summary_data['ea_err_fall']:.2f}"
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
            f.write("      <p>　全測定データにおける膜厚値の表にまとめる。<br>　1列目から順に，データファイル名，立ち上がり始めの点，立ち上がり終わりの点，立ち上がりの膜厚値，立ち下がり始めの点，立ち下がり終わりの点，立ち下がりの膜厚値である。</p>\n")
            f.write("      <table>\n")
            f.write("         <thead>\n")
            f.write("            <tr><th>File Name</th><th>Pre-Rise Value (Å)</th><th>Post-Rise Value (Å)</th><th>Rise Height (Å)</th><th>Pre-Fall Value (Å)</th><th>Post-Fall Value (Å)</th><th>Fall Height (Å)</th></tr>\n")
            f.write("         </thead>\n")
            f.write("         <tbody>\n")
            for r in results:
                pre_rise_val  = r['PreRiseValue'] if r['PreRiseValue'] is not None else '---'
                post_rise_val = r['PostRiseValue'] if r['PostRiseValue'] is not None else '---'
                pre_fall_val  = r['PreFallValue'] if r['PreFallValue'] is not None else '---'
                post_fall_val = r['PostFallValue'] if r['PostFallValue'] is not None else '---'
                rise_height   = r['RiseHeight'] if r['RiseHeight'] is not None else '---'
                fall_height   = r['FallHeight'] if r['FallHeight'] is not None else '---'

                f.write(f"            <tr><td>{r['InputFile']}</td>"
                        f"<td>{pre_rise_val:.2f}</td>"
                        f"<td>{post_rise_val:.2f}</td>"
                        f"<td>{rise_height:.2f}</td>"
                        f"<td>{pre_fall_val:.2f}</td>"
                        f"<td>{post_fall_val:.2f}</td>"
                        f"<td>{fall_height:.2f}</td>"
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
        print(f"\033[32m\nHTML saved to: {output_path}\033[0m")

    # LaTeX (PDF) 出力
    elif output_format == 'pdf':
        output_path = os.path.join(output_dir, "thickness_results.tex")
        with open(output_path, 'w') as f:
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
            f.write(r"理研CRの表面粗さ計（触針型）による測定データから，試料の厚みを計算する。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"粗さデータセットの波形は矩形波上になっており，本解析では，その矩形波の凸の上昇部分を「立ち上がり」，下降部分を「立ち下がり」と呼ぶ。この立ち上がりと立ち下がりの波高値をそれぞれ計算することにより，試料の厚みを評価する。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"波高値の計算方法は，まず素データをスムージングした後，隣り合うデータ点の差を取り，その最大値と最小値をそれぞれ立ち上がり点および立ち下がり点と決定する。次に，立ち上がり点からデータ始点，立ち下がり点からデータ終点までについてスキャンし，それぞれデータ点が増加する点を立ち上がり始点と立ち下がり始点とする。さらに，立ち上がり点から立ち下がり点，立ち下がり点から立ち上がり点までについてスキャンし，それぞれデータ点が減少する点を立ち上がり終点と立ち下がり終点とする。最後に，立ち上がりおよび立ち下がりの始点と終点の差を計算することで，試料の厚みをそれぞれ評価する。" + "\n")
            f.write(r"\clearpage" + "\n")
            f.write("\n")

            f.write(r"\section*{Thickness Results}" + "\n")
            f.write(r"\subsection*{Correlation Plots}" + "\n")
            f.write(r"基板ごとに，測定位置による相関をプロットする。以下，上の図は基板の測定位置（Ea番号）の詳細。 相関を見る位置は，(Ea1, Ea3, Ea5)と(Ea2, Ea3, Ea4)の2つとし，以下の下の図にまとめる。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"誤差は，各データについて，Measurement Indexにおける0から立ち上がり始めまでのデータ点から，不偏標準偏差で計算する。" + "\n")
            f.write(r"\begin{figure}[H]" + "\n")
            f.write(r"    \centering" + "\n")
            f.write(r"    \includegraphics[width=0.5\textwidth]{BoardSchematic.png}" + "\n")
            f.write(r"    \label{fig:BoardSchematic}" + "\n")
            f.write(r"\end{figure}" + "\n")
            ea_test = {}
            for group_key, file_path in ea_correlations.items():
                match = re.search(r'(test\d+)', group_key)
                if match:
                    test_num = match.group(1)
                    if test_num not in ea_test:
                        ea_test[test_num] = []
                    ea_test[test_num].append((group_key, file_path))
            sorted_test_nums = sorted(ea_test.keys(), key=lambda x: int(x.replace('test','')))
            
            for test_num in sorted_test_nums:
                f.write(r"\begin{figure}[H]" + "\n")
                f.write(r"    \centering" + "\n")
                group_data = ea_test[test_num]
                sorted_group_data = sorted(group_data, key=lambda x: 0 if 'Ea1-3-5' in x[0] else 1)
                for i, (group_key, plot_file) in enumerate(sorted_group_data):
                    f.write(r"    \begin{minipage}[t]{0.48\textwidth}" + "\n")
                    f.write(r"        \centering" + "\n")
                    f.write(r"        \includegraphics[width=\linewidth]{" + plot_file + "}" + "\n")
                    f.write(r"    \end{minipage}" + "\n")
                    if i % 2 == 0 and i < len(sorted_group_data) - 1:
                        f.write(r"    \hfill" + "\n")
                f.write(r"    \label{fig:" + test_num + "_ea_correlations}" + "\n")
                f.write(r"\end{figure}" + "\n")
            f.write(r"\clearpage" + "\n")
            f.write("\n")

            f.write(r"\subsection*{Average Thickness Table}" + "\n")
            f.write(r"各基板の位置ごとに，立ち上がり時間と立ち下がり時間による膜厚値の結果（平均値）を表にまとめる。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"誤差は，平均値の不偏標準偏差で計算する。" + "\n")
            f.write(r"\begin{table}[H]" + "\n")
            f.write(r"    \centering" + "\n")
            f.write(r"    \begin{tabular}{ccD{+}{\,\pm\,}{6,6}D{+}{\,\pm\,}{6,6}D{+}{\,\pm\,}{6,6}D{+}{\,\pm\,}{6,6}D{+}{\,\pm\,}{6,6}}" + "\n")
            f.write(r"        \toprule" + "\n")
            test_nums = sorted(list(set([k.split('_')[0] for k in ea_summary.keys() if 'test' in k])))
            ea_nums   = sorted(list(set([k.split('_')[1] for k in ea_summary.keys() if 'Ea' in k])), key=lambda x: int(re.search(r'Ea(\d+)', x).group(1)))
            f.write(r"        & & " + r" & ".join([f"\\multicolumn{{1}}{{c}}{{{x}}}" for x in ea_nums]) + r" \\" + "\n")
            f.write(r"        \midrule" + "\n")
            for i, test_name in enumerate(test_nums):
                if i > 0:
                    f.write(r"        \addlinespace" + "\n")
                f.write(r"        \multirow{2}{*}{" + test_name + r"} & Rise & ")
                rise_val = []
                for ea_name in ea_nums:
                    key = f"{test_name}_{ea_name}"
                    summary_data = ea_summary.get(key)
                    rise_val.append(f"{summary_data['ea_ave_rise']:.2f}+{summary_data['ea_err_rise']:.2f}")
                f.write(' & '.join(rise_val) + r" \\" + "\n")
                f.write(r"        & Fall & ")
                fall_val = []
                for ea_name in ea_nums:
                    key = f"{test_name}_{ea_name}"
                    summary_data = ea_summary.get(key)
                    fall_val.append(f"{summary_data['ea_ave_fall']:.2f}+{summary_data['ea_err_fall']:.2f}")
                f.write(' & '.join(fall_val) + r" \\" + "\n")
            f.write(r"        \bottomrule" + "\n")
            f.write(r"    \end{tabular}" + "\n")
            f.write(r"\end{table}" + "\n")
            f.write(r"\clearpage" + "\n")
            f.write("\n")

            f.write(r"\section*{Thickness Results Table}" + "\n")
            f.write(r"全測定データにおける膜厚値の表にまとめる。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"1列目から順に，データファイル名，立ち上がり始めの点，立ち上がり終わりの点，立ち上がりの膜厚値，立ち下がり始めの点，立ち下がり終わりの点，立ち下がりの膜厚値である。" + "\n")
            f.write(r"\begin{table}[H]" + "\n")
            f.write(r"    \resizebox{\textwidth}{!}{" + "\n")
            f.write(r"        \begin{tabular}{lD{.}{.}{4.2}D{.}{.}{4.2}D{.}{.}{4.2}D{.}{.}{4.2}D{.}{.}{4.2}D{.}{.}{4.2}}" + "\n")
            f.write(r"            \toprule" + "\n")
            f.write(r"            \multicolumn{1}{c}{File Name} & \multicolumn{1}{c}{Pre-Rise Value ($\rm\AA$)} & \multicolumn{1}{c}{Post-Rise Value ($\rm\AA$)} & \multicolumn{1}{c}{Rise Height ($\rm\AA$)} & \multicolumn{1}{c}{Pre-Fall Value ($\rm\AA$)} & \multicolumn{1}{c}{Post-Fall Value ($\rm\AA$)} & \multicolumn{1}{c}{Fall Height ($\rm\AA$)} \\" + "\n")
            f.write(r"            \midrule" + "\n")
            for r in results:
                file = r['InputFile'].replace('_', r'\_')
                label = os.path.splitext(r['InputFile'])[0].replace('_', '')
                pre_rise_val = r['PreRiseValue'] if r['PreRiseValue'] is not None else '---'
                post_rise_val = r['PostRiseValue'] if r['PostRiseValue'] is not None else '---'
                pre_fall_val = r['PreFallValue'] if r['PreFallValue'] is not None else '---'
                post_fall_val = r['PostFallValue'] if r['PostFallValue'] is not None else '---'
                rise_height = r['RiseHeight'] if r['RiseHeight'] is not None else '---'
                fall_height = r['FallHeight'] if r['FallHeight'] is not None else '---'
                f.write(f"            \\hyperref[fig:{label}]{{{file}}} & {pre_rise_val:.2f} & {post_rise_val:.2f} & {rise_height:.2f} & {pre_fall_val:.2f} & {post_fall_val:.2f} & {fall_height:.2f} \\\\" + "\n")
            f.write(r"            \bottomrule" + "\n")
            f.write(r"        \end{tabular}" + "\n")
            f.write(r"    }" + "\n")
            f.write(r"\end{table}" + "\n")
            f.write(r"\clearpage" + "\n")
            f.write("\n")

            f.write(r"\section*{Correlation Results}" + "\n")
            f.write(r"基板ごとに，立ち上がりと立ち下がりによる膜厚値の全計算結果をまとめてプロットする。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"横軸は測定点で，縦軸が膜厚値である。また，赤色の誤差付きマーカーは立ち上がりの膜厚値，青色の誤差付きマーカーは立ち下がりの膜厚値を表す。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"誤差は，各データについて，Measurement Indexにおける0から立ち上がり始めまでのデータ点から，不偏標準偏差で計算する。" + "\n")
            for group_name, correlation_file in test_correlations.items():
                f.write(r"\begin{figure}[H]" + "\n")
                f.write(r"    \centering" + "\n")
                f.write(r"    \includegraphics[width=0.8\textwidth]{" + correlation_file + "}" + "\n")
                f.write(r"    \label{fig:Correlation" + group_name.replace(' ', '') + "}" + "\n")
                f.write(r"\end{figure}" + "\n")
            f.write(r"\clearpage" + "\n")
            f.write("\n")

            f.write(r"\section*{Each Plot}" + "\n")
            f.write(r"全測定データを各々プロットする。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"横軸は測定位置のインデックスで，縦軸が膜厚値である。青色のラインは素データ，橙色のラインはuniform\_filter1dを用いてスムージングしたデータを表す。また，緑色のマーカーは立ち上がり始めの点，赤色のマーカーは立ち上がり終わりの点，シアン色のマーカーは立ち下がりは始めの点，マゼンタ色のマーカーは立ち下がり終わりの点を表す。" + "\n")
            for r in results:
                plot_png_path = r['EachPlotResult']
                f.write(r"\begin{figure}[H]" + "\n")
                f.write(r"    \centering" + "\n")
                f.write(r"    \includegraphics[width=\textwidth]{" + plot_png_path + "}" + "\n")
                f.write(r"    \label{fig:" + os.path.splitext(r['InputFile'])[0].replace('_', '') + "}" + "\n")
                f.write(r"\end{figure}" + "\n")
            f.write("\n")
            f.write(r"\end{document}" + "\n")
        print(f"\033[32m\nLaTeX saved to: {output_path}\033[0m")
        print("Please run the LaTeX file separately.")

    else:
        print("\nNo data processed or results found.")


# ----------------------------------------------------
# メインプログラム
# ----------------------------------------------------
if __name__ == "__main__":
    # 出力先ディレクトリを "(解析に指定したファイルの親ディレクトリ名)_results" に設定する関数
    def get_default_output_dir(input_files):
        input_file_path = os.path.abspath(input_files[0])
        parent_dir_path = os.path.dirname(input_file_path)
        return os.path.basename(parent_dir_path) + "_results"
    
    # コマンドライン引数の設定
    parser = argparse.ArgumentParser(description="Analyze roughness from a TXT file measured with a RIKEN CR surface roughness tester (stylus type).", 
                                     epilog="TXT files to be analyzed should be placed in a directory at the same level.\n" \
                                            "BoardSchematic.png should also be placed in a directory at the same level.\n" \
                                            "Run Example: python measure.py LOR3A_TEST/*.txt --output_format html",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_files', nargs='+', 
                        help="Input TXT files (e.g., 'LOR3A_TEST/*.txt').\n" \
                             "It is recommended to place all TXT files to be analyzed in the same directory and select them at once using a wildcard.")
    parser.add_argument('--exclude_files', nargs='*', default=[],
                        help="Space-separated list of TXT files to exclude from analysis (e.g., 'LOR3A_test1_250526_test.txt').")
    parser.add_argument('--output_dir', default=None,
                        help="Directory to save plots and results table. (Default: 'parent_directory_name'_results)")
    parser.add_argument('--data_column', choices=['Raw', 'RawLevel', 'Normal', 'Rough', 'Wavi'], default='Normal',
                        help="Column name to use for measurement data. (Default: 'Normal')")
    parser.add_argument('--output_format', choices=['pdf', 'html'], default='pdf',
                        help="Output format for the results table: 'pdf' (LaTeX) or 'html'. (Default: 'pdf')")
    parser.add_argument('--window_size', type=int, default=50,
                        help="Window size for data smoothing (moving average). (Default: 50)")
    args = parser.parse_args()

    # 出力先ディレクトリを作成
    if args.output_dir is None:
        args.output_dir = get_default_output_dir(args.input_files)
    os.makedirs(args.output_dir, exist_ok=True)

    # 解析ファイルを抽出
    excluded_basenames = {os.path.basename(f) for f in args.exclude_files}
    files = sorted([f for f in args.input_files if os.path.basename(f) not in excluded_basenames])

    if not files:
        print("\nNo files to process after exclusion. Exiting.")
        exit()

    results = []
    for file in files:
        print(f"\nProcessing file: {file}")
        try:
            # TXTファイルから指定の列（オプション引数 --data_column）を抽出
            with open(file, 'r') as f:
                header = 0
                for line in f:
                    header += 1
                    if args.data_column in line:
                        break
            df = pd.read_csv(file, sep=r'\s+', skiprows=header-1, skipinitialspace=True)

            if args.data_column not in df.columns:
                print(f"\nError: Column '{args.data_column}' not found in {file}. Available columns: {df.columns.tolist()}")
                continue
            data = df[args.data_column]

            # 立ち上がりおよび立ち下がりのインデックスを特定
            detects = get_indices(data,
                                  window_size=args.window_size
            )
            
            # 各データセットの結果をプロット
            plots, desc_plot = plot_results(file,
                                            data,
                                            detects['smoothed_data'], 
                                            args.output_dir,
                                            args.window_size,
                                            detects['pre_rise_idx'],
                                            detects['post_rise_idx'],
                                            detects['pre_fall_idx'],
                                            detects['post_fall_idx'],
                                            detects['pre_rise_val'],
                                            detects['post_rise_val'],
                                            detects['pre_fall_val'],
                                            detects['post_fall_val'],
                                            detects['rise_height'],
                                            detects['fall_height']
            )

            # 結果を格納
            results.append({
                'InputFile': os.path.basename(file),
                'DataList': data,
                'PreRiseIndex': detects['pre_rise_idx'],
                'PostRiseIndex': detects['post_rise_idx'],
                'PreFallIndex': detects['pre_fall_idx'],
                'PostFallIndex': detects['post_fall_idx'],
                'PreRiseValue': detects['pre_rise_val'],
                'PostRiseValue': detects['post_rise_val'],
                'PreFallValue': detects['pre_fall_val'],
                'PostFallValue': detects['post_fall_val'],
                'RiseHeight': detects['rise_height'],
                'FallHeight': detects['fall_height'],
                'EachPlotResult': plots
            })

            print(f"  Detected Pre-Rise      : {detects['pre_rise_idx']} ({detects['pre_rise_val']:.2f} Å)")
            print(f"  Detected Post-Rise     : {detects['post_rise_idx']} ({detects['post_rise_val']:.2f} Å)")
            print(f"  Detected Pre-Fall      : {detects['pre_fall_idx']} ({detects['pre_fall_val']:.2f} Å)")
            print(f"  Detected Post-Fall     : {detects['post_fall_idx']} ({detects['post_fall_val']:.2f} Å)")
            print(f"  \033[30;47mCalculated Rise Height : {detects['rise_height']:.2f} Å\033[0m")
            print(f"  \033[30;47mCalculated Fall Height : {detects['fall_height']:.2f} Å\033[0m")
            print(f"Finished processing: {file}")

        except Exception as e:
            print(f"Failed to process {file}: {e}")
            import traceback
            traceback.print_exc()

    # 結果を基板ごとにグループ化
    results_group = {}
    for r in results:
        match = re.search(r'test(\d+)', r['InputFile'])
        if match:
            group_key = f"test{match.group(1)}"
            if group_key not in results_group:
                results_group[group_key] = []
            results_group[group_key].append(r)
        else:
            # test(数字)形式でないファイルも考慮する場合
            if 'other' not in results_group:
                results_group['other'] = []
            results_group['other'].append(r)

    # 各基板およびその基板位置ごとに相関プロットを生成
    test_correlations = {}
    ea_summary        = {}
    ea_correlations   = {}
    ea_group_filters  = {'Ea1-3-5': ['Ea1', 'Ea3', 'Ea5'],
                        'Ea2-3-4': ['Ea2', 'Ea3', 'Ea4']
    }
    if results_group:
        print("\n\n--- Thickness Analysis Summary ---")
        for group_name, results_list in results_group.items():
            if not results_list:
                continue

            # test番号ごとのプロットを生成
            test_correlation_file, test_ave_rise, test_ave_fall, test_err_rise, test_err_fall, desc_test = plot_test_correlation(args.output_dir, results_list, group_name)
            test_correlations[group_name] = test_correlation_file
            print(f"\n[Summary for {group_name}]")
            print(f"  - Overall Average Rise Thickness: {test_ave_rise:.2f} ± {test_err_rise} (Å)")
            print(f"  - Overall Average Fall Thickness: {test_ave_fall:.2f} ± {test_err_fall} (Å)")

            # Ea番号ごとのデータをグループ化
            ea_groups = {}
            for r in results_list:
                match = re.search(r'Ea(\d+)', r['InputFile'])
                if match:
                    ea_key = f"{group_name}_Ea{match.group(1)}"
                    if ea_key not in ea_groups:
                        ea_groups[ea_key] = []
                    ea_groups[ea_key].append(r)
            
            # Ea番号ごとのプロットを生成
            for group_suffix, ea_list in ea_group_filters.items():
                ea_correlation_file, desc_ea = plot_ea_correlation(args.output_dir, results_list, ea_list, group_name, group_suffix)
                ea_correlations[f"{group_name}_{group_suffix}"] = ea_correlation_file

            # Ea番号ごとに厚み結果を計算
            print("\n  [Ea-specific Summary]")
            for ea_key, ea_results in ea_groups.items():
                ea_ave_rise, ea_ave_fall, ea_err_rise, ea_err_fall = cal_ea_thickness(ea_results)
                ea_summary[ea_key] = {
                    'ea_ave_rise': ea_ave_rise,
                    'ea_err_rise': ea_err_rise,
                    'ea_ave_fall': ea_ave_fall,
                    'ea_err_fall': ea_err_fall
                }
                print(f"    - {ea_key} Average Rise Thickness: {ea_ave_rise:.2f} ± {ea_err_rise:.2f} (Å)")
                print(f"    - {ea_key} Average Fall Thickness: {ea_ave_fall:.2f} ± {ea_err_fall:.2f} (Å)")
        
        # 解析方法や図の説明についての情報をまとめる
        description  = [desc_summary, desc_plot, desc_test, desc_ea]

        # 結果をオプション引数 --output_format で書き出し
        export_results(results, test_correlations, ea_correlations, ea_summary, args.output_dir, args.output_format, description)
        print("\nCompleted!\n")
    else:
        print("\nNo data processed or results found.\n")