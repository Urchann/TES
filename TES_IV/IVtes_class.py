#!/usr/bin/env python

# -*- coding: utf-8 -*-

r"""
 ___ __     __ _                               
|_ _|\ \   / /| |_   ___  ___     _ __   _   _ 
 | |  \ \ / / | __| / _ \/ __|   | '_ \ | | | |
 | |   \ V /  | |_ |  __/\__ \ _ | |_) || |_| |
|___|   \_/    \__| \___||___/(_)| .__/  \__, |
                                 |_|     |___/ 
"""

__version__ = '1.0'
__author__  = 'Ryota Fukuda' # mailto:25la018c@rikkyo.ac.jp
__credits__ = 'Tasuku Hayashi' # mailto:tasuku.hayashi@riken.jp
__url__     = 'https://colab.research.google.com/drive/1rRyx1eIyx2i56KETObivz26Km9cdnt5h?authuser=1'

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'
import lmfit as lmf
import argparse
import os
import subprocess
import shutil
import re
from jinja2 import Environment, FileSystemLoader


# ----------------------------------------------------
# TESの物理量を計算する関数群
# ----------------------------------------------------
def cal_Ites(Vout, Xi):
    """"
    TESに流れる電流（Ites）を計算する。

    Args:
        Vout (float) : 測定データ（出力電圧）[V]
        Xi (float)   : 電流電圧変換係数
    
    Returns:
        Ites (float) : TESに流れる電流（Ites）
    """
    Ites = Vout / Xi
    return Ites

def cal_Rtes(Itb, Ites, Rsh):
    """"
    TESの抵抗（Rtes）を計算する。

    Args:
        Itb (float)  : バイアス電流 [A]
        Ites (float) : TESに流れる電流 [A]
        Rsh (float)  : シャント抵抗 [Ω]
    
    Returns:
        Rtes (float) : TESの抵抗 [Ω]
    """
    Rtes = (Itb - Ites) * Rsh / Ites
    return Rtes

def cal_Vtes(Rtes, Ites):
    """"
    TES両端にかかる電圧（Vtes）を計算する。

    Args:
        Rtes (float) : TESの抵抗 [Ω]
        Ites (float) : TESに流れる電流 [A]
    
    Returns:
        Vtes (float) : TES両端にかかる電圧 [V]
    """
    Vtes = Rtes * Ites
    return Vtes

def pfunc(Tb, Tc, G0, n):
    """"
    TESのJoule発熱（Ptes）のフィッティング関数を作る。

    Args:
        Tb (float) : 熱浴温度 [K]
        Tc (float) : TESの温度 [K]
        G0 (float) : 熱浴温度が 1K のときの熱伝導度 [W/K]
        n (float)  : べき定数
    
    Returns:
        Pfunc (float) : TESのJoule発熱（フィッティング用）[W]
    """
    Pfunc = (G0 / n) * (Tc**n - Tb**n)
    return Pfunc

def cal_Ttes(Tb, G0, n, Ptes):
    """"
    TESの温度（Ttes）を計算する。

    Args:
        Tb (float)   : 熱浴温度 [K]
        G0 (float)   : 熱浴温度が 1K のときの熱伝導度 [W/K]
        n (float)    : べき定数
        Ptes (float) : TESのJoule発熱 [W]
    
    Returns:
        Ttes (float) : TESの温度 [K]
    """
    Ttes = ((Tb**n) + ((n * Ptes) / G0))**(1/n)
    return Ttes

def cal_Gtes(G0, n, Ttes):
    """"
    TESの熱伝導度（Gtes）を計算する。

    Args:
        G0 (float)   : 熱浴温度が 1K のときの熱伝導度 [W/K]
        n (float)    : べき定数
        Ttes (float) : TESの温度 [K]
    
    Returns:
        Gtes (float) : TESの熱伝導度 [W/K]
    """
    Gtes = G0 * (Ttes**(n-1))
    return Gtes

def cal_alpha(T, R):
    """"
    TESの感度（alpha）を計算する。

    Args:
        T (float) : 温度 [K]
        R (float) : 抵抗 [Ω]
    
    Returns:
        alpha (float) : 感度
    """
    dR = np.diff(R)
    dT = np.diff(T)
    T_ave = (T[1:] + T[:-1])/2
    R_ave = (R[1:] + R[:-1])/2
    alpha = (T_ave * dR) / (R_ave * dT)
    return alpha

def get_Tbath(files):
    """"
    熱浴温度のリストを作成する。

    Args:
        fles (str) : 得られた測定データのパス
    
    Returns:
        Tbath (list) : 全熱浴温度（Tb [K]）のリスト
    """
    Tbath = []
    for tb in files:
        # 各データセットのパスを取得
        tb_ = os.path.basename(tb)
        num = 0
        for i, name in enumerate(tb_.split('_')):
            # パス名をアンダーバーで区切り、末尾に 'mK' がある場所で熱浴温度を判定
            if 'mK' in name:
                num = i
                break
        
        # 熱浴温度の値部分のみを抽出（単位を mK から K に変換）
        tb_key = tb_.split('_')[num]
        tbath = float(tb_key[:-2]) / 1e3
        Tbath.append(tbath)
    Tbath = np.asarray(Tbath)
    print(f"Tbath values detected: {Tbath}")
    return Tbath


# ----------------------------------------------------
# フィッティング関数
# ----------------------------------------------------
def fit_pfunc(output_dir, Tbath, Rtes, Ptes, Tcmin, Tcmax, Tcfit, G0min, G0max, G0fit, nmin, nmax, nfit):
    """"
    TESのJoule発熱と熱浴温度の関係をフィッティングする。
    フィッティングから、熱浴温度が 1K のときの熱伝導度 G0 , べき定数 n , TESの温度 Tc の各パラメータを得る。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (lise)     : 全熱浴温度（Tb [K]）のリスト
        Itb (float)      : バイアス電流 [A]
        Ttes (float)     : TESの温度 [K]
        Rtes(float)      : TESの抵抗 [Ω]
        Tb (int)         : 感度を計算する対象の熱浴温度 [K]

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        Pb (list)         : 各データにおける熱浴温度の代表値のリスト
        Pberr (list)      : 熱浴温度の代表値 Pb の誤差（標準偏差）
        Tc (float)        : TESの温度 Tc の最適値 [K]
        G0 (float)        : 熱浴温度が 1K のときの熱伝導度 G0 の最適値 [W/K]
        n (float)         : べき定数 n の最適値
        Ttes (float)      : TESの温度 [K]
    """
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.abspath(__file__)
    output_file = os.path.dirname(filename) + "_fitting.png" 
    output_path = os.path.join(output_dir, output_file)

    # 熱浴温度の代表値および誤差を取得
    Pb    = []
    Pberr = []
    for Tb in Tbath:
        rtes_   = Rtes[Tb][Rtes[Tb]==Rtes[Tb]]
        rtesmax = rtes_.max()
        rtesmin = rtes_.min()
        Rn      = rtesmax - rtesmin
        rtes    = Rtes[Tb]/Rn
        Prep    = (0.2<rtes) & (rtes<0.3)
        Pb.append(np.average(Ptes[Tb][Prep]))
        Pberr.append(np.std(Ptes[Tb][Prep], ddof=1))
    Pb    = np.asarray(Pb)
    Pberr = np.asarray(Pberr)

    # フィッティング
    model = lmf.Model(pfunc)
    model.make_params(verbose=True)
    model.set_param_hint('Tc', min=Tcmin, max=Tcmax)
    model.set_param_hint('G0', min=G0min, max=G0max)
    model.set_param_hint('n' , min=nmin , max=nmax)
    result = model.fit(Pb, Tb=Tb, Tc=Tcfit, G0=G0fit, n=nfit)
    print(result.fit_report())
    Tc = result.best_values['Tc']
    G0 = result.best_values['G0']
    n  = result.best_values['n']

    # χ2の計算
    chi2_min = sum(((Pb - pfunc(Tb=Tb, Tc=Tc, G0=G0, n=n))/Pberr)**2)
    print(chi2_min)

    # フィッティング結果をプロット
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.errorbar(Tbath*1e3, Pb*1e12, yerr=Pberr*1e12, color='black', capsize=4, ms=5, marker='o', mec='black', mfc='black', ls='None')
    Tb_min = Tbath.min() - 0.010
    Tb_max = Tbath.max() + 0.010
    x = np.linspace(Tb_min, Tb_max, 100)
    ax.plot(x*1e3, pfunc(x, Tc, G0, n)*1e12, 'r--')
    ax.set_title(r"$P_{\rm b}\ vs.\ T_{\rm bath}$")
    ax.set_xlabel(r"$T_{\rm bath}\ ({\rm mK})$", fontsize=14)
    ax.set_ylabel(r"$P_{\rm b}\ ({\rm pW})$", fontsize=14)
    ax.legend(['data', r"$T_{\rm c}=$"+rf"${Tc:.3}$"+r"${\rm K},\ G_0=$"+rf"${G0*1e9:.3}$"+r"${\rm nW/K},\ n=$"+rf"${n:.2f}$"], fontsize=14, frameon=False)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save Fitting to: {output_path}")
    return output_file, Pb, Pberr, chi2_min, Tc, G0, n


# ----------------------------------------------------
# TESパラメータの特性をプロットする関数群
# ----------------------------------------------------
def plot_IVproperty(output_dir, Tbath, Itb, Vout, Ites, Vtes):
    """"
    I-V特性グラフを作成する。
    出力電圧 vs. バイアス電流 & TESの両端にかかる電圧 vs. TESに流れる電流の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (lise)     : 全熱浴温度（Tb [K]）のリスト
        Itb (float)      : バイアス電流 [A]
        Vout (float)     : 出力電圧 [V]
        Ites (float)     : TESに流れる電流 [A]
        Vtes (float)     : TESの両端にかかる電圧 [V]
    
    Returns:
        output_file (str) : プロット結果の出力ファイル名
    """
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.abspath(__file__)
    output_file = os.path.dirname(filename) + "_IVproperty.png" 
    output_path = os.path.join(output_dir, output_file)

    _, ax = plt.subplots(1, 2, figsize=(12, 4))
    # 出力電圧 vs. バイアス電流
    for Tb in Tbath:
        ax[0].plot(Itb[Tb]*1e6, Vout[Tb], '.', label=r"$T_{\rm bath}=$"+str(int(Tb*1e3))+r"${\rm mK}$")
    ax[0].set_title(r"$V_{\rm out}\ vs.\ I_{\rm TES_{bias}}$")
    ax[0].set_xlabel(r"$I_{\rm TES_{bias}}\ ({\rm\mu A})$", fontsize=14)
    ax[0].set_ylabel(r"$V_{\rm out}\ ({\rm V})$", fontsize=14)
    ax[0].legend(fontsize=12,frameon=False, ncols=2)
    ax[0].grid(ls='--')

    # TESの両端にかかる電圧 vs. TESに流れる電流
    for Tb in Tbath:
        ax[1].plot(Vtes[Tb]*1e6, Ites[Tb]*1e6, '.', label=r"$T_{\rm bath}=$"+str(int(Tb*1e3))+r"${\rm mK}$")
    ax[1].set_title(r"$I_{\rm TES}\ vs.\ V_{\rm TES}$")
    ax[1].set_xlabel(r"$V_{\rm TES}\ ({\rm\mu V})$", fontsize=14)
    ax[1].set_ylabel(r"$I_{\rm TES}\ ({\rm\mu A})$", fontsize=14)
    ax[1].legend(fontsize=12,frameon=False, ncols=2)
    ax[1].grid(ls='--')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save I-V property to: {output_path}")
    return output_file

def plot_PRproperty(output_dir, Tbath, Rtes, Ptes, xmin, xmax, ymin, ymax):
    """"
    P-R特性グラフを作成する。
    TESのJoule発熱 vs. TESの抵抗の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (lise)     : 全熱浴温度（Tb [K]）のリスト
        Rtes (float)     : TESのJoule発熱 [W]
        Ptes (float)     : TESの抵抗 [Ω]
        xmin (float)     : 描画するx軸の最小範囲
        xmax (float)     : 描画するx軸の最大範囲
        ymin (float)     : 描画するy軸の最小範囲
        ymax (float)     : 描画するy軸の最大範囲
    
    Returns:
        output_file (str) : プロット結果の出力ファイル名
    """
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.abspath(__file__)
    output_file = os.path.dirname(filename) + "_PRproperty.png" 
    output_path = os.path.join(output_dir, output_file)

    # TESのJoule発熱 vs. TESの抵抗
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    for Tb in Tbath:
        rtes    = Rtes[Tb][Rtes[Tb]==Rtes[Tb]]
        rtesmax = rtes.max()
        rtesmin = rtes.min()
        Rn      = rtesmax - rtesmin
        ax.plot((Rtes[Tb]-rtesmin)/Rn, Ptes[Tb]*1e9, '.', label=r"$T_{\rm bath}=$"+str(int(Tb*1e3))+r"${\rm mK}$")
    # ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(r"$P_{\rm TES}\ vs.\ R_{\rm TES_{bias}}$")
    ax.set_xlabel(r"$R_{\rm TES}/R_N\ ({\rm\mu A})$", fontsize=14)
    ax.set_ylabel(r"$P_{\rm TES}\ ({\rm nW})$", fontsize=14)
    ax.legend(fontsize=12,frameon=False, ncols=2)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save R-P property to: {output_path}")
    return output_file

def plot_RTproperty(output_dir, Tbath, Rtes, Ptes, G0, n):
    """"
    R-T特性グラフを作成する。
    TESの抵抗 vs. TESの温度の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (lise)     : 全熱浴温度（Tb [K]）のリスト
        Rtes (float)     : TESのJoule発熱 [W]
        Ptes (float)     : TESの抵抗 [Ω]
        G0 (float)       : 熱浴温度が 1K のときの熱伝導度 [W/K]
        n (float)        : べき定数

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        Ttes (float)      : TESの温度 [K]
    """
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.abspath(__file__)
    output_file = os.path.dirname(filename) + "_RTproperty.png" 
    output_path = os.path.join(output_dir, output_file)

    # TESの抵抗 vs. TESの温度
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    for Tb in Tbath:
        Ttes[Tb] = cal_Ttes(Tb, G0, n, Ptes[Tb])
        ax.plot(Ttes[Tb]*1e3, Rtes[Tb], '.', label=r"$T_{\rm bath}=$"+str(int(Tb*1e3))+r"${\rm mK}$")
    ax.set_title(r"$R_{\rm TES}\ vs.\ T_{\rm TES_}$")
    ax.set_xlabel(r"$T_{\rm TES}\ ({\rm mK})$", fontsize=14)
    ax.set_ylabel(r"$R_{\rm TES}\ ({\rm\Omega})$", fontsize=14)
    ax.legend(fontsize=12,frameon=False, ncols=2)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save R-T property to: {output_path}")
    return output_file, Ttes

def plot_GTproperty(output_dir, Tbath, Ttes, G0, n):
    """"
    G-T特性グラフを作成する。
    TESの熱伝導度 vs. TESの温度の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (lise)     : 全熱浴温度（Tb [K]）のリスト
        Ttes (float)     : TESの温度 [K]
        G0 (float)       : 熱浴温度が 1K のときの熱伝導度 [W/K]
        n (float)        : べき定数

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        Gtes (float)      : TESの熱伝導度 [W/K]
    """
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.abspath(__file__)
    output_file = os.path.dirname(filename) + "_GTproperty.png" 
    output_path = os.path.join(output_dir, output_file)

    # TESの熱伝導度 vs. TESの温度
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    for Tb in Tbath:
        Gtes = cal_Gtes(G0, n, Ttes[Tb])
        ax.plot(Ttes[Tb]*1e3, Gtes[Tb], '.', label=r"$T_{\rm bath}=$"+str(int(Tb*1e3))+r"${\rm mK}$")
    ax.set_title(r"$G_{\rm TES}\ vs.\ T_{\rm TES_}$")
    ax.set_xlabel(r"$T_{\rm TES}\ ({\rm mK})$", fontsize=14)
    ax.set_ylabel(r"$G_{\rm TES}\ ({\rm W/K})$", fontsize=14)
    ax.legend(fontsize=12,frameon=False, ncols=2)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save G-T property to: {output_path}")
    return output_file, Gtes

def plot_alpha(output_dir, Itb, Ttes, Rtes, Tb, xmin=0., xmax=200., ymin=150., ymax=300.):
    """"
    TESの感度グラフを作成する。
    TESの感度 vs. バイアス電流の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Itb (float)      : バイアス電流 [A]
        Ttes (float)     : TESの温度 [K]
        Rtes(float)      : TESの抵抗 [Ω]
        Tb (int)         : 感度を計算する対象の熱浴温度 [K]
        xmin (float)     : 描画するx軸の最小範囲
        xmax (float)     : 描画するx軸の最大範囲
        ymin (float)     : 描画するy軸の最小範囲
        ymax (float)     : 描画するy軸の最大範囲

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        alpha (float)     : TESの感度
    """
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.abspath(__file__)
    output_file = os.path.dirname(filename) + "_alpha.png" 
    output_path = os.path.join(output_dir, output_file)

    # TESの感度 vs. バイアス電流
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    alpha = cal_alpha(Ttes[Tb], Rtes[Tb])
    Ibias = (Itb[Tb][1:] + Itb[Tb][:-1]) / 2
    ax.plot(Ibias*1e6, alpha)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title("Sensitivity "+r"$\alpha$")
    ax.set_xlabel(r"$I_{\rm bias}\ ({\rm\mu A})$", fontsize=14)
    ax.set_ylabel(r"$\alpha$", fontsize=14)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save TES Sensitivity to: {output_path}")
    return output_file, alpha


# ----------------------------------------------------
# コントアを計算する関数
# ----------------------------------------------------
def cal_contour(output_dir, Tbath, Pb, Perr, Tc, G0, n, chi2_min, xn1, xn2, xg0, yg0, ytc1, ytc2, cslevels):
    """"
    フィッティングパラメータのコントアを計算する。
    G0 vs. n  &  Tc vs. n  &  Tc vs. G0 のそれぞれについて計算する。

    Args:
        file (str)       : 得られた測定データのパス
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (lise)     : 全熱浴温度（Tb [K]）のリスト
        Pb (list)        : 各データにおける熱浴温度の代表値のリスト
        Pberr (list)     : Pb の誤差（標準偏差）
        Tc (float)       : TESの温度 Tc の最適値 [K]
        G0 (float)       : 熱浴温度が 1K のときの熱伝導度 G0 の最適値 [W/K]
        n (float)        : べき定数 n の最適値
        chi2_min (float) : フィッティングにおける最小χ二乗
        xn (list)        : べき定数 n のメッシュ配列
        xg0 (list)       : 熱浴温度が 1K のときの熱伝導度 G0 のメッシュ配列
        yg0 (list)       : 熱浴温度が 1K のときの熱伝導度 G0 のメッシュ配列
        ytc (list)       : TESの温度 Tc のメッシュ配列

    Returns:
        output_file (str) : プロット結果の出力ファイル名
    """
    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.abspath(__file__)
    output_file = os.path.dirname(filename) + "_contour.png" 
    output_path = os.path.join(output_dir, output_file)

    _, ax = plt.subplots(1, 3, figsize=(18, 4))

    # G0 と n のコントアを計算
    Xn, YG0 = np.meshgrid(xn1, yg0)
    for mn, nn in enumerate(xn1):
        chi2 = np.array([])
        for mg0, ng0 in enumerate(yg0):
            ng0 = ng0*1e-9
            chi2_cont = sum(((Pb - pfunc(Tb=Tbath, Tc=Tc, G0=ng0, n=nn))/Perr)**2)
            print(f'\r{mn+1}/{len(Xn[0])},{mg0+1}/{len(YG0.T[0])}',end='',flush=True)
            #print(f'\r{mn+1}/{len(Xn[0])}',end='',flush=True)
            chi2 = np.hstack((chi2, chi2_cont))

        if mn == 0:
            chi2_array = chi2
        else:
            chi2_array = np.vstack((chi2_array, chi2))

    chi2_array -= chi2_min

    # プロット
    cont0 = ax[0].contour(Xn, YG0, chi2_array.T, levels=cslevels)
    cont0.clabel(fmt='%1.1f', fontsize=14)
    ax[0].plot(n, G0*1e9, 'r*', ms=10)
    ax[0].set_title(r"$G_0$"+" vs. "+r"$n$"+" Contour")
    ax[0].set_xlabel(r'$n$', fontsize=14)
    ax[0].set_ylabel(r'$G_0\ ({\rm nW/K})$', fontsize=14)

    # Tc と n のコントアを計算
    Xn, YTc = np.meshgrid(xn2, ytc1)
    for mn, nn in enumerate(xn2):
        chi2 = np.array([])
        for mtc, ntc in enumerate(ytc1):
            ntc = ntc*1e-3
            chi2_cont = sum(((Pb - pfunc(Tb=Tbath, Tc=ntc, Go=G0, n=nn))/Perr)**2)
            print(f'\r{mn+1}/{len(Xn[0])},{mtc+1}/{len(YTc.T[0])}',end='',flush=True)
            #print(f'\r{mn+1}/{len(X[0])}',end='',flush=True)
            chi2 = np.hstack((chi2, chi2_cont))

        if mn == 0:
            chi2_array = chi2
        else:
            chi2_array = np.vstack((chi2_array, chi2))

    chi2_array -= chi2_min

    # プロット
    cont1 = ax[1].contour(Xn, YTc, chi2_array.T, levels=cslevels)
    cont1.clabel(fmt='%1.1f', fontsize=14)
    ax[1].plot(n, Tc*1e3, 'r*', ms=10)
    ax[1].set_title(r"$T_{\rm c}$"+" vs. "+r"$n$"+" Contour")
    ax[1].set_xlabel(r'$n$', fontsize=14)
    ax[1].set_ylabel(r'$T_{\rm c}\ ({\rm mK})$', fontsize=14)

    # Tc と G0 のコントアを計算
    XG0, YTc = np.meshgrid(xg0, ytc2)
    for mg0, ng0 in enumerate(xg0):
        chi2 = np.array([])
        for mtc, ntc in enumerate(ytc2):
            ntc = ntc*1e-3
            chi2_cont = sum(((Pb - pfunc(Tb=Tbath, Tc=ntc, Go=ng0, n=n))/Perr)**2)
            print(f'\r{mg0+1}/{len(Xn[0])},{mtc+1}/{len(YTc.T[0])}',end='',flush=True)
            #print(f'\r{mg0+1}/{len(X[0])}',end='',flush=True)
            chi2 = np.hstack((chi2, chi2_cont))

        if mn == 0:
            chi2_array = chi2
        else:
            chi2_array = np.vstack((chi2_array, chi2))

    chi2_array -= chi2_min

    # プロット
    cont2 = ax[2].contour(XG0, YTc, chi2_array.T, levels=cslevels)
    cont2.clabel(fmt='%1.1f', fontsize=14)
    ax[2].plot(G0*1e9, Tc*1e3, 'r*', ms=10)
    ax[2].set_title(r"$T_{\rm c}$"+" vs. "+r"$G_0$"+" Contour")
    ax[2].set_xlabel(r'$G_0\ ({\rm nW/K})$', fontsize=14)
    ax[2].set_ylabel(r'$T_{\rm c}\ ({\rm mK})$', fontsize=14)

    # 全コントアの計算結果をプロットし、保存
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save Fitting Parameter Contour to: {output_path}")
    return output_file


# ----------------------------------------------------
# 解析結果をまとめる関数
# ----------------------------------------------------
def export_results(all_results, correlations, output_dir, output_format, info_dict):
    """
    解析結果をHTML形式ないしはLaTeX (PDF) 形式でまとめる。
    解析ファイル名、フィットパラメータ、プロット結果を順に表示する。
    """
    env = Environment(loader=FileSystemLoader('.'))

    # HTML出力
    if output_format == 'html':
        output_path = os.path.join(output_dir, "tes_analysis_report.html")
        html_template = env.get_template('html_report_template_tes.html')
        html_content = html_template.render(
            title=info_dict['title'],
            author=info_dict['author'],
            date=info_dict['date'],
            all_results=all_results, # 各ファイルの解析結果
            correlations=correlations, # 相関プロット
            plot_files_for_html=[os.path.basename(r['EachPlotResult']) for r in all_results if r['EachPlotResult']] + \
                                [os.path.basename(correlations[k]) for k in correlations if correlations[k]] + \
                                [os.path.basename(p) for p in info_dict['global_plot_files']] # 全てのプロットファイルをリストで渡す
        )
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"\n\033[32mHTML report saved to: {output_path}\033[0m")

    # LaTeX (PDF) 出力
    elif output_format == 'pdf':
        output_file_base = "tes_analysis_report"
        latex_path = os.path.join(output_dir, f"{output_file_base}.tex")
        pdf_path = os.path.join(output_dir, f"{output_file_base}.pdf")

        latex_template = env.get_template('latex_report_template_tes.tex')
        latex_content = latex_template.render(
            title=info_dict['title'],
            author=info_dict['author'],
            date=info_dict['date'],
            tc_best_global=f"{info_dict['Tc_best']:.3f}" if not np.isnan(info_dict['Tc_best']) else "N/A",
            go_best_global=f"{info_dict['Go_best']*1e9:.3f}" if not np.isnan(info_dict['Go_best']) else "N/A",
            n_best_global=f"{info_dict['n_best']:.2f}" if not np.isnan(info_dict['n_best']) else "N/A",
            chi2_min_global=f"{info_dict['chi2_min']:.2f}" if not np.isnan(info_dict['chi2_min']) else "N/A",
            all_results=all_results,
            correlations=correlations,
            global_plot_files=[os.path.basename(p) for p in info_dict['global_plot_files']]
        )

        with open(latex_path, 'w', encoding='utf-8') as f:
            f.write(latex_content)
        print(f"LaTeX file generated: {latex_path}")

        print("Attempting to compile PDF...")
        current_dir = os.getcwd()
        os.chdir(output_dir) # PDF生成のために出力ディレクトリに移動

        try:
            # pdflatexを2回実行して、相互参照を解決
            subprocess.run(['xelatex', f'{output_file_base}.tex', '-interaction=nonstopmode', '-output-directory=.'], check=True, capture_output=True, text=True)
            subprocess.run(['xelatex', f'{output_file_base}.tex', '-interaction=nonstopmode', '-output-directory=.'], check=True, capture_output=True, text=True)
            print(f"\033[32mPDF successfully generated: {pdf_path}\033[0m")

            # 一時ファイルのクリーンアップ
            temp_extensions = ['.aux', '.log', '.out', '.toc', '.lof', '.lot']
            for ext in temp_extensions:
                temp_file = f"{output_file_base}{ext}"
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                    # print(f"  Removed temporary file: {temp_file}")
            
        except subprocess.CalledProcessError as e:
            print(f"Error compiling LaTeX to PDF: {e}")
            print(f"STDOUT:\n{e.stdout}")
            print(f"STDERR:\n{e.stderr}")
            print("Please ensure xelatex is installed and in your PATH.")
        except FileNotFoundError:
            print("Error: 'xelatex' command not found.")
            print("Please ensure xelatex is installed and in your system's PATH.")
        finally:
            os.chdir(current_dir) # 元のディレクトリに戻る

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
    parser = argparse.ArgumentParser(description="Analyze TES IV curve data.", 
                                     epilog="TXT files to be analyzed should be placed in a directory at the same level.\n" \
                                            "Run Example: python IVtes.py 230918/*.txt --output_format pdf",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_files', nargs='+', 
                        help="Input TXT files (e.g., '230918/*.txt').\n" \
                             "It is recommended to place all TXT files to be analyzed in the same directory and select them at once using a wildcard.")
    parser.add_argument('--exclude_files', nargs='*', default=[],
                        help="Space-separated list of TXT files to exclude from analysis (e.g., 'IV_070mK_20230918.txt').")
    parser.add_argument('--output_dir', default=None,
                        help="Directory to save plots and results table. (Default: 'parent_directory_name'_results)")
    parser.add_argument('--output_format', choices=['pdf', 'html'], default='pdf',
                        help="Output format for the results table: 'pdf' (LaTeX) or 'html'. (Default: 'pdf')")
    parser.add_argument('--Tcmin', type=float, default=0.100,
                        help="Minimum TES temperature fitting. (Default: 0.100 K)")
    parser.add_argument('--Tcmax', type=float, default=0.320,
                        help="Maximum TES temperature fitting. (Default: 0.320 K)")
    parser.add_argument('--Tcfit', type=float, default=0.290,
                        help="Initial value for TES temperature fitting. (Default: 0.290 K)")
    parser.add_argument('--G0min', type=float, default=1.e-15,
                        help="Minimum thermal conductivity fitting at heat bath temperature 1K. (Default: 1 pW/K)")
    parser.add_argument('--G0max', type=float, default=1.e2,
                        help="Maximum thermal conductivity fitting at heat bath temperature 1K. (Default: 10 W/K)")
    parser.add_argument('--G0fit', type=float, default=1.e-8,
                        help="Initial value for thermal conductivity fitting at heat bath temperature 1K. (Default: 10 nW/K)")
    parser.add_argument('--nmin', type=float, default=3.,
                        help="Minimum power constant fitting. (Default: 3)")
    parser.add_argument('--nmax', type=float, default=4.,
                        help="Maximum power constant fitting. (Default: 4)")
    parser.add_argument('--nfit', type=float, default=3.5,
                        help="Initial value for power constant fitting. (Default: 3.5)")
    parser.add_argument('--PR_xmin', type=float, default=0.,
                        help="Minimum x-axis range for plotting P-R property. (Default: 0 µA)")
    parser.add_argument('--PR_xmax', type=float, default=300.,
                        help="Maximum x-axis range for plotting P-R property. (Default: 300 µA)")
    parser.add_argument('--PR_ymin', type=float, default=0.,
                        help="Minimum y-axis range for plotting P-R property. (Default: 0 nW/K)")
    parser.add_argument('--PR_ymax', type=float, default=0.1,
                        help="Maximum y-axis range for plotting P-R property. (Default: 0.1 nW/K)")
    parser.add_argument('--alpha_Tb', type=float, default=120.,
                        help="Tbath for which alpha vs. Ibias plot is generated. (Default: 120.0 mK)")
    parser.add_argument('--alpha_xmin', type=float, default=0.,
                        help="Minimum x-axis range for plotting alpha-Itb property. (Default: 0 µA)")
    parser.add_argument('--alpha_xmax', type=float, default=200.,
                        help="Maximum x-axis range for plotting alpha-Itb property. (Default: 200 µA)")
    parser.add_argument('--alpha_ymin', type=float, default=150.,
                        help="Minimum y-axis range for plotting alpha-Itb property. (Default: 150)")
    parser.add_argument('--alpha_ymax', type=float, default=300.,
                        help="Maximum y-axis range for plotting alpha-Itb property. (Default: 300)")
    parser.add_argument('--xn1', default=np.arange(2.5, 3.5, 0.01),
                        help="Scan range of x-axis (power constant) contour. (Default: np.arange(2.5, 3.5, 0.01))")
    parser.add_argument('--xn2', default=np.arange(2.95, 3.05, 0.01),
                        help="Scan range of x-axis (power constant) contour. (Default: np.arange(2.95, 3.05, 0.01))")
    parser.add_argument('--xg0', default=np.arange(20, 40, 1),
                        help="Scan range of x-axis (thermal conductivity at heat bath temperature 1K) contour. (Default: np.arange(20, 40, 1))")
    parser.add_argument('--yg0', default=np.arange(0, 60, 0.1),
                        help="Scan range of y-axis (thermal conductivity at heat bath temperature 1K) contour. (Default: np.arange(0, 60, 0.1))")
    parser.add_argument('--ytc1', default=np.arange(154, 156.5, 0.01),
                        help="Scan range of y-axis (TES temperature) contour. (Default: np.arange(154, 156.5, 0.01))")
    parser.add_argument('--ytc2', default=np.arange(153, 157, 0.1),
                        help="Scan range of y-axis (TES temperature) contour. (Default: np.arange(153, 157, 0.1))")
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

    # --- 定数 ---
    # コントア定数（自由度2のカイ二乗分布における68.27%, 90%, 99%誤差の累積確率）--> https://keisan.site/exec/system/1161228834
    CRange68 = 2.295815160785974337606
    CRange90 = 4.605170185988091368036
    CRange99 = 9.210340371976182736072
    cslevels = [CRange68, CRange90, CRange99]
    # SQUID定数（シャント抵抗 Rsh ，入力コイル相互インダクタンス Min ，フィードバックコイル相互インダクタンス Mfb ，フィードバック抵抗 Rfb ，電流電圧変換係数 Xi ）
    Rsh = 3.9e-3    # Ω
    # Min =  # H
    # Mfb =  # H
    Rfb = 100.e3    # Ω
    Xi  = Rfb # * Min / Mfb
    # 熱浴温度
    Tbath = get_Tbath(files) # K

    # TESの物理量（出力電圧 Vout , バイアス電流 Itb , TESに流れる電流 Ites , TESの抵抗 Rtes , TESの両端にかかる電圧 Vtes , TESのJoule発熱 Ptes , TESの温度 Ttes）を取得
    Vout = {}
    Itb  = {}
    Ites = {}
    Rtes = {}
    Vtes = {}
    Ptes = {}
    Ttes = {}
    for i, file in enumerate(files):
        Tb = Tbath[i]
        print(f"\nProcessing file: {os.path.basename(file)} (Tbath = {Tb} mK)")
        try:
            # TXTファイルから出力電圧 Vout とバイアス電流 Itb を抽出、およびその他のパラメータを計算
            data = np.loadtxt(file)
            num  = len(data[:,0])
            Vout[Tb] = data[:,1] - data[:,1].min()      # V
            Itb[Tb]  = data[:,0] / 1e6                  # A (μA to A)
            Ites[Tb] = cal_Ites(Vout[Tb], Xi)           # A
            Rtes[Tb] = cal_Rtes(Itb[Tb], Ites[Tb], Rsh) # Ω
            Vtes[Tb] = cal_Vtes(Rtes[Tb], Ites[Tb])     # V
            Ptes[Tb] = Ites[Tb] * Vtes[Tb]              # W

            if i == 0:
                # 超伝導での傾きを求めて、傾きが正になるように電流電圧変換係数 Xi を補正
                tilt = np.average(np.diff(Ites[Tb][num-10:num]))/np.average(np.diff(Vtes[Tb][num-10:num]))
                print(f'T    = {Tb}')
                print(f'Xi   = {Xi}')
                print(f'tilt = {tilt}')
                j = 0
                while tilt < 0:
                    j  += 1
                    Xi += 1
                    Ites[Tb] = cal_Ites(Vout[Tb], Xi)
                    Rtes[Tb] = cal_Rtes(Itb[Tb], Ites[Tb], Rsh)
                    Vtes[Tb] = cal_Vtes(Rtes[Tb], Ites[Tb])
                    tilt     = np.average(np.diff(Ites[Tb][num-10:num]))/np.average(np.diff(Vtes[Tb][num-10:num]))
                print(f'{j} modified Xi   = {Xi}')
                print(f'{j} modified tilt = {tilt}')
            
            # 最終的なTESパラメータを計算
            Ites[Tb] = cal_Ites(Vout[Tb], Xi)
            Rtes[Tb] = cal_Rtes(Itb[Tb], Ites[Tb], Rsh)
            Vtes[Tb] = cal_Vtes(Rtes[Tb], Ites[Tb])
            Ptes[Tb] = Ites[Tb] * Vtes[Tb]

        except Exception as e:
            print(f"Error processing file {file}: {e}")
            Vout[Tb] = np.array([])
            Itb[Tb]  = np.array([])
            Ites[Tb] = np.array([])
            Rtes[Tb] = np.array([])
            Vtes[Tb] = np.array([])
            Ptes[Tb] = np.array([])
            continue

    result_IVproperty = plot_IVproperty(args.output_dir,
                                        Tbath,
                                        Itb,
                                        Vout,
                                        Ites,
                                        Vtes
    )
    result_PRproperty = plot_PRproperty(args.output_dir,
                                        Tbath,
                                        Rtes,
                                        Ptes,
                                        args.PR_xmin,
                                        args.PR_xmax,
                                        args.PR_ymin,
                                        args.PR_ymax
    )
    result_fitting, Pb, Pberr, chi2_min, Tc, G0, n = fit_pfunc(args.output_dir,
                                                               Tbath,
                                                               Rtes,
                                                               Ptes,
                                                               args.Tcmin,
                                                               args.Tcmax,
                                                               args.Tcfit,
                                                               args.G0min,
                                                               args.G0max,
                                                               args.G0fit,
                                                               args.nmin,
                                                               args.nmax,
                                                               args.nfit
    )
    result_RTproperty, Ttes = plot_RTproperty(args.output_dir,
                                        Tbath,
                                        Rtes,
                                        Ptes,
                                        G0,
                                        n
    )
    result_GTproperty, Gtes = plot_GTproperty(args.output_dir,
                                              Tbath,
                                              Ttes,
                                              G0,
                                              n
    )
    result_alpha, alpha = plot_alpha(args.output_dir,
                                     Itb,
                                     Ttes,
                                     Rtes,
                                     args.alpha_Tb,
                                     args.alpha_xmin,
                                     args.alpha_xmax,
                                     args.alpha_ymin,
                                     args.alpha_ymax
    )
    result_contour = cal_contour(args.output_dir,
                                 Tbath,
                                 Pb,
                                 Pberr,
                                 Tc,
                                 G0,
                                 n,
                                 chi2_min,
                                 args.xn1,
                                 args.xn2,
                                 args.xg0,
                                 args.yg0,
                                 args.ytc1,
                                 args.ytc2,
                                 cslevels
    )

    results = []
    

    print("\nCompleted TES IV Curve Analysis!")