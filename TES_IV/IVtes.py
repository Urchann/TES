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

# メモ
__version__ = '1.0' # 2025.07.12
__version__ = '1.1' # 2025.08.03 --> export_results関数の調整を行いました。
__version__ = '1.2' # 2025.08.20 --> main関数とexport_results関数の調整を行いました。
__author__  = 'Ryota Fukuda' # mailto:25la018c@rikkyo.ac.jp
__credits__ = 'Tasuku Hayashi' # mailto:tasuku.hayashi@riken.jp
__url__     = 'https://colab.research.google.com/drive/1rRyx1eIyx2i56KETObivz26Km9cdnt5h?authuser=1'

# インポート
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.cm as cm
import lmfit as lmf
import argparse
import os
import subprocess
import shutil

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
        Pb (float) : TESのJoule発熱（フィッティング用）[W]
    """
    Pb = (G0 / n) * (Tc**n - Tb**n)
    return Pb

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
        Tbath (list) : 全熱浴温度（Tb [mK]）のリスト
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
        
        # 熱浴温度の値部分のみを抽出
        tb_key = tb_.split('_')[num]
        tbath  = float(tb_key[:-2])
        Tbath.append(tbath)
    Tbath = np.asarray(Tbath)
    print(f"Tbath values detected: {Tbath}")
    return Tbath


# ----------------------------------------------------
# フィッティング関数
# ----------------------------------------------------
def fit_pfunc(output_dir, Tbath, Rtes, Ptes, repmin, repmax, Tcmin, Tcmax, Tcfit, G0min, G0max, G0fit, nmin, nmax, nfit):
    """"
    TESのJoule発熱と熱浴温度の関係をフィッティングする。
    フィッティングから、熱浴温度が 1K のときの熱伝導度 G0 , べき定数 n , TESの温度 Tc の各パラメータを得る。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (list)     : 全熱浴温度（Tb [mK]）のリスト
        Rtes (dict)      : 全TESの抵抗（Rtes [Ω]）の辞書型配列
        Ptes (dict)      : 全TESのJoule発熱（Ptes [W]）の辞書型配列
        repmin (float)   : TESのJoule発熱代表値の最小基準
        repmax (float)   : TESのJoule発熱代表値の最大基準
        Tcmin (float)    : TESの温度のフィッティング最小値 [K]
        Tcmax (float)    : TESの温度のフィッティング最大値 [K]
        Tcfit (float)    : TESの温度のフィッティング初期値 [K]
        G0min (float)    : 熱浴温度が 1K のときの熱伝導度のフィッティング最小値 [W/K]
        G0max (float)    : 熱浴温度が 1K のときの熱伝導度のフィッティング最大値 [W/K]
        G0fit (float)    : 熱浴温度が 1K のときの熱伝導度のフィッティング初期値 [W/K]
        nmin (float)     : べき定数のフィッティング最小値
        nmax (float)     : べき定数のフィッティング最大値
        nfit (float)     : べき定数のフィッティング初期値

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        Tb (list)         : 全熱浴温度（Tb [K]）のリスト
        Pb (list)         : 各データにおける熱浴温度の代表値（Pb [W]）のリスト
        Pberr (list)      : 熱浴温度の代表値 Pb の誤差（標準偏差）のリスト
        Tc (float)        : TESの温度 Tc の最適値 [K]
        G0 (float)        : 熱浴温度が 1K のときの熱伝導度 G0 の最適値 [W/K]
        n (float)         : べき定数 n の最適値
        description (str) : 図と計算方法の説明
    """
    description = rf"""            <p>　TESのJoule発熱\(\, P_{{\rm TES}}\, \)と熱浴温度\(\, T_{{\rm bath}}\, \)の関係</p>
            \[
            P_{{\rm TES}} = \frac{{G_0}}{{n}}(T_{{\rm TES}}^n - T_{{\rm bath}}^n)
            \]
            <p>をフィッティングする。\(P_{{\rm TES}}\, \)は，PR特性グラフからTESの超伝導転移端の\(\, P_{{\rm TES}}\)（Joule発熱が一定の領域）を平均した代表値を使う。
            フィッティングにより，TESの超伝導転移端\(\, T_{{\rm c}}\, \)，熱浴温度が\(\, 1\ {{\rm K}}\, \)のときの熱伝導度\(\, G_0\, \)，べき定数\(\, n\, \)のそれぞれの最適値が得られる。
            <br>　図は\(\, P_{{\rm TES}}\, \)（誤差付き）と\(\, T_{{\rm bath}}\, \)の関係をプロットしたものと，そのフィッティング結果。\(\, P_{{\rm TES}}\, \)の誤差は不偏標準偏差で計算した。</p>
            """

    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.basename(__file__)
    output_file = os.path.splitext(filename)[0] + "_fitting.png" 
    output_path = os.path.join(output_dir, output_file)

    # 熱浴温度の代表値および誤差を取得
    Pb    = []
    Pberr = []
    for tb in Tbath:
        rtes_   = Rtes[tb][Rtes[tb]==Rtes[tb]]
        rtesmax = rtes_.max()
        rtesmin = rtes_.min()
        Rn      = rtesmax - rtesmin
        rtes    = (Rtes[tb]-rtesmin)/Rn
        Prep    = (repmin<rtes) & (rtes<repmax)
        Pb.append(np.average(Ptes[tb][Prep]))
        Pberr.append(np.std(Ptes[tb][Prep], ddof=1))
    Tb    = Tbath / 1e3 
    Pb    = np.asarray(Pb)
    Pberr = np.asarray(Pberr)
    print(Pb)
    print(Pberr)

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
    ax.errorbar(Tb*1e3, Pb*1e12, yerr=Pberr*1e12, color='black', capsize=4, ms=5, marker='o', mec='black', mfc='black', ls='None')
    Tb_min = Tbath.min()/1e3 - 0.010
    Tb_max = Tbath.max()/1e3 + 0.010
    x = np.linspace(Tb_min, Tb_max, 100)
    ax.plot(x*1e3, pfunc(x, Tc, G0, n)*1e12, 'r--')
    ax.set_title(r"$P_{\rm b}\ {\rm vs}.\ T_{\rm bath}$")
    ax.set_xlabel(r"$T_{\rm bath}\ ({\rm mK})$", fontsize=14)
    ax.set_ylabel(r"$P_{\rm b}\ ({\rm pW})$", fontsize=14)
    ax.legend(['data', r"$T_{\rm c}=$"+rf"${Tc:.3}$"+r"${\rm K},\ G_0=$"+rf"${G0*1e9:.3}$"+r"${\rm nW/K},\ n=$"+rf"${n:.2f}$"], frameon=False)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save Fitting to: {output_path}")
    return output_file, Tb, Pb, Pberr, chi2_min, Tc, G0, n, description


# ----------------------------------------------------
# TESパラメータの特性をプロットする関数群
# ----------------------------------------------------
def plot_IVproperty(output_dir, Tbath, Itb, Vout, Ites, Vtes):
    """"
    I-V特性グラフを作成する。
    出力電圧 vs. バイアス電流 & TESの両端にかかる電圧 vs. TESに流れる電流の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (list)     : 全熱浴温度（Tb [mK]）のリスト
        Itb (dict)       : 全バイアス電流（Itb [A]）の辞書型配列
        Vout (dict)      : 全出力電圧（Vout [V]）の辞書型配列
        Ites (dict)      : 全TESに流れる電流（Ites [A]）の辞書型配列
        Vtes (dict)      : 全TESの両端にかかる電圧（Vtes [V]）の辞書型配列
    
    Returns:
        output_file (str) : プロット結果の出力ファイル名
        description (str) : 図と計算方法の説明
    """
    description = rf"""            <p>　TESに流れる電流\(\, I_{{\rm TES}}\, \)とTESの両端にかかる電圧\(\, V_{{\rm TES}}\, \)の関係</p>
            \[
            V_{{\rm TES}} = \frac{{V_{{\rm TES}}}}{{R_{{\rm TES}}}}
            \]
            <p>をプロットする。\(R_{{\rm TES}}\, \)はTESの抵抗である。このとき，\(I_{{\rm TES}}\, \)は出力電圧\(\, V_{{\rm out}}\, \)と電流電圧変換係数\(\, \Xi\, \)を用いて</p>
            \[
            I_{{\rm TES}} \simeq \frac{{1}}{{\Xi}}V_{{\rm out}}
            \]
            <p>で求められる。\(\, \Xi\, \)は，SQUIDの入力コイル相互インダクタンス\(\, M_{{\rm in}}\, \)，フィードバックコイル相互インダクタンス\(\, M_{{\rm FB}}\, \)，フィードバック抵抗値\(\, R_{{\rm FB}}\, \)により</p>
            \[
            \Xi = \frac{{M_{{\rm in}}}}{{M_{{\rm FB}}}}R_{{\rm FB}}
            \]
            <p>で表される係数。また，\(R_{{\rm TES}}\, \)は測定バイアス電流\(\, I_{{\rm bias}}\, \)シャント抵抗\(\, R_{{\rm sh}}\, \)を用いて</p>
            \[
            R_{{\rm TES}} = \left(\frac{{I_{{\rm bias}}}}{{I_{{\rm TES}}}}-1\right)R_{{\rm sh}}
            \]
            <p>で書ける。<br>　左上図は各熱浴温度\(\, T_{{\rm bath}}\, \)に対する測定結果の\(\, V_{{\rm out}}\, \)と\(\, I_{{\rm bias}}\, \)の関係をプロットしたもので，右上図はそのグラフを超伝導転移端付近で拡大したものである。
            また，左下図は各熱浴温度\(\, T_{{\rm bath}}\, \)に対する計算結果の\(\, V_{{\rm TES}}\, \)と\(\, I_{{\rm TES}}\, \)の関係をプロットしたもので，右下図はそのグラフを超伝導転移端付近で拡大したものである。</p>"""

    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.basename(__file__)
    output_file = os.path.splitext(filename)[0] + "_IVproperty.png" 
    output_path = os.path.join(output_dir, output_file)

    _, ax = plt.subplots(2, 2, figsize=(12, 16))
    cmap  = plt.get_cmap('jet', len(Tbath))
    # 出力電圧 vs. バイアス電流
    for i, tb in enumerate(Tbath):
        ax[0,0].plot(Itb[tb]*1e6, Vout[tb], '.', color=cmap(i), label=r"$T_{\rm bath}=$"+str(int(tb))+r"${\rm mK}$")
    ax[0,0].set_title(r"$V_{\rm out}\ {\rm vs}.\ I_{\rm TES_{bias}}$")
    ax[0,0].set_xlabel(r"$I_{\rm TES_{bias}}\ ({\rm\mu A})$", fontsize=14)
    ax[0,0].set_ylabel(r"$V_{\rm out}\ ({\rm V})$", fontsize=14)
    ax[0,0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=3, frameon=False)
    ax[0,0].grid(ls='--')

    # TESの両端にかかる電圧 vs. TESに流れる電流
    for i, tb in enumerate(Tbath):
        ax[1,0].plot(Vtes[tb]*1e6, Ites[tb]*1e6, '.', color=cmap(i), label=r"$T_{\rm bath}=$"+str(int(tb))+r"${\rm mK}$")
    ax[1,0].set_title(r"$I_{\rm TES}\ {\rm vs.}\ V_{\rm TES}$")
    ax[1,0].set_xlabel(r"$V_{\rm TES}\ ({\rm\mu V})$", fontsize=14)
    ax[1,0].set_ylabel(r"$I_{\rm TES}\ ({\rm\mu A})$", fontsize=14)
    ax[1,0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=3, frameon=False)
    ax[1,0].grid(ls='--')

    # 出力電圧 vs. バイアス電流 (超伝導転移部分の拡大)
    for i, tb in enumerate(Tbath):
        ax[0,1].plot(Itb[tb]*1e6, Vout[tb], '.', color=cmap(i), label=r"$T_{\rm bath}=$"+str(int(tb))+r"${\rm mK}$")
    ax[0,1].set_xlim(50., 200.)
    ax[0,1].set_ylim(0., 3.)
    ax[0,1].set_title(r"$V_{\rm out}\ {\rm vs}.\ I_{\rm TES_{bias}}$")
    ax[0,1].set_xlabel(r"$I_{\rm TES_{bias}}\ ({\rm\mu A})$", fontsize=14)
    ax[0,1].set_ylabel(r"$V_{\rm out}\ ({\rm V})$", fontsize=14)
    ax[0,1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=3, frameon=False)
    ax[0,1].grid(ls='--')

    # TESの両端にかかる電圧 vs. TESに流れる電流 (超伝導転移部分の拡大)
    for i, tb in enumerate(Tbath):
        ax[1,1].plot(Vtes[tb]*1e6, Ites[tb]*1e6, '.', color=cmap(i), label=r"$T_{\rm bath}=$"+str(int(tb))+r"${\rm mK}$")
    ax[1,1].set_xlim(-.1, .5)
    ax[1,1].set_ylim(40., 150.)
    ax[1,1].set_title(r"$I_{\rm TES}\ {\rm vs.}\ V_{\rm TES}$")
    ax[1,1].set_xlabel(r"$V_{\rm TES}\ ({\rm\mu V})$", fontsize=14)
    ax[1,1].set_ylabel(r"$I_{\rm TES}\ ({\rm\mu A})$", fontsize=14)
    ax[1,1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=3, frameon=False)
    ax[1,1].grid(ls='--')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save I-V property to: {output_path}")
    return output_file, description

def plot_PRproperty(output_dir, Tbath, Rtes, Ptes, repmin, repmax, xmin, xmax, ymin, ymax):
    """"
    P-R特性グラフを作成する。
    TESのJoule発熱 vs. TESの抵抗の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (list)     : 全熱浴温度（Tb [mK]）のリスト
        Rtes (dict)      : 全TESの抵抗（Rtes [Ω]）の辞書型配列
        Ptes (dict)      : 全TESのJoule発熱（Ptes [W]）の辞書型配列
        repmin (float)   : TESのJoule発熱代表値の最小基準
        repmax (float)   : TESのJoule発熱代表値の最大基準
        xmin (float)     : 描画するx軸の最小範囲
        xmax (float)     : 描画するx軸の最大範囲
        ymin (float)     : 描画するy軸の最小範囲
        ymax (float)     : 描画するy軸の最大範囲
    
    Returns:
        output_file (str) : プロット結果の出力ファイル名
        description (str) : 図と計算方法の説明
    """
    description = rf"""            <p>　TESの抵抗\(\, R_{{\rm TES }}\, \)とTESのJoule発熱\(\, P_{{\rm TES}}\, \)の関係をプロットする。ここでの\(R_{{\rm TES}}\, \)は，正規化したTESの抵抗とする。
            <br>　左図は\(\, R_{{\rm TES }}\, \)と\(\, P_{{\rm TES}}\, \)の関係をプロットしたもので，右図はそのグラフを超伝導転移端付近で拡大したものである。</p>"""

    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.basename(__file__)
    output_file = os.path.splitext(filename)[0] + "_PRproperty.png" 
    output_path = os.path.join(output_dir, output_file)

    _, ax = plt.subplots(1, 2, figsize=(12, 8))
    cmap  = plt.get_cmap('jet', len(Tbath))
    # TESのJoule発熱 vs. TESの抵抗
    Pb = []
    for i, tb in enumerate(Tbath):
        rtes    = Rtes[tb][Rtes[tb]==Rtes[tb]]
        rtesmax = rtes.max()
        rtesmin = rtes.min()
        Rn      = rtesmax - rtesmin
        Prep    = ()
        ax[0].plot((Rtes[tb]-rtesmin)/Rn, Ptes[tb]*1e9, '.', color=cmap(i), label=r"$T_{\rm bath}=$"+str(int(tb))+r"${\rm mK}$")
    # ax[0].set_xlim(xmin, xmax)
    ax[0].set_ylim(ymin, ymax)
    ax[0].set_title(r"$P_{\rm TES}\ {\rm vs}.\ R_{\rm TES_{bias}}$")
    ax[0].set_xlabel(r"$R_{\rm TES}/R_N$", fontsize=14)
    ax[0].set_ylabel(r"$P_{\rm TES}\ ({\rm nW})$", fontsize=14)
    ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=3, frameon=False)
    ax[0].grid(ls='--')

    # TESのJoule発熱 vs. TESの抵抗
    for i, tb in enumerate(Tbath):
        rtes    = Rtes[tb][Rtes[tb]==Rtes[tb]]
        rtesmax = rtes.max()
        rtesmin = rtes.min()
        Rn      = rtesmax - rtesmin
        ax[1].plot((Rtes[tb]-rtesmin)/Rn, Ptes[tb]*1e9, '.', color=cmap(i), label=r"$T_{\rm bath}=$"+str(int(tb))+r"${\rm mK}$")
    ax[1].set_xlim(repmin, repmax)
    ax[1].set_ylim(0., 0.02)
    ax[1].set_title(r"$P_{\rm TES}\ {\rm vs}.\ R_{\rm TES_{bias}}$")
    ax[1].set_xlabel(r"$R_{\rm TES}/R_N$", fontsize=14)
    ax[1].set_ylabel(r"$P_{\rm TES}\ ({\rm nW})$", fontsize=14)
    ax[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=3, frameon=False)
    ax[1].grid(ls='--')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save R-P property to: {output_path}")
    return output_file, description

def plot_RTproperty(output_dir, Tbath, Rtes, Ptes, G0, n):
    """"
    R-T特性グラフを作成する。
    TESの抵抗 vs. TESの温度の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (list)     : 全熱浴温度（Tb [mK]）のリスト
        Rtes (dict)      : 全TESの抵抗（Rtes [Ω]）の辞書型配列
        Ptes (dict)      : 全TESのJoule発熱（Ptes [W]）の辞書型配列
        G0 (float)       : 熱浴温度が 1K のときの熱伝導度 [W/K]
        n (float)        : べき定数

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        Ttes (float)      : TESの温度 [K]
        description (str) : 図と計算方法の説明
    """
    description = rf"""            <p>　TESの抵抗\(\, R_{{\rm TES }}\, \)とTESの温度\(\, T_{{\rm TES}}\, \)の関係をプロットする。\(T_{{\rm TES}}\, \)は，フィッティング関数を逆算して</p>
            \[
            T_{{\rm TES}} = \left(T_{{\rm bath}}^n + \frac{{n\cdot P_{{\rm TES}}}}{{G_0}}\right)^{{1/n}}
            \]
            <p>で計算できる。<br>　図は\(\, R_{{\rm TES }}\, \)と\(\, T_{{\rm TES}}\, \)の関係をプロットしたものである。</p>"""

    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.basename(__file__)
    output_file = os.path.splitext(filename)[0] + "_RTproperty.png" 
    output_path = os.path.join(output_dir, output_file)

    # TESの抵抗 vs. TESの温度
    _, ax = plt.subplots(1, 1, figsize=(6, 8))
    cmap  = plt.get_cmap('jet', len(Tbath))
    Ttes = {}
    for i, tb in enumerate(Tbath):
        Ttes[tb] = cal_Ttes(tb/1e3, G0, n, Ptes[tb])
        ax.plot(Ttes[tb]*1e3, Rtes[tb], '.', color=cmap(i), label=r"$T_{\rm bath}=$"+str(int(tb))+r"${\rm mK}$")
    ax.set_title(r"$R_{\rm TES}\ {\rm vs}.\ T_{\rm TES}$")
    ax.set_xlabel(r"$T_{\rm TES}\ ({\rm mK})$", fontsize=14)
    ax.set_ylabel(r"$R_{\rm TES}\ ({\rm\Omega})$", fontsize=14)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=3, frameon=False)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save R-T property to: {output_path}")
    return output_file, Ttes, description

def plot_GTproperty(output_dir, Tbath, Ttes, G0, n):
    """"
    G-T特性グラフを作成する。
    TESの熱伝導度 vs. TESの温度の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tbath (list)     : 全熱浴温度（Tb [mK]）のリスト
        Ttes (dict)      : 全TESの温度（Ttes [K]）の辞書型配列
        G0 (float)       : 熱浴温度が 1K のときの熱伝導度 [W/K]
        n (float)        : べき定数

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        Gtes (float)      : TESの熱伝導度 [W/K]
        description (str) : 図と計算方法の説明
    """
    description = rf"""            <p>　TESの熱伝導度\(\, G_{{\rm TES }}\, \)とTESの温度\(\, T_{{\rm TES}}\, \)の関係</p>
            \[
            G_{{\rm TES}} = G_0T_{{\rm TES}}^{{n-1}}
            \]
            </p>をプロットする。\(G_0\, \)と\(\, n\, \)は，フィッティング結果の値を用いる。
            <br>　図は\(\, G_{{\rm TES }}\, \)と\(\, T_{{\rm TES}}\, \)の関係をプロットしたものである。</p>"""

    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.basename(__file__)
    output_file = os.path.splitext(filename)[0] + "_GTproperty.png" 
    output_path = os.path.join(output_dir, output_file)

    # TESの熱伝導度 vs. TESの温度
    _, ax = plt.subplots(1, 1, figsize=(6, 8))
    cmap  = plt.get_cmap('jet', len(Tbath))
    Gtes = {}
    for i, tb in enumerate(Tbath):
        Gtes[tb] = cal_Gtes(G0, n, Ttes[tb])
        ax.plot(Ttes[tb]*1e3, Gtes[tb], '.', color=cmap(i), label=r"$T_{\rm bath}=$"+str(int(tb))+r"${\rm mK}$")
    ax.set_title(r"$G_{\rm TES}\ {\rm vs}.\ T_{\rm TES}$")
    ax.set_xlabel(r"$T_{\rm TES}\ ({\rm mK})$", fontsize=14)
    ax.set_ylabel(r"$G_{\rm TES}\ ({\rm W/K})$", fontsize=14)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncols=3, frameon=False)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save G-T property to: {output_path}")
    return output_file, Gtes, description

def plot_alpha(output_dir, Itb, Ttes, Rtes, alpha_Tb, xmin, xmax, ymin, ymax):
    """"
    TESの感度グラフを作成する。
    TESの感度 vs. バイアス電流の関係をプロットする。

    Args:
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Itb (dict)       : 全バイアス電流（Itb [A]）の辞書型配列
        Ttes (dict)      : 全TESの温度（Ttes [K]）の辞書型配列
        Rtes (dict)      : 全TESの抵抗（Rtes [Ω]）の辞書型配列
        alpha_Tb (int)   : 感度を計算する対象の熱浴温度 [K]
        xmin (float)     : 描画するx軸の最小範囲
        xmax (float)     : 描画するx軸の最大範囲
        ymin (float)     : 描画するy軸の最小範囲
        ymax (float)     : 描画するy軸の最大範囲

    Returns:
        output_file (str) : プロット結果の出力ファイル名
        alpha (float)     : TESの感度
        description (str) : 図と計算方法の説明
    """
    description = rf"""            <p>　TESの感度\(\, \alpha\, \)</p>
            \[
            \alpha = \frac{{T}}{{R}}\frac{{{{\rm d}}R}}{{{{\rm d}}T}}
            \]
            </p>をプロットする。TESの温度とTESの抵抗はともに\(\, T_{{\rm bath}}={alpha_Tb:.0f}\ {{\rm mK}}\, \)のときの値を用いた。
            <br>　図は\(\, T_{{\rm bath}}={alpha_Tb:.0f}\ {{\rm mK}}\, \)のときの\(\, \alpha\, \)と\(\, I_{{\rm bias}}\, \)の関係をプロットしたものである。</p>"""

    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.basename(__file__)
    output_file = os.path.splitext(filename)[0] + "_alpha.png" 
    output_path = os.path.join(output_dir, output_file)

    # TESの感度 vs. バイアス電流
    _, ax = plt.subplots(1, 1, figsize=(6, 4))
    alpha = cal_alpha(Ttes[alpha_Tb], Rtes[alpha_Tb])
    Ibias = (Itb[alpha_Tb][1:] + Itb[alpha_Tb][:-1]) / 2
    ax.plot(Ibias*1e6, alpha)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title("Sensitivity "+r"$\alpha$"+" at "+rf"{alpha_Tb:.0f}"+r"${\rm mK}$")
    ax.set_xlabel(r"$I_{\rm bias}\ ({\rm\mu A})$", fontsize=14)
    ax.set_ylabel(r"$\alpha$", fontsize=14)
    ax.grid(ls='--')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save TES Sensitivity to: {output_path}")
    return output_file, alpha, description


# ----------------------------------------------------
# コントアを計算する関数
# ----------------------------------------------------
def cal_contour(output_dir, Tb, Pb, Perr, Tc, G0, n, chi2_min, xn1, xn2, xg0, yg0, ytc1, ytc2, cslevels):
    """"
    フィッティングパラメータのコントアを計算する。
    G0 vs. n  &  Tc vs. n  &  Tc vs. G0 のそれぞれについて計算する。

    Args:
        file (str)       : 得られた測定データのパス
        output_dir (str) : プロット結果の出力先ディレクトリの名前
        Tb (list)        : 全熱浴温度（Tb [K]）のリスト
        Pb (list)        : 各データにおける熱浴温度の代表値（Pb [W]）のリスト
        Pberr (list)     : 熱浴温度の代表値 Pb の誤差（標準偏差）のリスト
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
        description (str) : 図と計算方法の説明
    """
    description = rf"""            <p>　フィッティングパラメータのコントアをプロットする。信頼範囲の定義として，自由度2のカイ二乗分布における累積確率の値からパーセント点を求める。
            累積確率は\(\, 68.27\ \%,\ 90\ \%,\ 99\ \%\, \)の3点のP値を考えることとする。<br>　図は\(\, n,\ G_0,\ T_{{\rm c}}\, \)のそれぞれを2つのパラメータでコントアをプロットしたものである。</p>"""

    # プロット結果の出力先ディレクトリおよびファイル名を設定
    filename = os.path.basename(__file__)
    output_file = os.path.splitext(filename)[0] + "_contour.png" 
    output_path = os.path.join(output_dir, output_file)

    _, ax = plt.subplots(1, 3, figsize=(18, 4))

    # G0 と n のコントアを計算
    Xn1, YG0 = np.meshgrid(xn1, yg0)
    for mn1, nn1 in enumerate(xn1):
        chi2 = np.array([])
        for mg01, ng01 in enumerate(yg0):
            ng01 = ng01*1e-9
            chi2_cont = sum(((Pb - pfunc(Tb=Tb, Tc=Tc, G0=ng01, n=nn1))/Perr)**2)
            print(f'\r{mn1+1}/{len(Xn1[0])},{mg01+1}/{len(YG0.T[0])}',end='',flush=True)
            #print(f'\r{mn1+1}/{len(Xn1[0])}',end='',flush=True)
            chi2 = np.hstack((chi2, chi2_cont))

        if mn1 == 0:
            chi2_array = chi2
        else:
            chi2_array = np.vstack((chi2_array, chi2))
    print()
    chi2_array -= chi2_min

    # プロット
    cont0 = ax[0].contour(Xn1, YG0, chi2_array.T, levels=cslevels)
    cont0.clabel(fmt='%1.1f', fontsize=14)
    ax[0].plot(n, G0*1e9, 'r*', ms=10)
    ax[0].set_title(r"$G_0\ {\rm vs}.\ n$"+" Contour")
    ax[0].set_xlabel(r'$n$', fontsize=14)
    ax[0].set_ylabel(r'$G_0\ ({\rm nW/K})$', fontsize=14)

    # Tc と n のコントアを計算
    Xn2, YTc1 = np.meshgrid(xn2, ytc1)
    for mn2, nn2 in enumerate(xn2):
        chi2 = np.array([])
        for mtc1, ntc1 in enumerate(ytc1):
            ntc1 = ntc1*1e-3
            chi2_cont = sum(((Pb - pfunc(Tb=Tb, Tc=ntc1, G0=G0, n=nn2))/Perr)**2)
            print(f'\r{mn2+1}/{len(Xn2[0])},{mtc1+1}/{len(YTc1.T[0])}',end='',flush=True)
            #print(f'\r{mn2+1}/{len(Xn2[0])}',end='',flush=True)
            chi2 = np.hstack((chi2, chi2_cont))

        if mn2 == 0:
            chi2_array = chi2
        else:
            chi2_array = np.vstack((chi2_array, chi2))
    print()
    chi2_array -= chi2_min

    # プロット
    cont1 = ax[1].contour(Xn2, YTc1, chi2_array.T, levels=cslevels)
    cont1.clabel(fmt='%1.1f', fontsize=14)
    ax[1].plot(n, Tc*1e3, 'r*', ms=10)
    ax[1].set_title(r"$T_{\rm c}\ {\rm vs}.\ n$"+" Contour")
    ax[1].set_xlabel(r'$n$', fontsize=14)
    ax[1].set_ylabel(r'$T_{\rm c}\ ({\rm mK})$', fontsize=14)

    # Tc と G0 のコントアを計算
    XG0, YTc2 = np.meshgrid(xg0, ytc2)
    for mg02, ng02 in enumerate(xg0):
        chi2 = np.array([])
        ng02 = ng02*1e-9
        for mtc2, ntc2 in enumerate(ytc2):
            ntc2 = ntc2*1e-3
            chi2_cont = sum(((Pb - pfunc(Tb=Tb, Tc=ntc2, G0=ng02, n=n))/Perr)**2)
            print(f'\r{mg02+1}/{len(XG0[0])},{mtc2+1}/{len(YTc2.T[0])}',end='',flush=True)
            #print(f'\r{mg0+1}/{len(XG0[0])}',end='',flush=True)
            chi2 = np.hstack((chi2, chi2_cont))

        if mg02 == 0:
            chi2_array = chi2
        else:
            chi2_array = np.vstack((chi2_array, chi2))
    print()
    chi2_array -= chi2_min

    # プロット
    cont2 = ax[2].contour(XG0, YTc2, chi2_array.T, levels=cslevels)
    cont2.clabel(fmt='%1.1f', fontsize=14)
    ax[2].plot(G0*1e9, Tc*1e3, 'r*', ms=10)
    ax[2].set_title(r"$T_{\rm c}\ {\rm vs}.\ G_0$"+" Contour")
    ax[2].set_xlabel(r'$G_0\ ({\rm nW/K})$', fontsize=14)
    ax[2].set_ylabel(r'$T_{\rm c}\ ({\rm mK})$', fontsize=14)

    # 全コントアの計算結果をプロットし、保存
    plt.tight_layout()
    plt.savefig(output_path, dpi=300) 
    plt.close()
    print(f"Save Fitting Parameter Contour to: {output_path}")
    return output_file, description


# ----------------------------------------------------
# 解析結果をまとめる関数
# ----------------------------------------------------
def export_results(info, output_dir, output_format):
    """
    解析結果をHTML形式ないしはLaTeX (PDF) 形式でまとめる。
    TESパラメータの特性グラフのプロット結果を順に表示する。

    Args:
        info (dict)         : 解析結果の情報まとめ
        output_dir (str)    : 解析結果の出力先ディレクトリ
        output_format (str) : 出力フォーマット（HTML形式またはLaTeX (PDF) 形式）
    """
    # HTML出力
    if output_format == 'html':
        output_path = os.path.join(output_dir, "IVtes_result.html")
        with open(output_path, 'w') as f:
            f.write("<!DOCTYPE html>\n")
            f.write("<html lang='ja'>\n")
            f.write("<head>\n")
            f.write("    <meta charset='UTF-8'>\n")
            f.write("    <meta name='viewport' content='width=device-width, initial-scale=1.0'>\n")
            f.write("    <title>TES IV Curve Analysis Results</title>\n")
            f.write("    <style>\n")
            f.write("        body { font-family: 'Times New Roman', Times, serif; margin: 20px; line-height: 1.6; }\n")
            f.write("        h1 { font-size: 2.5rem; color: #333; }\n")
            f.write("        h2 { font-size: 2rem; color: #333; }\n")
            f.write("        h3 { font-size: 1.5rem; color: #333; }\n")
            f.write("        table { width: 100%; border-collapse: collapse; margin-top: 20px; }\n")
            f.write("        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n")
            f.write("        th { background-color: #f2f2f2; }\n")
            f.write("        .section { margin-top: 40px; border-top: 1px solid #ccc; padding-top: 20px; }\n")
            f.write("        .section img { max-width: 50%; height: auto; display: block; margin: auto; border: 1px solid #eee; box-shadow: 2px 2px 5px rgba(0,0,0,0.1); }\n")
            f.write("        .section a { display: inline-block; }\n")
            f.write("        .section .plot-container { text-align: center; }\n")
            f.write("        .footer { text-align: center; margin-top: 50px; font-size: 0.8em; color: #777; }\n")
            f.write("        .parameter-box {\n")
            f.write("            border: 1px solid #eee;\n")
            f.write("            padding: 15px;\n")
            f.write("            margin-bottom: 20px;\n")
            f.write("            background-color: #f9f9f9;\n")
            f.write("            border-radius: 5px;\n")
            f.write("        }\n")
            f.write("        .parameter-box h3 { margin-top: 0; }\n")
            f.write("        .parameter-box p { margin: 5px 0; }\n")
            f.write("        hr {\n")
            f.write("            border: none;\n")
            f.write("            height: 1px;\n")
            f.write("            background-color: #ccc;\n")
            f.write("            margin: 20px 0;\n")
            f.write("        }\n")
            f.write("    </style>\n")
            f.write("    <script src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>\n")
            f.write("    <link rel='stylesheet' href='css/smart-style.css'>\n")
            f.write("    <link rel='stylesheet' href='css/pc-style.css' media='screen and (min-width: 768px)'>\n")
            f.write("</head>\n")
            f.write("<body>\n")
            f.write(f"    <h1>{info['title']} for {'_'.join(output_dir.split('_')[:-1])}</h1>\n")
            f.write(f"    <p><strong>Author:</strong> {info['author']}</p>\n")
            f.write(info['summary'])

            f.write("    <div class='section'>\n")
            f.write("        <h2>Fitting Results</h2>\n")
            f.write(info['parameter'])
            f.write("        <div class='parameter-box'>\n")
            f.write(f"            <p><strong>Tc (TES Temperature):</strong> {info['Tc_best']*1e3:.3f} mK</p>\n")
            f.write(f"            <p><strong>G0 (Thermal Conductivity at 1K):</strong> {info['G0_best']*1e9:.3f} nW/K</p>\n")
            f.write(f"            <p><strong>n (Power Constant):</strong> {info['n_best']:.2f}</p>\n")
            f.write(f"            <p><strong>Chi-squared (Minimum):</strong> {info['chi2_min']:.2f}</p>\n")
            f.write("        </div>\n")
            f.write("    </div>\n")

            f.write("    <div class='section'>\n")
            f.write("        <h2>Result Plots</h2>\n")
            for i, (plot_file, description) in enumerate(zip(info['plots'], info['descriptions'])):
                if plot_file and description:
                    if i != 0:
                        f.write("        <hr>\n")
                    f.write(f"        <h3>{os.path.basename(plot_file).replace('.png', '').split('_', 1)[1]}</h3>\n")
                    f.write(description + "\n")
                    f.write("            <div class='plot-container'>\n")
                    f.write(f"                <a href='./{os.path.basename(plot_file)}' target='_blank'><img src='./{os.path.basename(plot_file)}' alt='{os.path.basename(plot_file)}'></a>\n")
                    f.write("            </div>\n")
            f.write("    </div>\n")
            f.write("</body>\n")
            f.write("</html>\n")
        print(f"\n\033[32mHTML report saved to: {output_path}\033[0m")

    # LaTeX (PDF) 出力
    elif output_format == 'pdf':
        output_path = os.path.join(output_dir, "IVtes_results.tex")
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

            f.write(r"\title{" + info['title'].replace('_', r'\_') + "}" + "\n")
            f.write(r"\author{}" + "\n")
            f.write(r"\vspace{-5zw}" + "\n")
            f.write(r"\date{\today}" + "\n")
            f.write("\n")
            f.write(r"\begin{document}" + "\n")
            f.write(r"\maketitle" + "\n")
            f.write(r"本解析では，TESの$I-V$測定で得られたデータから，TESカロリメータの諸パラメータを決定する。" + "\n")
            f.write(r"\par" + "\n")
            f.write(r"測定では，バイアス電流$I_{\rm bias}$に対する出力電圧$V_{\rm out}$を調べる。これについて，まずは電流電圧変換係数$\Xi$を用いてTESの$I-V$特性を得る。さらに，異なる熱浴温度$T_{\rm bath}$における$I-V$特性を求めることで，フィッティングにより熱伝導率$G$や温度依存性のべき定数$n$，TESの温度$T_{\rm TES}$が決定できる。また，TESの$R-T$特性を調べることで転移温度$T_{\rm c}$がわかり，温度感度$\alpha$が計算できる。" + "\n")
            f.write("\n")

            f.write(r"\section*{Fitting Parameters}" + "\n")
            f.write(r"フィッティング結果のパラメータをまとめる。" + "\n")
            f.write(r"\begin{itemize}" + "\n")
            f.write(r"    \item \textbf{Tc (TES Temperature)}: $" + f"{info['Tc_best']:.3f}" + r"\ {\rm K}$" + "\n")
            f.write(r"    \item \textbf{G0 (Thermal Conductivity at 1K)}: $" + f"{info['G0_best']*1e9:.3f}" + r"\ {\rm nW/K}$" + "\n")
            f.write(r"    \item \textbf{n (Power Constant)}: $" + f"{info['n_best']:.2f}" + r"$" + "\n")
            f.write(r"    \item \textbf{Chi-squared (Minimum)}: $" + f"{info['chi2_min']:.2f}" + r"$" + "\n")
            f.write(r"\end{itemize}" + "\n")
            f.write(r"\clearpage" + "\n")
            f.write("\n")
            
            f.write(r"\section*{Plot Results}" + "\n")
            for plot_file, description in zip(info['plots'], info['descriptions']):
                if plot_file and description:
                    latex_desc = description.replace(r'<p>　', '').replace(r'<p>', '').replace(r'</p>', r'').replace(r'<br>', r'\\').replace(r'<hr>', r'').replace(r'</hr>', r'')
                    latex_desc = latex_desc.replace(r'\(\, ', r'$').replace(r'\(', r'$').replace(r'\, \)', r'$').replace(r'\)', r'$')
                    latex_desc = latex_desc.replace(r'            ', r'')
                    
                    if plot_file != info['plots'][-1]:
                        f.write(r"\subsection*{" + os.path.basename(plot_file).replace('.png', '').replace('_', r'\_') + "}" + "\n")
                        f.write(latex_desc + "\n")
                        f.write(r"\begin{figure}[H]" + "\n")
                        f.write(r"    \centering" + "\n")
                        f.write(r"    \includegraphics[width=0.8\textwidth]{" + os.path.basename(plot_file) + "}" + "\n")
                        f.write(r"    \label{fig:" + os.path.basename(plot_file).split('.')[0].replace('_', '') + "}" + "\n")
                        f.write(r"\end{figure}" + "\n")
                        f.write(r"\clearpage" + "\n")
                        f.write("\n")
                    else:
                        f.write(r"\subsection*{" + os.path.basename(plot_file).replace('.png', '').replace('_', r'\_') + "}" + "\n")
                        f.write(latex_desc + "\n")
                        f.write(r"\begin{figure}[H]" + "\n")
                        f.write(r"    \centering" + "\n")
                        f.write(r"    \includegraphics[width=0.8\textwidth]{" + os.path.basename(plot_file) + "}" + "\n")
                        f.write(r"    \label{fig:" + os.path.basename(plot_file).split('.')[0].replace('_', '') + "}" + "\n")
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
    
    # 'start,stop,step' 形式の文字列をnp.arange()の引数に変換してNumPy配列を生成する関数
    def parse_arange_str(arange_str):
        try:
            parts = arange_str.split(',')
            if len(parts) != 3:
                raise ValueError("Invalid format: Must be 'start,stop,step'")
            start = float(parts[0].strip())
            stop = float(parts[1].strip())
            step = float(parts[2].strip())
            return np.arange(start, stop, step)
        
        except ValueError as e:
            raise argparse.ArgumentTypeError(f"Invalid arange format for '{arange_str}': {e}. Expected 'start,stop,step'.")

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
    parser.add_argument('--repmin', type=float, default=0.4,
                        help="Minimum standard for typical Joule heating of TES. (Default: 0.4 µA)")
    parser.add_argument('--repmax', type=float, default=0.6,
                        help="Maximum standard for typical Joule heating of TES. (Default: 0.6 µA)")
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
    parser.add_argument('--nmin', type=float, default=1.,
                        help="Minimum power constant fitting. (Default: 1)")
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
    parser.add_argument('--alpha_Tb', type=float, default=100.,
                        help="Tbath for which alpha vs. Ibias plot is generated. (Default: 100.0 mK)")
    parser.add_argument('--alpha_xmin', type=float, default=0.,
                        help="Minimum x-axis range for plotting alpha-Itb property. (Default: 0 µA)")
    parser.add_argument('--alpha_xmax', type=float, default=200.,
                        help="Maximum x-axis range for plotting alpha-Itb property. (Default: 200 µA)")
    parser.add_argument('--alpha_ymin', type=float, default=-20000.,
                        help="Minimum y-axis range for plotting alpha-Itb property. (Default: -20000)")
    parser.add_argument('--alpha_ymax', type=float, default=600000.,
                        help="Maximum y-axis range for plotting alpha-Itb property. (Default: 600000)")
    parser.add_argument('--xn1', type=str, default="3.5,4.5,0.01",
                        help="Scan range of x-axis (power constant) contour. (Default: \"3.5,4.5,0.01\")")
    parser.add_argument('--xn2', type=str, default="3.9,4.1,0.01",
                        help="Scan range of x-axis (power constant) contour. (Default: \"3.9,4.1,0.01\")")
    parser.add_argument('--xg0', type=str, default="141,201,1",
                        help="Scan range of x-axis (thermal conductivity at heat bath temperature 1K) contour. (Default: \"141,201,1\")")
    parser.add_argument('--yg0', type=str, default="50,350,0.1",
                        help="Scan range of y-axis (thermal conductivity at heat bath temperature 1K) contour. (Default: \"50,350,0.1\")")
    parser.add_argument('--ytc1', type=str, default="140,152,0.1",
                        help="Scan range of y-axis (TES temperature) contour. (Default: \"140,152,0.1\")")
    parser.add_argument('--ytc2', type=str, default="140,152,0.1",
                        help="Scan range of y-axis (TES temperature) contour. (Default: \"140,152,0.1\")")
    args = parser.parse_args()

    # オプション引数 --xn1 ~ --ytc2 までのメッシュ範囲配列について、文字列をNumPy配列に変換
    args.xn1 = parse_arange_str(args.xn1)
    args.xn2 = parse_arange_str(args.xn2)
    args.xg0 = parse_arange_str(args.xg0)
    args.yg0 = parse_arange_str(args.yg0)
    args.ytc1 = parse_arange_str(args.ytc1)
    args.ytc2 = parse_arange_str(args.ytc2)

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
    # Min = 1.009e-10 # H
    # Mfb = 8.912e-11 # H
    # Rfb = 30.e3     # Ω
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
        tb = Tbath[i]
        print(f"\nProcessing file: {os.path.basename(file)} (Tbath = {tb} mK)")
        try:
            # TXTファイルから出力電圧 Vout とバイアス電流 Itb を抽出、およびその他のパラメータを計算
            data = np.loadtxt(file)
            num  = len(data[:,0])
            Vout[tb] = data[:,1] - data[:,1].min()      # V
            Itb[tb]  = data[:,0] / 1e6                  # A (μA to A)
            Ites[tb] = cal_Ites(Vout[tb], Xi)           # A
            Rtes[tb] = cal_Rtes(Itb[tb], Ites[tb], Rsh) # Ω
            Vtes[tb] = cal_Vtes(Rtes[tb], Ites[tb])     # V
            Ptes[tb] = Ites[tb] * Vtes[tb]              # W

            if i == 0:
                # 超伝導での傾きを求めて、傾きが正になるように電流電圧変換係数 Xi を補正
                tilt = np.average(np.diff(Ites[tb][num-11:num-1]))/np.average(np.diff(Vtes[tb][num-11:num-1]))
                print(f'T    = {tb}')
                print(f'Xi   = {Xi}')
                print(f'tilt = {tilt}')
                j = 0
                while tilt > 0:
                    j  += 1
                    Xi -= 1
                    Ites[tb] = cal_Ites(Vout[tb], Xi)
                    Rtes[tb] = cal_Rtes(Itb[tb], Ites[tb], Rsh)
                    Vtes[tb] = cal_Vtes(Rtes[tb], Ites[tb])
                    tilt     = np.average(np.diff(Ites[tb][num-11:num-1]))/np.average(np.diff(Vtes[tb][num-11:num-1]))
                while tilt < 0:
                    j  += 1
                    Xi += 1
                    Ites[tb] = cal_Ites(Vout[tb], Xi)
                    Rtes[tb] = cal_Rtes(Itb[tb], Ites[tb], Rsh)
                    Vtes[tb] = cal_Vtes(Rtes[tb], Ites[tb])
                    tilt     = np.average(np.diff(Ites[tb][num-11:num-1]))/np.average(np.diff(Vtes[tb][num-11:num-1]))
                print(f'{j} modified Xi   = {Xi}')
                print(f'{j} modified tilt = {tilt}')
            
            # 最終的なTESパラメータを計算
            Ites[tb] = cal_Ites(Vout[tb], Xi)
            Rtes[tb] = cal_Rtes(Itb[tb], Ites[tb], Rsh)
            Vtes[tb] = cal_Vtes(Rtes[tb], Ites[tb])
            Ptes[tb] = Ites[tb] * Vtes[tb]

        except Exception as e:
            print(f"Error processing file {file}: {e}")
            Vout[tb] = np.array([])
            Itb[tb]  = np.array([])
            Ites[tb] = np.array([])
            Rtes[tb] = np.array([])
            Vtes[tb] = np.array([])
            Ptes[tb] = np.array([])
            continue

    result_IVproperty, desc_IVproperty = plot_IVproperty(output_dir=args.output_dir,
                                        Tbath=Tbath,
                                        Itb=Itb,
                                        Vout=Vout,
                                        Ites=Ites,
                                        Vtes=Vtes
    )
    result_PRproperty, desc_PRproperty = plot_PRproperty(output_dir=args.output_dir,
                                        Tbath=Tbath,
                                        Rtes=Rtes,
                                        Ptes=Ptes,
                                        repmin=args.repmin,
                                        repmax=args.repmax,
                                        xmin=args.PR_xmin,
                                        xmax=args.PR_xmax,
                                        ymin=args.PR_ymin,
                                        ymax=args.PR_ymax
    )
    result_fitting, Tb, Pb, Pberr, chi2_min, Tc, G0, n, desc_fitting = fit_pfunc(output_dir=args.output_dir,
                                                                   Tbath=Tbath,
                                                                   Rtes=Rtes,
                                                                   Ptes=Ptes,
                                                                   repmin=args.repmin,
                                                                   repmax=args.repmax,
                                                                   Tcmin=args.Tcmin,
                                                                   Tcmax=args.Tcmax,
                                                                   Tcfit=args.Tcfit,
                                                                   G0min=args.G0min,
                                                                   G0max=args.G0max,
                                                                   G0fit=args.G0fit,
                                                                   nmin=args.nmin,
                                                                   nmax=args.nmax,
                                                                   nfit=args.nfit
    )
    result_RTproperty, Ttes, desc_RTproperty = plot_RTproperty(output_dir=args.output_dir,
                                              Tbath=Tbath,
                                              Rtes=Rtes,
                                              Ptes=Ptes,
                                              G0=G0,
                                              n=n
    )
    result_GTproperty, Gtes, desc_GTproperty = plot_GTproperty(output_dir=args.output_dir,
                                              Tbath=Tbath,
                                              Ttes=Ttes,
                                              G0=G0,
                                              n=n
    )
    result_alpha, alpha, desc_alpha = plot_alpha(output_dir=args.output_dir,
                                     Itb=Itb,
                                     Ttes=Ttes,
                                     Rtes=Rtes,
                                     alpha_Tb=args.alpha_Tb,
                                     xmin=args.alpha_xmin,
                                     xmax=args.alpha_xmax,
                                     ymin=args.alpha_ymin,
                                     ymax=args.alpha_ymax
    )
    result_contour, desc_contour = cal_contour(output_dir=args.output_dir,
                                 Tb=Tb,
                                 Pb=Pb,
                                 Perr=Pberr,
                                 Tc=Tc,
                                 G0=G0,
                                 n=n,
                                 chi2_min=chi2_min,
                                 xn1=args.xn1,
                                 xn2=args.xn2,
                                 xg0=args.xg0,
                                 yg0=args.yg0,
                                 ytc1=args.ytc1,
                                 ytc2=args.ytc2,
                                 cslevels=cslevels
    )

    # 解析結果の情報をまとめる
    summary = rf"""            <p>　本解析では，TESの\(\, I-V\, \)測定で得られたデータから，TESカロリメータの諸パラメータを決定する。
            <br>　測定では，バイアス電流\(\, I_{{{{\rm bias}}}}\, \)に対する出力電圧\(\, V_{{{{\rm out}}}}\, \)を調べる。これについて，まずは電流電圧変換係数\(\, \Xi\, \)を用いてTESの\(\, I-V\, \)特性を得る。
            さらに，異なる熱浴温度\(\, T_{{{{\rm bath}}}}\, \)における\(\, I-V\, \)特性を求めることで，フィッティングにより熱伝導率\(\, G\, \)や温度依存性のべき定数\(\, n\, \)，TESの温度\(\, T_{{{{\rm TES}}}}\, \)が決定できる。
            また，TESの\(\, R-T\, \)特性を調べることで転移温度\(\, T_{{{{\rm c}}}}\, \)がわかり，温度感度\(\, \alpha\, \)が計算できる。</p>"""
    parameter = rf"""            <p>　フィッティング結果のパラメータをまとめる。</p>"""
    info = {
        'title'     : "TES IV Curve Analysis Results",
        'author'    : __author__,
        'summary'   : summary,
        'parameter' :parameter,
        'Tc_best'   : Tc,
        'G0_best'   : G0,
        'n_best'    : n,
        'chi2_min'  : chi2_min,
        'plots'     : [result_IVproperty,
                     result_PRproperty,
                     result_fitting,
                     result_RTproperty,
                     result_GTproperty,
                     result_alpha,
                     result_contour
        ],
        'descriptions': [desc_IVproperty,
                         desc_PRproperty,
                         desc_fitting,
                         desc_RTproperty,
                         desc_GTproperty,
                         desc_alpha,
                         desc_contour
        ]
    }

    # 結果をオプション引数 --output_format で書き出し
    export_results(info, args.output_dir, args.output_format)

    print("\nCompleted TES IV Curve Analysis!\n")
