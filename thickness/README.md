<div id="top"></div>

# 膜厚解析ツール

<p align="right">最終更新：2025.08.28</br>作成者：福田 凌大&ensp;&nbsp;&nbsp;</p>


## 概要
<div id="abstract"></div>

&emsp;本ディレクトリは，理研Clean Roomの表面粗さ計（触針型）による測定データから，試料の厚みを計算するためのものです [^1] 。<br/>
&emsp;このREADMEファイル [^2] にて，本ディレクトリの構成および使い方などを説明します。


## 目次
<div id="tableofcount"></div>

- [膜厚解析ツール](#膜厚解析ツール)
  - [概要](#概要)
  - [目次](#目次)
  - [ディレクトリ構成](#ディレクトリ構成)
  - [ディレクトリ用法](#ディレクトリ用法)
    - [解析データ](#解析データ)
    - [measure.py](#measurepy)
  - [実行環境](#実行環境)


## ディレクトリ構成
<div id="section1"></div>

&emsp;本ディレクトリの全体構成は以下の折りたたみにある通りです。ディレクトリ内には，別に2つのディレクトリと，このREADMEを含めた5つのファイルがあります。下に少しまとめます。チェックありはactive，チェックなしはlegacyかoboleteを表します。

- [x] **BoardSchematic.png**：LOR3A_TESTディレクトリに格納されている解析データについて，その測定基板の位置情報が描画された図。measure.pyを走らせると，自動的に解析結果をまとめるディレクトリにコピーされます。
- [x] **LOR3A_TEST**：理研CRの表面粗さ計（触針型）で測定されたLOR3A基板の膜厚データをまとめたディレクトリ。
- [x] **LOR3A_TEST_results**：LOR3A_TESTディレクトリをmeasure.py(version1.3)にかけて解析し，その結果をまとめたディレクトリ。
- [x] **measure.py**：理研CRの表面粗さ計（触針型）で測定された膜厚測定データを解析するためのスクリプトファイル。後述 $\rightharpoonup$ [measure.py](#measurepy)
- [x] **README.md**：本ファイル。上記のディレクトリおよびファイルの説明書。
- [ ] **test.png**： obsoleteです。test.pyを実行して得られた図です。
- [ ] **test.py**：obsoleteです。膜厚解析のシミュレーション用に作成したスクリプトファイルです。

<details>
    <summary><span style="background-color: #20B2AA; ">ディレクトリ構成（クリックで折りたたみを展開）</span></summary>

    .
    ├── BoardSchematic.png
    ├── LOR3A_TEST
    │   ├── LOR3A_test1_250526_LtoR_Ea1_1.txt
    │   ├── LOR3A_test1_250526_LtoR_Ea3_1.txt
    │   ├── LOR3A_test1_250526_LtoR_Ea3_2.txt
    │   ├── LOR3A_test1_250526_LtoR_Ea4_1.txt
    │   ├── LOR3A_test1_250526_LtoR_Ea5_1.txt
    │   ├── LOR3A_test1_250526_RtoL_Ea1_1.txt
    │   ├── LOR3A_test1_250526_RtoL_Ea1_2.txt
    │   ├── LOR3A_test1_250526_RtoL_Ea2_1.txt
    │   ├── LOR3A_test1_250526_RtoL_Ea3_1.txt
    │   ├── LOR3A_test1_250526_RtoL_Ea4_1.txt
    │   ├── LOR3A_test1_250526_RtoL_Ea4_2.txt
    │   ├── LOR3A_test1_250526_RtoL_Ea5_1.txt
    │   ├── LOR3A_test1_250526_test.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea1_1.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea1_2.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea2_1.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea3_1.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea3_2.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea3_3.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea4_1.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea4_2.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea5_1.txt
    │   ├── LOR3A_test2_250526_LtoR_Ea5_2.txt
    │   ├── LOR3A_test2_250526_RtoL_Ea1_1.txt
    │   ├── LOR3A_test2_250526_RtoL_Ea1_2.txt
    │   ├── LOR3A_test2_250526_RtoL_Ea1_3.txt
    │   ├── LOR3A_test2_250526_RtoL_Ea2_1.txt
    │   ├── LOR3A_test2_250526_RtoL_Ea3_1.txt
    │   ├── LOR3A_test2_250526_RtoL_Ea3_2.txt
    │   ├── LOR3A_test2_250526_RtoL_Ea4_1.txt
    │   ├── LOR3A_test2_250526_RtoL_Ea5_1.txt
    │   ├── LOR3A_test3_250526_LtoR_Ea1_1.txt
    │   ├── LOR3A_test3_250526_LtoR_Ea2_1.txt
    │   ├── LOR3A_test3_250526_LtoR_Ea2_2.txt
    │   ├── LOR3A_test3_250526_LtoR_Ea3_1.txt
    │   ├── LOR3A_test3_250526_LtoR_Ea3_2.txt
    │   ├── LOR3A_test3_250526_LtoR_Ea4_1.txt
    │   ├── LOR3A_test3_250526_LtoR_Ea5_1.txt
    │   ├── LOR3A_test3_250526_RtoL_Ea1_1.txt
    │   ├── LOR3A_test3_250526_RtoL_Ea2_1.txt
    │   ├── LOR3A_test3_250526_RtoL_Ea3_1.txt
    │   ├── LOR3A_test3_250526_RtoL_Ea3_2.txt
    │   ├── LOR3A_test3_250526_RtoL_Ea4_1.txt
    │   └── LOR3A_test3_250526_RtoL_Ea5_1.txt
    ├── LOR3A_TEST_results
    │   ├── BoardSchematic.png
    │   ├── LOR3A_test1_250526_LtoR_Ea1_1.png
    │   ├── LOR3A_test1_250526_LtoR_Ea3_1.png
    │   ├── LOR3A_test1_250526_LtoR_Ea3_2.png
    │   ├── LOR3A_test1_250526_LtoR_Ea4_1.png
    │   ├── LOR3A_test1_250526_LtoR_Ea5_1.png
    │   ├── LOR3A_test1_250526_RtoL_Ea1_1.png
    │   ├── LOR3A_test1_250526_RtoL_Ea1_2.png
    │   ├── LOR3A_test1_250526_RtoL_Ea2_1.png
    │   ├── LOR3A_test1_250526_RtoL_Ea3_1.png
    │   ├── LOR3A_test1_250526_RtoL_Ea4_1.png
    │   ├── LOR3A_test1_250526_RtoL_Ea4_2.png
    │   ├── LOR3A_test1_250526_RtoL_Ea5_1.png
    │   ├── LOR3A_test1_250526_test.png
    │   ├── LOR3A_test2_250526_LtoR_Ea1_1.png
    │   ├── LOR3A_test2_250526_LtoR_Ea1_2.png
    │   ├── LOR3A_test2_250526_LtoR_Ea2_1.png
    │   ├── LOR3A_test2_250526_LtoR_Ea3_1.png
    │   ├── LOR3A_test2_250526_LtoR_Ea3_2.png
    │   ├── LOR3A_test2_250526_LtoR_Ea3_3.png
    │   ├── LOR3A_test2_250526_LtoR_Ea4_1.png
    │   ├── LOR3A_test2_250526_LtoR_Ea4_2.png
    │   ├── LOR3A_test2_250526_LtoR_Ea5_1.png
    │   ├── LOR3A_test2_250526_LtoR_Ea5_2.png
    │   ├── LOR3A_test2_250526_RtoL_Ea1_1.png
    │   ├── LOR3A_test2_250526_RtoL_Ea1_2.png
    │   ├── LOR3A_test2_250526_RtoL_Ea1_3.png
    │   ├── LOR3A_test2_250526_RtoL_Ea2_1.png
    │   ├── LOR3A_test2_250526_RtoL_Ea3_1.png
    │   ├── LOR3A_test2_250526_RtoL_Ea3_2.png
    │   ├── LOR3A_test2_250526_RtoL_Ea4_1.png
    │   ├── LOR3A_test2_250526_RtoL_Ea5_1.png
    │   ├── LOR3A_test3_250526_LtoR_Ea1_1.png
    │   ├── LOR3A_test3_250526_LtoR_Ea2_1.png
    │   ├── LOR3A_test3_250526_LtoR_Ea2_2.png
    │   ├── LOR3A_test3_250526_LtoR_Ea3_1.png
    │   ├── LOR3A_test3_250526_LtoR_Ea3_2.png
    │   ├── LOR3A_test3_250526_LtoR_Ea4_1.png
    │   ├── LOR3A_test3_250526_LtoR_Ea5_1.png
    │   ├── LOR3A_test3_250526_RtoL_Ea1_1.png
    │   ├── LOR3A_test3_250526_RtoL_Ea2_1.png
    │   ├── LOR3A_test3_250526_RtoL_Ea3_1.png
    │   ├── LOR3A_test3_250526_RtoL_Ea3_2.png
    │   ├── LOR3A_test3_250526_RtoL_Ea4_1.png
    │   ├── LOR3A_test3_250526_RtoL_Ea5_1.png
    │   ├── thickness_correlation_test1_Ea1-3-5.png
    │   ├── thickness_correlation_test1_Ea2-3-4.png
    │   ├── thickness_correlation_test1.png
    │   ├── thickness_correlation_test2_Ea1-3-5.png
    │   ├── thickness_correlation_test2_Ea2-3-4.png
    │   ├── thickness_correlation_test2.png
    │   ├── thickness_correlation_test3_Ea1-3-5.png
    │   ├── thickness_correlation_test3_Ea2-3-4.png
    │   ├── thickness_correlation_test3.png
    │   ├── thickness_results.html
    │   ├── thickness_results.pdf
    │   └── thickness_results.tex
    ├── measure.py
    ├── README.md
    ├── test.png
    └── test.py

<p align="right"><a href="#section1"><strong>セクションに戻る»</strong></a></p>
</details>


## ディレクトリ用法
<div id="section2"></div>

&emsp;本ディレクトリ内のファイルおよび子ディレクトリを使うのに際して，**単に膜厚データの解析をしたい場合** は，解析したいデータが格納されたディレクトリを別途用意し，そのディレクトリと同階層にmeasure.pyファイルをコピーしてください。**膜厚データ解析の勉強をしたい場合** は，measure.pyファイルおよびLOR3Aディレクトリを同階層にコピーしてください。簡単な注意点と使い方の説明として，次のセクションも併せて確認してください。

</br>


### 解析データ
<div id="section2.1"></div>

&emsp;解析データの名前は，すべて <span style="color: #FF0000; ">"[基板名]\_test[基板番号]\_[測定年月日]\_[測定方向（LtoR or RtoL）]\_[測定位置（Ea番号）]\_[測定番号].txt"</span> のように付けてください（例：LOR3A_test1_20250526_LtoR_Ea1_1.txt）[^3] 。その上で，同じ基板名の膜厚測定データをすべてまとめて任意の名前のディレクトリに格納してください。作成したディレクトリは，解析スクリプトと同じ階層に置いてください。

> [!WARNING]
> &emsp;拡張子がtxt形式でない場合，txt形式に変えるか解析ファイルを一部修正する必要があります。

</br>


### measure.py
<div id="section2.2"></div>

&emsp;測定データの膜厚について個々とその相関を解析し，一括に結果を出力できるようにしたスクリプトです。このファイルを解析したいデータが格納されているディレクトリと同階層に置き，その階層にてターミナルで

```bash:実行コマンド
python measure.py [解析データが格納されたディレクトリ名]/*.txt
```

を入力・実行すれば解析が進行します。このとき，同階層に新しく **[解析データが格納されたディレクトリ名]\_results** という名前のディレクトリが自動で作られ，その中に解析結果がすべてpng形式で保存されます。また，結果をPDFでまとめて見るための $\TeX$ ファイルも同時に生成されます。"**thickness_results.tex**" という名前で保存されるので，別途コンパイルをしてください。</br>
&emsp;さらに，このスクリプトファイルは以下のオプションを用意してあります。

- `--exclude_files`：解析に使用しないファイルを指定します。デフォルトは未指定です。

    ```bash:使用例
    python measure.py LOR3A/*.txt --exclude_files LOR3A_test1_250526_test.txt
    ```
- `--output_dir`：解析結果の出力先ディレクトリを指定します。デフォルトは "[解析データが格納されたディレクトリ名]\_results" です。

    ```bash:使用例
    python measure.py LOR3A/*.txt --output_dir thickness_LOR3A
    ```
- `data_column`：測定データのtxtファイルのうち，解析に用いる数値の列をRaw，RawLevel，Normal，Rough，Waviの中から選択します。デフォルトはNormalです。

    ```bash:使用例
    python measure.py LOR3A/*.txt --data_column Raw
    ```
    > [!IMPORTANT]
    > &emsp;理研CRの表面粗さ計（触針型）で測定したデータの解析で必要になる選択ですが，通常，他の列は使用しません（他の列の機能は現時点で確認できていません）。</br>
    > &emsp;また，理研CRの表面粗さ計（触針型）以外の機器で測定されたデータを使用する場合は，main関数の設定変更が必要になりますので注意してください。
- `output_format`：解析結果をまとめる出力形式をPDFまたはHTMLで指定します。デフォルトはPDFで，thickness_results.texというファイルが生成されます。HTMLを指定すると，thickness_results.htmlというファイルが生成されます。

    ```bash:使用例
    python measure.py LOR3A/*.txt --output_format html
    ```
- `window_size`：uniform_filter1dにより素データのプロットをスムージングするときのwindow size（何個のデータポイントを平均化するか）を指定します。デフォルトは50です。

    ```bash:使用例
    python measure.py LOR3A/*.txt --window_size 100
    ```

> [!NOTE]
> &emsp;出力される解析結果は以下の通りです。
> - **各々の測定データのプロット**： "[測定データ名].png" という名前で出力されます。
> - **基板ごとに測定位置の厚みの相関**： "thickness\_correlation\_test[基板番号]\_[相関を見る位置の組み合わせ].png" という名前で出力されます。
> - **基板ごとの厚みの相関**： "thickness\_correlation\_test[基板番号].png" という名前で出力されます。
> - **解析結果をまとめたレポートファイル**： "thickness_result.html" または "thickness_result.tex" という名前で出力されます。
> 
> &emsp;thickness_result.texは別途コンパイルし，PDFに変換して使用してください。

</br>


## 実行環境
<div id="section3"></div>

&emsp;measure.pyファイルを用いて解析する場合，当然ですがPythonの実行環境が必須です。このとき，必要なライブラリとしてはNumPy，Scipy，Matplotlib，Pandasがあります。これらのインストールをしてから動作させてください。さらに，`--output_format pdf` オプションを使用してPDFレポートを生成する場合は，追加で $\LaTeX$ の実行環境と関連するパッケージが必要になります。これらの環境構築を事前に済ませてください。


<p align="right"><strong><a href="#top">トップへ»</a></strong></p>


[^1]: 本ツールは，理研CRの表面粗さ計（触針型）で測定したものに限らず，測定プロットが矩形波状になっているデータであれば使用可能です。
[^2]: このREADMEファイルは，GitHub仕様に微修正しています。立教サーバの方を見るときは，VS Codeで開くことを推奨します。ただし，"Markdown All in One" および "Markdown Named CodeBlocks" の2つの拡張機能をインストールしてから使用してください。
[^3]: 現状，測定データ名に対する汎用性が低いため，改良次第アップロードします。
