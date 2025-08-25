<div id="top"></div>

# TES IVカーブ解析ツール

<p align="right">最終更新：2025.08.25</br>作成者：福田 凌大&ensp;&nbsp;&nbsp;</p>


## 概要
<div id="abstract"></div>

&emsp;本ディレクトリは，超伝導転移端センサ（Transition Edge Sensor; TES）の $I-V$ 測定から得られたデータを解析するためのものです。<br/>
&emsp;このREADMEファイル [^1]にて，本ディレクトリの構成および使い方などを説明します。


## 目次
<div id="tableofcount"></div>

- [TES IVカーブ解析ツール](#tes-ivカーブ解析ツール)
  - [概要](#概要)
  - [目次](#目次)
  - [ディレクトリ構成](#ディレクトリ構成)
  - [ディレクトリ用法](#ディレクトリ用法)
    - [解析データ](#解析データ)
    - [ivana.ipynb](#ivanaipynb)
      - [ivana.ipynbファイル](#ivanaipynbファイル)
      - [ivana\_re.ipynbファイル](#ivana_reipynbファイル)
    - [IVtes.py](#ivtespy)
    - [IVtes\_class.py](#ivtes_classpy)
  - [実行環境](#実行環境)


## ディレクトリ構成
<div id="section1"></div>

&emsp;本ディレクトリの全体構成は以下の折りたたみにある通りです。ディレクトリ内には，別に6つのディレクトリと，このREADMEを含めた5つのファイルがあります。下に少しまとめます。チェックありはactive，チェックなしはlegacyかoboleteを表します。

- [ ] **230198**：産総研にて測定されたJAXA120 Sa4の $I-V$ データをまとめたディレクトリ。
- [ ] **230918_results**：230918ディレクトリをIVtes.py(version1.0)にかけて解析し，その結果をまとめたディレクトリ。
- [x] **230918_subset**：230918ディレクトリの全データから，超伝導転移近傍のデータを抽出してまとめたディレクトリ。
- [x] **230918_subset_results**：230918_subsetディレクトリをIVtes.py(version1.2)にかけて解析し，その結果をまとめたディレクトリ。
- [ ] **ivana_re.ipynb**：ivana.ipynbファイルをVS Code環境用に再構築したノートブックファイル。後述 $\rightharpoonup$ [ivana.ipynb](#ivanaipynb)
- [ ] **ivana.ipynb**： $I-V$ データ解析のためのノートブックファイル。後述 $\rightharpoonup$ [ivana.ipynb](#ivanaipynb)
- [ ] **IVdata_A5**：宇宙研にて測定されたJAXA120 Ea4 A5の $I-V$ データをまとめたディレクトリ。
- [ ] **IVdata_A5_results**：IVdata_A5ディレクトリをIVtes.py(version1.0)にかけて解析し，その結果をまとめたディレクトリ。
- [ ] **IVtes_class.py**：IVtes.pyをオブジェクト指向で書き直した解析ファイル。後述 $\rightharpoonup$ [IVtes_class.py](#ivtes_classpy)
- [x] **IVtes.py**： $I-V$ データ解析のためのスクリプトファイル。後述 $\rightharpoonup$ [IVtes.py](#ivtespy)
- [x] **README.md**：本ファイル。上記のディレクトリおよびファイルの説明書。

<details>
    <summary><span style="background-color: #20B2AA; "><b>ディレクトリ構成（クリックで折りたたみを展開）</b></span></summary>

    .
    ├── 230918
    │   ├── IV_070mK_20230918.txt
    │   ├── IV_075mK_20230918.txt
    │   ├── IV_080mK_20230918.txt
    │   ├── IV_085mK_20230918.txt
    │   ├── IV_090mK_20230918.txt
    │   ├── IV_095mK_20230918.txt
    │   ├── IV_100mK_20230918.txt
    │   ├── IV_105mK_20230918.txt
    │   ├── IV_110mK_20230918.txt
    │   ├── IV_115mK_20230918.txt
    │   ├── IV_120mK_20230918.txt
    │   ├── IV_125mK_20230918.txt
    │   ├── IV_130mK_20230918.txt
    │   ├── IV_135mK_20230918.txt
    │   ├── IV_140mK_20230918.txt
    │   ├── IV_145mK_20230918.txt
    │   ├── IV_150mK_20230918.txt
    │   ├── IV_155mK_20230918.txt
    │   ├── IV_160mK_20230918.txt
    │   ├── IV_165mK_20230918.txt
    │   ├── IV_170mK_20230918.txt
    │   ├── IV_175mK_20230918.txt
    │   ├── IV_180mK_20230919.txt
    │   ├── IV_185mK_20230919.txt
    │   ├── IV_190mK_20230919.txt
    │   ├── IV_195mK_20230919.txt
    │   ├── IV_200mK_20230919.txt
    │   ├── IV_205mK_20230919.txt
    │   ├── IV_210mK_20230919.txt
    │   ├── IV_215mK_20230919.txt
    │   ├── IV_220mK_20230919.txt
    │   ├── IV_225mK_20230919.txt
    │   ├── IV_230mK_20230919.txt
    │   └── IV_235mK_20230919.txt
    ├── 230918_results
    │   ├── IVtes_alpha.png
    │   ├── IVtes_contour.png
    │   ├── IVtes_fitting.png
    │   ├── IVtes_GTproperty.png
    │   ├── IVtes_IVproperty.png
    │   ├── IVtes_PRproperty.png
    │   ├── IVtes_result.html
    │   ├── IVtes_results.pdf
    │   ├── IVtes_results.tex
    │   └── IVtes_RTproperty.png
    ├── 230918_subset
    │   ├── IV_070mK_20230918.txt
    │   ├── IV_075mK_20230918.txt
    │   ├── IV_080mK_20230918.txt
    │   ├── IV_085mK_20230918.txt
    │   ├── IV_090mK_20230918.txt
    │   ├── IV_095mK_20230918.txt
    │   ├── IV_100mK_20230918.txt
    │   ├── IV_105mK_20230918.txt
    │   ├── IV_110mK_20230918.txt
    │   ├── IV_115mK_20230918.txt
    │   └── IV_120mK_20230918.txt
    ├── 230918_subset_results
    │   ├── IVtes_alpha.png
    │   ├── IVtes_contour.png
    │   ├── IVtes_fitting.png
    │   ├── IVtes_GTproperty.png
    │   ├── IVtes_IVproperty.png
    │   ├── IVtes_PRproperty.png
    │   ├── IVtes_result.html
    │   ├── IVtes_results.pdf
    │   ├── IVtes_results.tex
    │   └── IVtes_RTproperty.png
    ├── ivana_re.ipynb
    ├── ivana.ipynb
    ├── IVdata_A5
    │   ├── IV_110mK_20230331.txt
    │   ├── IV_115mK_20230331.txt
    │   ├── IV_120mK_20230331.txt
    │   ├── IV_125mK_20230331.txt
    │   ├── IV_130mK_20230331.txt
    │   ├── IV_135mK_20230331.txt
    │   ├── IV_140mK_20230331.txt
    │   ├── IV_145mK_20230331.txt
    │   └── IV_150mK_20230331.txt
    ├── IVdata_A5_results
    │   ├── IVtes_alpha.png
    │   ├── IVtes_contour.png
    │   ├── IVtes_fitting.png
    │   ├── IVtes_GTproperty.png
    │   ├── IVtes_IVproperty.png
    │   ├── IVtes_PRproperty.png
    │   ├── IVtes_result.html
    │   └── IVtes_RTproperty.png
    ├── IVtes_class.py
    ├── IVtes.py
    └── README.md

<p align="right"><a href="#section1"><strong>セクションに戻る»</strong></a></p>
</details>


## ディレクトリ用法
<div id="section2"></div>

&emsp;本ディレクトリ内のファイルおよび子ディレクトリを使うのに際して，**単に $I-V$ データの解析をしたい場合** と **$I-V$ データ解析の勉強をしたい場合** で必要なものが異なるため，自分の用途に応じて以下の通り読み進めてください。

- **単に $I-V$ データの解析をしたい**

    &emsp;解析したいデータが格納されたディレクトリ（別途）とIVtes.pyファイルをコピーしてください。次のセクションも併せて確認してください。
    - 簡単な注意点として，[解析データ](#解析データ) を読んでください。
    - 使い方の説明として，[IVtes.py](#ivtespy) は読んでください。
    ---
- **$I-V$ データ解析の勉強をしたい**
 
    &emsp;IVdata_A5ディレクトリとivana.ipynbファイルないしivana_re.ipynbファイルをコピーしてください。次のセクションも併せて確認してください。
    - 簡単な注意点として，[解析データ](#解析データ) を読んでください。
    - 使い方の説明として，[ivana.ipynb](#ivanaipynb) は読んでください。

</br>


### 解析データ
<div id="section2.1"></div>

&emsp;解析データの名前は，すべて <span style="color: #FF0000; ">"IV_温度(整数)mK_測定年月日.txt"</span> のように付けてください（例：IV_100mK_20250825.txt）。その上で，同じ素子で $I-V$ 測定して得られた全データをまとめて任意の名前のディレクトリに格納してください。作成したディレクトリは，解析スクリプトと同じ階層に置いてください。ただし，ivana.ipynbを用いて解析するのであればその限りではありません [^2] 。

> [!WARNING]
> &emsp;拡張子がtxt形式でない場合，txt形式に変えるか解析ファイルを一部修正する必要があります。

</br>


### ivana.ipynb
<div id="section2.2"></div>

#### ivana.ipynbファイル
&emsp;私が学部生のときに林さんから頂いて，勉強に使っていた解析ファイルです。[Google Colab](https://colab.research.google.com/drive/1rRyx1eIyx2i56KETObivz26Km9cdnt5h?authuser=1) 環境を前提とした構成になっています。当時の備忘録を [Qiita](https://qiita.com/Urchan/private/58672d874ca034e2a7a5) にまとめたので，こちらを参照してください。使い方としてここに載っている以上の説明はありませんが，大雑把には解析データが格納されたディレクトリを自分のGoogleドライブに保存して，あとは上のブロックから順に実行していくだけです。このとき，<span style="color: #FF0000; ">5ブロック目のglobのPATHは修正してください</span>。
> [!WARNING]
> &emsp;8ブロック目にあるSQUIDのパラメータは，JAXA120 Ea4 A5（IVdata_A5ディレクトリのデータ）用の値になっています。別の $I-V$ データを用いる際にはここの変更を忘れずに行ってください。

#### ivana_re.ipynbファイル
&emsp;ivana.ipynbファイルとは別に，ivana_re.ipynbファイルがあります。こちらは，Jupyter NotebookやVS Codeの環境で動かせるようにivana.ipynbファイルを微調整したファイルです。Google Colab環境が肌に合わない場合はこれを使ってください。同階層に解析データが格納されたディレクトリを置けば，あとはivana.ipynbファイルと同様に上から順にブロックを実行すれば使えます。また，**3ブロック目のglobと5ブロック目のSQUIDパラメータは，JAXA120 Sa4（230918ディレクトリのデータ）用になっています**。別のデータを使う際は修正を加えてください。
> [!IMPORTANT]
> &emsp;ivana.ipynbファイルとivana_re.ipynbファイルのいずれもobsoleteです。いくつかコードに不備が含まれる可能性があるので，[Qiita](https://qiita.com/Urchan/private/58672d874ca034e2a7a5) またはIVtes.pyを参考にして，調整が必要な箇所があれば修正してください。

</br>


### IVtes.py
<div id="section2.3"></div>

&emsp;ivana.ipynbにある解析内容をまとめて，一括に結果を出力できるようにしたスクリプトです。このファイルを解析したいデータが格納されているディレクトリと同階層に置き，その階層にてターミナルで

```bash:実行コマンド
python IVtes.py [解析データが格納されたディレクトリ名]/*.txt
```

を入力・実行すれば解析が進行します。このとき，同階層に新しく **\[解析データが格納されたディレクトリ名\]_results** という名前のディレクトリが自動で作られ，その中に解析結果がすべてpng形式で保存されます。また，結果をPDFでまとめて見るための $\TeX$ ファイルも同時に生成されます。"**IVtes_results.tex**" という名前で保存されるので，別途コンパイルをしてください。</br>
&emsp;さらに，このスクリプトファイルはいくつかオプションを用意してあります。そのうち主なものを以下にまとめます。

- `--exclude_files`：解析に使用しないファイルを指定します。デフォルトは未指定です。

    ```bash:使用例
    python IVtes.py 230918/*.txt --exclude_files IV_070mK_20230918.txt
    ```
- `--output_dir`：解析結果の出力先ディレクトリを指定します。デフォルトは "\[解析データが格納されたディレクトリ名\]_results" です。

    ```bash:使用例
    python IVtes.py 230918/*.txt --output_dir IV_230918
    ```
- `output_format`：解析結果をまとめる出力形式をPDFまたはHTMLで指定します。デフォルトはPDFで，IV_results.texというファイルが生成されます。HTMLを指定すると，IV_results.htmlというファイルが生成されます。

    ```bash:使用例
    python IVtes.py 230918/*.txt --output_format html
    ```

&emsp;その他のオプションについては，`python IVtes.py -h` のコマンドを入力・実行して，ヘルプメッセージで参照してください。

> [!NOTE]
> &emsp;出力される解析結果は以下の通りです。
> - **IVtes_IVproperty.png**：TESの $I-V$ 特性グラフ
> - **IVtes_PRproperty.png**：TESの $P-R$ 特性グラフ
> - **IVtes_RTproperty.png**：TESの $R-T$ 特性グラフ
> - **IVtes_GTproperty.png**：TESの $G-T$ 特性グラフ
> - **IVtes_alpha.png**：TESの $I-\alpha$ グラフ
> - **IVtes_fitting.png**：フィッティング結果のグラフ
> - **IVtes_contour.png**：フィッティングパラメータのコントア図
> - **IVtes_result.html** または **IVtes_result.pdf**：解析結果をまとめたレポートファイル
> 
> &emsp;それぞれの詳細については，ivana.ipynbファイルで解析する内容と同じのため，[Qiita](https://qiita.com/Urchan/private/58672d874ca034e2a7a5) を参照してください。

</br>


### IVtes_class.py
<div id="section2.4"></div>

&emsp;IVtes.pyファイルをオブジェクト指向で書き直したスクリプトファイルです。ChatGPTに通しただけで，一切修正を加えていないlegacyとなっています。このままでは動かず，まだ使用できません。

</br>


## 実行環境
<div id="section3"></div>

&emsp;IVtes.pyファイルを用いて解析する場合，当然ですがPythonの実行環境が必須です。このとき，必要なライブラリとしてはNumPy，Matplotlib，lmfitがあります。これらのインストールをしてから動作させてください。さらに，`--output_format pdf` オプションを使用してPDFレポートを生成する場合は，追加で $\LaTeX$ の実行環境と関連するパッケージが必要になります。これらの環境構築を事前に済ませてください。


<p align="right"><strong><a href="#top">トップへ»</a></strong></p>


[^1]: このREADMEファイルは，GitHub仕様に微修正しています。立教サーバの方を見るときは，VS Codeで開くことを推奨します。ただし，"Markdown All in One" および "Markdown Named CodeBlocks" の2つの拡張機能をインストールしてから使用してください。
[^2]: ivana_re.ipynbを用いて解析する場合は，3ブロック目にあるglobのPATHに注意すれば，同階層に置かなくても動きます。
