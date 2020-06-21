# RubyGROWTH

## はじめに
RubyGROWTHは冬季雷の高エネルギー現象を観測するGamma-Ray Observation of Winter Thunderclouds (GROWTH) 実験において、定型的なデータ解析を行うために用意されたRubyのライブラリです。

本ライブラリを使用するには[RubyFITS](https://github.com/yuasatakayuki/RubyFits)と[RubyROOT](https://github.com/odakahirokazu/RubyROOT)のインストールが必要です。これらのインストール方法については各ページを参照していただくか、Macの場合は[こちら](https://blog.yuuki-wd.space/archives/17)もご参照ください。なおCERNの開発するROOTに沿ったライブラリとなっているため、必要に応じてROOTの各種ドキュメントも参照してください。

本ライブラリのインストールは任意のディレクトリで

`git clone https://github.com/YuukiWada/RubyGROWTH.git`

として、~/.bashrcなどに

`export RUBYLIB=/path-to-ripository/RubyGROWTH/lib:$RUBYLIB`

と記述してください。

本ライブラリを使用して発生した不利益・損害などについて、作者は一切の責任を負いません。

## 使い方

このライブラリで出来ることは

- スペクトルのquick look (ql_spec)
- ライトカーブのquicklook (ql_lightcurve)
- 任意の時間帯のライトカーブの抽出 (lightcuirve, lc_unixtime)
- 任意の時間帯のスペクトルの抽出 (spec)
- バックグラウンドを引いたスペクトルの抽出 (spec_src)
- 波形の最大値と最小値の散布図の抽出 (scatter)
- 任意の時間におけるΔt分布の作成 (delta_time)
- 文字列からunix timeへの変換 (unixtime) <- おまけです

です。カッコの中はメソッド名です。

### FITSファイルの読み込み

インスタンスを作成するには

    require "GROWTH"
    growth=Growth.new(fits_file)

とします。以下、インスタンス名はgrowthとしておきます。fits_fileには読み込みたいFITSファイルのパスを文字列で、もしくは文字列の配列で与えます。配列で与えた場合は、複数のFITSファイルを読み込ませることが可能です。ただし復数のファイルが連続したものでない場合 (ファイル間で観測時間のギャップがある) や、正しい順番で並んでいない場合は警告が出ます。その場合は意図したファイルが入力されているか確認してください。

読み込ませるFITSファイルはLv.2のパイプライン処理後のもののみを受け付け、時刻・エネルギー校正がされていないと正常に読み込まれません。

### unix timeの取得

FITSファイルに格納されている光子イベントの時刻はunix timeで管理されています。文字列からunix timeを取得するには

    unix_time=growth.unixtime("2017-12-10 10:19:10", "Asia/Tokyo")

とします。返り値の型はfloatです。第1引数はRubyのTimeクラスで読める形式であればよく、比較的自由なフォーマットに対応しています。第2引数は省略できますが、その場合はデフォルトでAsia/Tokyoが選択されます。第2引数のフォーマットは[こちら](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones)のTZ database nameです。代表的なものをあげておくと

- Asia/Tokyo
- Asia/Seoul
- Asia/Shanghai
- Asia/Singapore
- Asia/Jakarta
- Asia/Dubai
- Europe/Moscow
- Europe/Athens
- Europe/Paris
- Europe/London
- America/New_York
- America/Chicago
- America/Denver
- America/Los_Angeles
- Pacific/Honolulu
- Etc/UTC

です。

### エネルギースペクトルのquick look

ファイルに格納されているスペクトルをとりあえず見てみたいという場合に便利です。GROWTH実験では30分で1ファイルとなっていますので、30分の積算スペクトルを見ることができます。

    spec=growth.ql_spec(adc, min, max, num, log, name)

で取得できます。adcはADCボードの入力チャンネル (0,1,2,3: int型)、minとmaxはヒストグラムの最小値と最大値 (MeV単位、float型)、numはヒストグラムの分割数 (int型)、logはbool型でtrueだとlog表示で等間隔になるように、falseだとlinear表示で等間隔になるようにビンが分割されます。nameはstringでヒストグラムのオブジェクト名を与えます。1つだけのヒストグラムを作成する場合は何でも良いですが、2つ以上を作成する場合は、同じにならないように気をつけてください。これはROOTがC++で動いているため、設定が必要です。返り値はROOTのTH1Dオブジェクトです。

quick lookは1つのファイルのみを読み込むことを前提としているため、インスタンスを作成する際に復数のFITSファイルを読み込ませていると、動作しません。

### ライトカーブのquick look

スペクトルと同様にライトカーブもquick lookできます。同様に1つのFITSファイルのみを読み込ませてください。

    lc=growth.ql_lightcurve(adc, min, max, width, name)

adc、min、max、nameはql_specと同じです。widthはライトカーブのビン幅 (秒、float型) です。X軸は時刻で表示されます。返り値はROOTのTH1Dオブジェクトです。

### ライトカーブの抽出

時間を指定してライトカーブを抽出します。

    lc=growth.lightcurve(adc, center, before, after, min, max, width, name)

centerは抽出した時間帯の中心をunix timeで与えます。unixtimeメソッドで取得すると便利です。beforeとafterはcenterから前後に抽出する時間です (秒、float型、いずれも正の値)。このメソッドではデッドタイムが補正されます。ただしデッドタイム補正のアルゴリズム上、あまりに短すぎるビン幅や狭すぎるエネルギー帯域を指定した場合には、正常な補正ができない可能性があります。

返り値はTH1Fです。またX軸はcenterで与えた時刻からの秒数となります。

### ライトカーブの抽出 (X軸が時刻)

ライトカーブの横軸をunix timeで出力したい場合は、こちらのメソッドを利用します。

    lc=growth.lc_unixtime(adc, center, before, after, min, max, width, name)

引数や返り値の仕様はlightcurveと同じです。

### スペクトルの抽出

時刻を指定してスペクトルを抽出します。

    spec=growth.spectrum(adc, start, duration, min, max, num, log, name)

startにはスペクトルを抽出する時間領域の始点のunix timeをfloat型で、durationは抽出する時間幅をfloat型・秒数で与えます。その他の引数はql_specと同じです。返り値はTH1Dです。デッドタイム補正が行われます。

なお抽出する時間帯が読み込ませたFITSファイルに全て含まれているかどうか、注意してください。時間帯が見切れている場合、スペクトルのnormalizationを間違えることになります。それを避けるため、前後のFITSファイルを読み込ませ、時間幅に余裕をもたせることをおすすめします。

### バックグラウンドを引いたスペクトルの抽出

バックグラウンドを引いたスペクトルを抽出します。

    spec=growth.spec_src(adc, start, duration, bgd_start, bgd_duration, min, max, num, log, name)

bgd_startにはバックグラウンドを抽出する時間領域の始点をunix timeで、bgd_durationには抽出する時間幅をfloat型・秒数で与えます。その他はspectrumと同じです。返り値はTH1Dです。デッドタイム補正が行われます。

### 波形の最小値と最大値の散布図の抽出

主にショートバーストの解析向けですが、アンダーシュートを確認するための散布図を抽出します。

  scat=growth.scatter(adc=0, center, before, after)

centerは抽出したい時間の中心のunix time、beforeとafterは抽出したい前後の秒数です。返り値は配列で最初の要素には最大値のTGraph、2番目の要素には最小値のTGraphが格納されています。3・4・5番目の配列にはTGraphを構成する値の配列が格納されており、3番目は時刻、4番目と5番目はそれぞれ波形の最大値と最小値の配列です。

### Δt分布の取得

時間を指定してΔt分布を取得します。

    dt=growth.delta_time(adc, start, duration, min, max, num, name)

minとmaxはヒストグラムの最小値と最大値を与えます (float型、秒)。numはヒストグラムのビン数です。その他の引数はspectrumと同じです。返り値はTH1Dです。

## サンプルプログラム

取得したTH1F/TH1Dオブジェクトを表示する、ROOTファイルに書き出す、QDP形式でテキスト書き出しするサンプルプログラムはそれぞれ

- samples/hist_draw.rb
- samples/hist_root.rb
- samplea/hist_qdp.rb

にあります。また散布図の表示、ROOTファイルへの書き出し、CSVでの書き出しのサンプルプログラムは

- samples/scatter_draw.rb
- samples/scatter_root.rb
- samples/scatter_csv.rb

にあります。
