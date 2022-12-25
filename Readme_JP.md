# Costas Loop Search

# 使い方

最初に受信したい信号のIQ信号を記録したファイルを準備します．
受信時間はおおむね１秒としてください．(長いと表示に時間がかかります)
下図はGNURadioでIQファイルを作成する一例です．

![IQ_record](/image/IQ_record.png)

ターミナルで起動します．引数でIQファイル名を指定します．
~~~
$ python costas_loop_search.py iqfile_name
~~~
信号の確認のためのスペクトラムの表示のあとにシンボル同期を行ったコンスタレーションを表示します．

![PSD](/image/psd.png)
![SymbolSync](/image/after_symbolsync.png)

円形になっていますが周波数が微妙にずれているためです．
その後にコスタスループによるキャリアリカバリが行われます．

alphaとbetaの２つのパラメータがありますがそれぞれ５通りに変化させて計算して結果を表示します．しばらく時間がかかるでしょう．

最初のプロットは周波数オフセットです．同期するまでの速さと定常時のジッターの大きさを見ます．
見やすくするために最初の３０００シンボルについて表示しています．

![freq_log](/image/freq_log.png)

次のプロットは同様にパラメータを変化させてコンスタレーションをプロットします．
ここではBPSKなので２つの点の集まりが現れます．
見やすくするためにロックインしていない最初の方のシンボルを除いて表示します．

２つの点の集まりがはっきり分かれているか見ます．

![constellation](/image/constellation.png)

次に２つのプロットからもっとも良かったパラメータを選んで入力します．

また同期するまでのシンボル数を整数で入力します．
目分量です．良い方法を知っている人がいたら教えてください．

入力したパラメータによって復調された周波数オフセットとコンスタレーションが表示されます．

![freq_log2](/image/freq_log2.png)

![constellation2](/image/constellation2.png)

もっと別なパラメータで探したい場合はプログラムの中にalpha とbetaの候補のリストがあります．
これを書き換えれば新しいパラメータでプロットします．

実際の受信機でパラメータを決定する場合はいくつかのIQファイルで確認することが必要と思います．
周波数オフセットは変化したりばらついたりするでしょうから．

# 背景

私はSDRやPythonでBPSK受信機を作っていますがコスタスループの良いパラメータの求め方が分かりませんでした．
そこで色々なパラメータについてロックインタイムやコンスタレーションをプロットして良さそうな値を選ぶことにしました．

PySDRのドキュメント（https://pysdr.org/content/rds.html　）　はたいへん役にたちました．感謝します．