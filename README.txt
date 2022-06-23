0. 説明

check_validity は巡回路ファイルが順序制約を守っている巡回路を表現
しているファイルであるかどうかをチェックするプログラム. 

1. コンパイル

$ gcc check_validity.c -o check_validity.out -lm

2. 実行の仕方

$ ./check_validity.out <dat filename> <tsp filename>

3. 実行例

$ ./solve.out in_data2/d198.tsp
34965.6
$ ./check_validity.out tour1.dat in_data2/d198.tsp
Valid dat file. 34965.607960




