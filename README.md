# 1DLax
Hydrodynamics calculation code with 1D Lax scheme using MPI.

mpirunで実行すると各コア毎にファイルが生成されます。
シェルスクリプトなどでそれらのファイルを結合し、gnuplotで表示させると、
ある時刻での計算領域全体の様子を確認できます。

今度、簡単なシェルスクリプトも作ろうかと思います。
