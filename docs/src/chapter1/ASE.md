# Atomic Simulation Environment (ASE)

## ASEとは？
Atomic Simulation Environment (ASE)というパッケージがあります。これは、第一原理計算や分子動力学法などの様々なパッケージをPythonによって統一的に操作できるフリーのソフトウエアです。

https://wiki.fysik.dtu.dk/ase/index.html

面白いのは、ある原子配置が与えられた時に、その時のエネルギーと原子間に働く力さえ返せる外部パッケージであれば、様々なものをPythonで扱うことができるようになります。有名どころでは、第一原理計算パッケージのVASPやQuantumEspressoなどがあります。

また、分子動力学法ではなく、普通の第一原理計算ももちろん可能です。
ここでは、MateriApps Live!上でASEをインストールする方法について述べます。

## いくつかのやり方
ASEを使うにはいくつかのやり方があります。

- 一つは、ターミナル（左下のツバメのボタンを押して「System Tools」から「LXTerminal」を選ぶ）から
- Jupyter Notebookの使用（MatLabやMathematica的な操作感）
- Virtual BoxにホストOSからSSH接続し、使う

などがあります。Jupyter notebookが一番簡単かもしれないので、Jupyter notebookを使ってやってみましょう。

[Jupyter notebookの使用](@ref)

以下はターミナルを用いたやり方です。あとでJupyter Notebookを使ったやり方に修正します。


##  ASEのインストール
MateriApps LIVE!にはPythonが入っていますので、

```sh
pip install ase
```
とすればASEを簡単に入れることができます。このPythonのバージョンは2.7です。
ここで注意ですが、以下のサンプルコードはPython3系では動きませんでした。Python2系では動きました。

Python3で動かすためには、site-package/are/espresso.pyの292行目のmapをlist()で囲めば良いようです。

## サンプルコードの実行

NaCl結晶の格子定数を求めてみましょう。
実験的な格子定数は
http://crystalbase.co.jp/index/b/nacl.html
で、関連する日本語論文は
https://www.jstage.jst.go.jp/article/jsssj/28/3/28_3_144/_pdf
にあります。


https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html
にあるサンプルコードはそのままでは動きません。
```rocksalt.set_calculator(calc)```
が足りません。また、擬ポテンシャルをちゃんと持っているものを指定しなければなりません。
ですので、

```python
from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.constraints import UnitCellFilter
from ase.optimize import LBFGS


pseudopotentials = {'Na': 'Na.pbe-spn-kjpaw_psl.1.0.0.UPF',
                    'Cl': 'Cl.pbe-n-rrkjus_psl.1.0.0.UPF'}

rocksalt = bulk('NaCl', crystalstructure='rocksalt', a=6.0)
calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(3, 3, 3))

rocksalt.set_calculator(calc)

ucf = UnitCellFilter(rocksalt)
opt = LBFGS(ucf)
opt.run(fmax=0.005)

# cubic lattic constant
print((8*rocksalt.get_volume()/len(rocksalt))**(1.0/3.0))

```
とします。

ここで、MateriApps LIVE!でのQuantum Espressoの擬ポテンシャルは
```/usr/share/espresso/pseudo```
にありますので、

```bash
cd /usr/share/espresso/pseudo
sudo wget https://www.quantum-espresso.org/upf_files/Na.pbe-spn-kjpaw_psl.1.0.0.UPF
sudo wget https://www.quantum-espresso.org/upf_files/Cl.pbe-n-rrkjus_psl.1.0.0.UPF
```
として、NaとClの擬ポテンシャルをダウンロードします。なお、他の擬ポテンシャルは
https://www.quantum-espresso.org/pseudopotentials
にあります。

コードでは、格子定数を最初は6として、設定したあと、構造最適化をしています。最後に、構造最適化して出てきた格子定数をプリントしています。
サンプルコードをファイルに保存して、実行すると、

```console
user@malive:~$ python test.py 
       Step     Time          Energy         fmax
LBFGS:    0 06:20:31    -1960.772207        0.5120
LBFGS:    1 06:20:36    -1960.780670        0.4990
LBFGS:    2 06:20:41    -1960.819410        0.3478
LBFGS:    3 06:20:46    -1960.836529        0.1518
LBFGS:    4 06:20:50    -1960.831494        0.0331
LBFGS:    5 06:20:55    -1960.834038        0.0063
LBFGS:    6 06:21:00    -1960.834046        0.0009
5.659067832564933
```
こんな感じで出力されて、第一原理計算を用いてNaClの結晶構造の最適化を行うことができました。

## Julia言語での実装
PythonよりもJuliaが好きなので、Juliaでの実装についても述べておきます。

最後に、Julia言語での実装も書いておきます。
MateriApps LIVE!のディスクイメージへのインストールは、

```sh
wget https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.0-linux-x86_64.tar.gz
tar xzvf julia-1.1.0-linux-x86_64.tar.gz 
```
となります。ここではバージョンは1.1.0ですが、現在(2020年3月28日現在、1.4.0が出ていますので、そちらで構いません)

```sh
~/julia-1.1.0/bin/julia
```
で立ち上げて、さっきまで使っていたPythonを使うように

```sh
ENV["PYTHON"]="/usr/bin/python"
```
とし、
]を押してパッケージモードにして、

```
add PyCall
```
としてPyCallをインストールします。

そして、コードは、

```julia
using PyCall
const bulk = pyimport("ase.build").bulk
const Espresso = pyimport("ase.calculators.espresso").Espresso
const UnitCellFilter = pyimport("ase.constraints").UnitCellFilter
const LBFGS = pyimport("ase.optimize").LBFGS


pseudopotentials =Dict("Na"=> "Na.pbe-spn-kjpaw_psl.1.0.0.UPF","Cl" => "Cl.pbe-n-rrkjus_psl.1.0.0.UPF")

rocksalt = bulk("NaCl", crystalstructure="rocksalt", a=6.0)
calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=true, tprnfor=true, kpts=(3, 3, 3))

rocksalt.set_calculator(calc)


ucf = UnitCellFilter(rocksalt)
opt = LBFGS(ucf)
opt.run(fmax=0.005)

# cubic lattic constant
println((8*rocksalt.get_volume()/length(rocksalt))^(1.0/3.0))
```
となります。



