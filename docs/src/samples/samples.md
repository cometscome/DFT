# ASEによる計算例

## 銅への窒素の吸着エネルギー

https://wiki.fysik.dtu.dk/ase/tutorials/surface.html
窒素分子が銅のスラブへ吸着する際のエネルギーを見積もってみましょう。このエネルギーは、「独立した窒素分子のエネルギー＋銅原子スラブのエネルギー」と「銅原子スラブに窒素分子が吸着しているエネルギー」の差となります。

コードは以下に示します。適宜コメントを入れておきました。

```python
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate

h = 1.85
d = 1.10

slab = fcc111('Cu', size=(4, 4, 2), vacuum=10.0) #銅原子スラブのセット

slab.set_calculator(EMT()) #銅原子スラブの計算にはEMTを使用
e_slab = slab.get_potential_energy() #スラブのポテンシャルエネルギーを計算

molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)]) #窒素分子のセット。(0,0,0)が一つ目のNの位置、(0,0,d)が二つ目のNの位置。
molecule.set_calculator(EMT()) #窒素分子の計算にはEMTを使用
e_N2 = molecule.get_potential_energy() #窒素分子のポテンシャルエネルギーの計算

add_adsorbate(slab, molecule, h, 'ontop') #窒素分子を上にのせる
constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab]) #拘束条件としては、計算を高速化するため、銅原子の位置を緩和させずに固定
slab.set_constraint(constraint) #拘束条件をセット
dyn = QuasiNewton(slab, trajectory='N2Cu.traj') #準ニュートン法を設定
dyn.run(fmax=0.05) #構造緩和スタート。全ての原子に働く力がfmax以下になるまで。

print('Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy()) #ポテンシャルエネルギーを計算し、先ほどの二つとの差を取る
```
となります。

```slab.set_calculator(EMT()) ```の部分をQuantum Espressoに置き換えれば第一原理計算で見積もることができます。

以下工事中。