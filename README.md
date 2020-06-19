# calculate dynamics for quantum many body systems

# contents

* [1D transverse field Ising model](#1d-transverse-field-ising-model)
  * [after sudden quench, h=0 to small](#1DTFI-h-0-to-small)
  * [after sudden quench, h=inf to small](#1DTFI-h-inf-to-small)
  * [slow quench](#1DTFI-slow-quench)
* [2D transverse field Ising model](#2d-transverse-field-ising-model)
  * [after sudden quench, h=inf to small](#2DTFI-h-inf-to-small)
  * [slow quench](#2DTFI-slow-quench)

----

# 1D transverse field Ising model
## after sudden quench, h=0 to small <a name="1DTFI-h-0-to-small"></a>
* Usage
  ```console
  foo@bar:~$ python quench_dynamics_1d_FM_TFIsing.py -N [system size]
  ```
  Model: H(t) = - (\sum\_{i} \sigma\_i^z \sigma\_{i+1}^z + h(t) \sum\_{i} \sigma\_i^x) <br>
  Field: h(t=0)=0 (Ising GS: all up) --> h(t>0)=0.5 <br>
  The result for the system size N=16 is shown as an example.

* References
  * https://doi.org/10.1103/PhysRevA.2.1075
  * https://doi.org/10.1103/PhysRevA.3.786
  * https://doi.org/10.1103/PhysRevA.3.2137
  * https://doi.org/10.1103/PhysRevLett.106.227203
  * https://doi.org/10.1088/1742-5468/2012/07/P07016
  * https://doi.org/10.1088/1742-5468/2012/07/P07022
  * https://doi.org/10.1007/978-3-642-33039-1
  * https://amslaurea.unibo.it/id/eprint/15866 (I followed this notation.)
  * https://arxiv.org/abs/2005.03104 (numerical linked cluster expansion + forward propagation of pure states)
  * https://doi.org/10.1088/1742-5468/2016/06/064007 (long time limit)

* Results
  * Analitycal and asymptotic results in the thermodynamic limit are also given.

![magnetization (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing__field_0_to_small/fig_mz.png "magnetization (Ising direction)")

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing__field_0_to_small/fig_mx.png "magnetization (field direction)")

![logarithmic Loschmidt echo](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing__field_0_to_small/fig_loschmidt_echo.png "logarithmic Loschmidt echo")


## after sudden quench, h=inf to small <a name="1DTFI-h-inf-to-small"></a>
* Usage
  ```console
  foo@bar:~$ python quench_dynamics_1d_FM_TFIsing.py -N [system size]
  ```
  Model: H(t) = - (\sum\_{i} \sigma\_i^z \sigma\_{i+1}^z + h(t) \sum\_{i} \sigma\_i^x) <br>
  Field: h(t=0)=inf --> h(t>0)=0.1 <br>
  The result for the system size N=10 is shown as an example.

* References
  * https://doi.org/10.21468/SciPostPhys.4.2.013 (See Fig.2 for comparison, note that h'/J'=0.05 <--> h/J=0.1 where J'/4=J and h'/2=h)

* Results

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing__field_inf_to_small/fig_mx_vs_t.png "magnetization (field direction)")

![nearest neighbor spin correlation (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing__field_inf_to_small/fig_mz0mz1_vs_t.png "nearest neighbor spin correlation (Ising direction)")


## slow quench <a name="1DTFI-slow-quench"></a>
* Usage
  ```console
  foo@bar:~$ python slow_dynamics_1d_FM_TFIsing.py -N [system size] -tau [total time]
  ```
  Model: H(t) = - (\sum\_{i} \sigma\_i^z \sigma\_{i+1}^z + h(t) \sum\_{i} \sigma\_i^x) <br>
  Field: h(t)=hi+(hf-hi)\*t/tau (hi=2\*hc, hf=0, hc=1) <br>
  The result for the system size N=10 and the total time tau=32 is shown as an example.

* References
  * https://doi.org/10.1103/PhysRevLett.95.035701
  * https://doi.org/10.1103/PhysRevLett.95.105701
  * https://doi.org/10.1103/PhysRevLett.95.245701
  * https://doi.org/10.1103/PhysRevB.72.161201
  * https://doi.org/10.21468/SciPostPhys.1.1.003
  * https://doi.org/10.1080/00018732.2010.514702
  * https://doi.org/10.1007/978-3-642-33039-1
  * https://doi.org/10.1088/1367-2630/aa65bc
  * https://people.mines.edu/lcarr/wp-content/uploads/sites/23/2018/11/jaschke_thesis_2018.pdf
  * https://doi.org/10.1088/1742-5468/2016/06/064007 (long time limit)

* Results  
  * The kink density at the final time (t=tau) is expected to be scaled as (kink density)~1/(tau)^0.5 in the thermodynamic limit.

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_1d_FM_TFIsing__field_large_to_0/fig_mx.png "magnetization (field direction)")

![kink density](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_1d_FM_TFIsing__field_large_to_0/fig_kink_density.png "kink density")

![scaling of kink density 1](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_1d_FM_TFIsing__field_large_to_0/dat_kinkdens_scaling/fig_kinkdens_vs_inversetau.png "scaling of kink density 1")

![scaling of kink density 2](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_1d_FM_TFIsing__field_large_to_0/dat_kinkdens_scaling/fig_kinkdens_vs_tau_loglog.png "scaling of kink density 2")


# 2D transverse field Ising model
## after sudden quench, h=inf to small <a name="2DTFI-h-inf-to-small"></a>
* Usage
  ```console
  foo@bar:~$ python quench_dynamics_2d_FM_TFIsing.py -Lx [size Lx] -Ly [size Ly] -hi [initial field] -hf [final field] -tau [total time]
  ```
  Model: H(t) = - (\sum\_{i} \sigma\_i^z \sigma\_{i+1}^z + h(t) \sum\_{i} \sigma\_i^x) <br>
  Field: h(t=0)=inf --> h(t>0)=0.1, hc/10, hc, 2\*hc (hc=3.04438) <br>
  The result for the system size N=4x4 is shown as an example.

* References (h(t=0)=inf --> h(t>0)=0.1)
  * https://doi.org/10.21468/SciPostPhys.4.2.013 (See Fig.2 for comparison, note that h'/J'=0.05 <--> h/J=0.1 where J'/4=J and h'/2=h) (neural networks)

* References (h(t=0)=inf --> h(t>0)=hc/10, hc, 2\*hc (hc=3.04438))
  * https://doi.org/10.1103/PhysRevB.99.035115 (iPEPS)
  * https://arxiv.org/abs/1912.08828 (neural networks, comparison with iPEPS)
  * https://arxiv.org/abs/2005.03104 (numerical linked cluster expansion + forward propagation of pure states)

* Other references
  * https://doi.org/10.1038/srep38185 (VMC)
  * https://arxiv.org/abs/1710.07696 (numerical linked cluster expansion)
  * https://arxiv.org/abs/1811.09275 (MPS)
  * https://arxiv.org/abs/1910.10726 (MPS, triangular TFI)
  * https://arxiv.org/abs/1912.08831 (neural networks, time evolution error)
  * (https://doi.org/10.21468/SciPostPhys.6.3.031 (iPEPS, Heisenberg+field))
  * (https://doi.org/10.21468/SciPostPhys.8.2.021 (iPEPS, t-J))

* Results (h(t=0)=inf --> h(t>0)=0.1)

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_2d_FM_TFIsing__field_inf_to_small/Hf_0.1/fig_mx.png "magnetization (field direction)")

![nearest neighbor spin correlation (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_2d_FM_TFIsing__field_inf_to_small/Hf_0.1/fig_mz0mz1.png "nearest neighbor spin correlation (Ising direction)")

* Results (h(t=0)=inf --> h(t>0)=hc/10 (hc=3.04438))

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_2d_FM_TFIsing__field_inf_to_small/Hf_Hc_x_0.1/fig_mx.png "magnetization (field direction)")

![nearest neighbor spin correlation (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_2d_FM_TFIsing__field_inf_to_small/Hf_Hc_x_0.1/fig_mz0mz1.png "nearest neighbor spin correlation (Ising direction)")

* Results (h(t=0)=inf --> h(t>0)=hc (hc=3.04438))

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_2d_FM_TFIsing__field_inf_to_small/Hf_Hc_x_1/fig_mx.png "magnetization (field direction)")

![nearest neighbor spin correlation (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_2d_FM_TFIsing__field_inf_to_small/Hf_Hc_x_1/fig_mz0mz1.png "nearest neighbor spin correlation (Ising direction)")

* Results (h(t=0)=inf --> h(t>0)=2\*hc (hc=3.04438))

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_2d_FM_TFIsing__field_inf_to_small/Hf_Hc_x_2/fig_mx.png "magnetization (field direction)")

![nearest neighbor spin correlation (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_2d_FM_TFIsing__field_inf_to_small/Hf_Hc_x_2/fig_mz0mz1.png "nearest neighbor spin correlation (Ising direction)")

## slow quench <a name="2DTFI-slow-quench"></a>
* Usage
  ```console
  foo@bar:~$ python slow_dynamics_2d_FM_TFIsing.py -Lx [size Lx] -Ly [size Ly] -hi [initial field] -hf [final field] -tau [total time]
  ```
  Model: H(t) = - (\sum\_{i} \sigma\_i^z \sigma\_{i+1}^z + h(t) \sum\_{i} \sigma\_i^x) <br>
  Field: h(t)=hi+(hf-hi)\*t/tau (hi=2\*hc, hf=0, hc=3.04438) <br>
  The result for the system size N=4x4 and the total time tau=32 is shown as an example.

* Related references
  * Kibble-Zurek
    * https://doi.org/10.1088/0305-4470/9/8/029
    * https://doi.org/10.1016/0370-1573(80)90091-5
    * https://doi.org/10.1038/317505a0
    * https://www.actaphys.uj.edu.pl/R/24/7/1301 https://www.actaphys.uj.edu.pl/R/24/7/1301/pdf
    * https://doi.org/10.1016/S0370-1573(96)00009-9
  * dynamics by QMC (real time dynamics ~ imaginary time dynamics when v=1/tau << 1)
    * http://indico.ictp.it/event/a14246/session/17/contribution/28/material/6/0.pdf http://physics.bu.edu/~sandvik/trieste15/sandvik4.pdf
    * https://doi.org/10.1103/PhysRevB.84.224303
    * https://doi.org/10.1088/0953-8984/25/40/404216
    * https://doi.org/10.1103/PhysRevB.87.174302
  * QA
    * https://arxiv.org/abs/2001.11637 (scaling in D-dimensions: kink density, higher order cumulants: see references therein)
    * http://hdl.handle.net/2433/142655
    * http://hdl.handle.net/2433/110379 (scaling of residual energy)
    * https://www.smapip.is.tohoku.ac.jp/~dex-smi/2006/Workshop200612/ExtendedAbstracts/SeiSuzuki.pdf

* Results

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_2d_FM_TFIsing__field_large_to_0/fig_mx.png "magnetization (field direction)")

![nearest neighbor spin correlation (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_2d_FM_TFIsing__field_large_to_0/fig_mz0mz1.png "nearest neighbor spin correlation (Ising direction)")
