# calculate dynamics for quantum many body systems

* 1D transverse field Ising model (after sudden quench)
  * Usage
    ```console
    foo@bar:~$ python quench_dynamics_1d_FM_TFIsing.py -N [system size]
    ```
    Model: H(t) = - (\sum_{i} \sigma_i^z \sigma_{i+1}^z + h(t) \sum_{i} \sigma_i^x) <br>
    Field: h(t=0)=0 (Ising GS: all up) --> h(t>0)=0.5 <br>
    The result for the system size N=16 is shown as an example.
  * Results
    * Analitycal and asymptotic results in the thermodynamic limit are also given.

![magnetization (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing__field_0_to_small/fig_mz.png "magnetization (Ising direction)")

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing__field_0_to_small/fig_mx.png "magnetization (field direction)")

![logarithmic Loschmidt echo](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing__field_0_to_small/fig_loschmidt_echo.png "logarithmic Loschmidt echo")

  * References
    * https://doi.org/10.1103/PhysRevLett.106.227203
    * https://doi.org/10.1088/1742-5468/2012/07/P07016
    * https://doi.org/10.1088/1742-5468/2012/07/P07022
    * https://doi.org/10.1007/978-3-642-33039-1
    * https://amslaurea.unibo.it/id/eprint/15866 (I followed this notation.)

* 1D transverse field Ising model (slow quench)
  * Usage
    ```console
    foo@bar:~$ python slow_dynamics_1d_FM_TFIsing.py -N [system size] -tau [total time]
    ```
    Model: H(t) = - (\sum_{i} \sigma_i^z \sigma_{i+1}^z + h(t) \sum_{i} \sigma_i^x) <br>
    Field: h(t)=hi+(hf-hi)\*t/tau (hi=2\*hc, hf=0, hc=1) <br>
    The result for the system size N=10 and the total time tau=32 is shown as an example.
  * Results

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_1d_FM_TFIsing__field_large_to_0/fig_mx.png "magnetization (field direction)")

![kink density](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_1d_FM_TFIsing__field_large_to_0/fig_kink_density.png "kink density")

![scaling of kink density 1](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_1d_FM_TFIsing__field_large_to_0/dat_kinkdens_scaling/fig_kinkdens_vs_inversetau.png "scaling of kink density 1")

![scaling of kink density 2](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/slow_dynamics_1d_FM_TFIsing__field_large_to_0/dat_kinkdens_scaling/fig_kinkdens_vs_tau_loglog.png "scaling of kink density 2") <br>
The kink density is expected to be scaled as (kink density)~1/(tau)^0.5 in the thermodynamic limit.

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
