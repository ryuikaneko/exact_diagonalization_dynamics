# calculate dynamics for quantum many body systems

* 1D transverse field Ising model (after sudden quench)
  * Usage
    ```console
    foo@bar:~$ python quench_dynamics_1d_FM_TFIsing.py -N [system size]
    ```
    The result for the system size N=16 is shown as an example.
  * Results
    * Analitycal and asymptotic results are also given.

![magnetization (Ising direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing/fig_mz.png "magnetization (Ising direction)")

![magnetization (field direction)](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing/fig_mx.png "magnetization (field direction)")

![logarithmic Loschmidt echo](https://raw.githubusercontent.com/ryuikaneko/exact_diagonalization_dynamics/master/quench_dynamics_1d_FM_TFIsing/fig_loschmidt_echo.png "logarithmic Loschmidt echo")

  * References
    * https://doi.org/10.1103/PhysRevLett.106.227203
    * https://doi.org/10.1088/1742-5468/2012/07/P07016
    * https://doi.org/10.1088/1742-5468/2012/07/P07022
    * https://doi.org/10.1007/978-3-642-33039-1
    * https://amslaurea.unibo.it/id/eprint/15866 (I followed this notation.)
