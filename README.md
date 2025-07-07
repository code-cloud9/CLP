## Introduction

The code is for reproducing the numerical experiments and data analysis results in 'Conformal Link Prediction with False Discovery Rate Control'.

## Guidelines for Result Reproduction

### Simulation Results

<!-- For tunning clp -->
1) **Figure 3 for Simulation 1:** Please go to folder `./simulation/tune clp/` and run `clp_c*_stgs.m`, where `*` is 05,1,15 for `c` values in {0.5,1,1.5}, to obtain the results. To save time, you may use the intermediate results saved in the folder `./results/` to calculate the empirical power and FDR.

<!-- For comparison in weighted networks -->
2) **Figure 4 for Simulation 2:** Please go to folder `./simulation/comparison/weighted network/` and run `*_unifc_stgs.m`, where `*` is `clp` or `cmc`. To save time, you may use the intermediate results saved in the folder `./results/` to calculate the empirical power and FDR.

<!-- For comparison in unweighted networks -->
3) **Figure 5 for Simulation 2:** Please go to folder `./simulation/comparison/unweighted network/` and run `*_unifc_stgs.m`, where `*` is `clp` or `cmc` or `dss`. To save time, you may use the intermediate results saved in the folder `./results/` to calculate the empirical power and FDR.

### Numerical Analysis of Trading Network
  
1) **Figure 6:** Please go to folder `./application/` and run `trade_main.m` for different significance levels.

### Additional Results in the Appendix

<!-- For comparison in weighted networks -->
1) **Figure B.1:** Please go to folder `./simulation/appendix/weighted network/` and run `*_signal_stgs.m`, where `*` is `clp` or `cmc`, to obtain the results. To save time, you may use the intermediate results saved in the folder `./results/` to calculate the empirical power and FDR.

<!-- For comparison in unweighted networks -->
2) **Figure B.2:** Please go to folder `./simulation/appendix/unweighted network/` and run `*_signal_stgs.m`, where `*` is `clp` or `cmc` or `dss`. To save time, you may use the intermediate results saved in the folder `./results/` to calculate the empirical power and FDR.

<!-- Numerical Analysis of Trading Network in the Appendix -->
3) **Figure C.1:** Please go to folder `./application/` and run `trade_appendix.m` for different null values.

## License

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

