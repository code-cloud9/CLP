## Introduction

The code is for reproducing the numerical experiments and data example results in 'Conformal Link Prediction with False Discovery Rate Control'.

## Guidelines for Result Replication

### Simulation Results in the Main File

1) For tunning clp, please go to folder `./simulation/tune clp/` and run `clp_c*_stgs.m` where `*` is 05,1,15 for `c` values in {0.5,1,1.5}. The results are saved in the folder `./results/`.
2) For comparison in weighted networks, please go to folder `./simulation/comparison/weighted network/` and run `*_unifc_stgs.m` where `*` is clp or cmc. The results are saved in the folder `./results/`.
3) For comparison in unweighted networks, please go to folder `./simulation/comparison/unweighted network/` and run `*_unifc_stgs.m` where `*` is clp or cmc or dss. The results are saved in the folder `./results/`.

### Numerical Analysis of Trading Network in the Main File
  
Please go to folder `./application/` and run `trade_main.m` for different significance levels.

### Simulation Results in the Appendix

1) For comparison in weighted networks, please go to folder `./simulation/appendix/weighted network/` and run `*_signal_stgs.m` where `*` is clp or cmc. The results are saved in the folder `./results/`.
2) For comparison in unweighted networks, please go to folder `./simulation/appendix/unweighted network/` and run `*_signal_stgs.m` where `*` is clp or cmc or dss. The results are saved in the folder `./results/`.

### Numerical Analysis of Trading Network in the Appendix
  
Please go to folder `./application/` and run `trade_appendix.m` for different null values.

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

