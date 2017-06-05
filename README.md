# DBMT
Dynamic Bayesian Multitaper Estimation MATLAB Codes

Description: This repository contains implementations of the algorithms developed in Dynamic Bayesian Multi Taper estimation paradigm. 

Copyright (c) 2017 Proloy Das All Rights Reserved 

Contact: proloy@umd.edu
(Das and Babadi, Dynamic Bayesian Multitaper Spectral Analysis,(under review) Submitted @ IEEE Trans. Signal Process. on 06/05/2017)

Date: June 5, 2017

Requirements:
  implemented in Matlab R2016b version, but should run on most versions.
  
Contents:
  1. main.m: Master script.
  2. TSpectrogram.m: genrates single taper sSpectrogram estimates.
  2. MTSpectrogram.m: genrates overlapped multi taper spectrogram estimates.
  3. DBMTSpectrogram.m: genrates DBMT estimates.
  4. log_DBMTSpectrogram.m: genrates log_DBMT estimates.
  5. DBMT_EM.m: EM step for DBMT algorithm.
  6. log_DBMT_EM.m: EM step for log-DBMT algorithm.
  7. Window_then_Taper.m: Segments the data and then multiplyes by taper (Used in DBMTSpectrogram).
  8. Post_mode_var.m: Calculates posterior mode and variance in log_DBMT algorithm.
  9. confidence_interval_CMT.m: plots 95% confidence interval of conventional multitaper estimates of a single window.
  10. confidence_interval_DBMT.m: plots 95% confidence interval for DBMT estimate of a single window.
  11. confidence_interval_log_DBMT.m: plots 95% confidence interval for log_DBMT estimate of a single window.


Instructions: Simple and easy.
  Download all the codes in a directory and run main.m, that will generate one toy_example, its three spectrogram estimates as well as
  the true estimate. It also compares estimates of a single window, while constructing the confidence intervals. To use 
  the functions individually please look at the function descriptions.
  
Additional Note: Make sure that 3 folders (named CMT, DBMT, log_DBMT) is created after running the main script, they are necessary for 
  constructing confidence intervals.
