Code for joint estimation of velocity, angle and range in multipath environment (JEVAR).

The detailed algorithm can be seen in the paper --- Z. Yang, R. Wang, Y. Jiang and J. Li, "Joint Estimation of Velocity, Angle-of-Arrival and Range (JEVAR) Using a Conjugate Pair of Zadoff-Chu Sequences," in IEEE Transactions on Signal Processing, vol. 69, pp. 6009-6022, 2021, doi: 10.1109/TSP.2021.3122907.



Introductions for functions



main_compare: compare the AP with the SAGE in SNR vs. RMSE

main_JEVAR: using the AP method to solve the JEVAR problem



gen_multipathSig: generate multipath signal

raised_cosine: generate raised cosine filter

JEVAR_ap: AP method

JEVAR_sage: SAGE method

est_single, est_single_sage, est_projection: estimate single-path parameters

coarse_est_2Dfft: use 2D-fft to obtain the initial estimate of frequency offset, DOA and time delay

fine_est_Newton, fine_est_Newton_projection: use the Newton's iteration method to obtain the refined estimate

crb_cal: calculate the CRBs
