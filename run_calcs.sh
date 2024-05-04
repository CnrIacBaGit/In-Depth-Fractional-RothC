#!/bin/sh
ulimit -c unlimited

rm -f -d -r results
mkdir results
mkdir results/3_9_1year2
mkdir results/3_10_1year2
mkdir results/5_9_1year2
mkdir results/3_9_1year2_noV

mkdir results/3_9
mkdir results/3_9_099_0
mkdir results/3_9_099_1
mkdir results/3_9_restore

time ./z_rothC NB 100 testing 3 Zoutstep 1 C_q 1.0 corrector_type 0 Tm 365 Om 1 LB 1 H0 -3 ls_eps 1e-10 ls_max_iter 200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth1.txt Ta_file Ta1.txt SoC_file SOC_inputs2.txt Kc_file Kc1.txt C00 2.78 C01 2.78 C02 2.78 C03 2.78
./move_to_folder.sh results/3_9_1year2
time ./z_rothC NB 100 testing 3 Zoutstep 1 C_q 1.0 corrector_type 0 Tm 365 Om 1 LB 1 H0 -3 ls_eps 1e-10 ls_max_iter 200 vgm vgm2.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth1.txt Ta_file Ta1.txt SoC_file SOC_inputs2.txt Kc_file Kc1.txt C00 2.78 C01 2.78 C02 2.78 C03 2.78
./move_to_folder.sh results/3_10_1year2
time ./z_rothC NB 100 testing 3 Zoutstep 1 C_q 1.0 corrector_type 0 Tm 365 Om 1 LB 1 H0 -3 ls_eps 1e-10 ls_max_iter 200 vgm vgm3.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth1.txt Ta_file Ta1.txt SoC_file SOC_inputs2.txt Kc_file Kc1.txt C00 2.78 C01 2.78 C02 2.78 C03 2.78
./move_to_folder.sh results/5_9_1year2
time ./z_rothC NB 100 Zoutstep 1 testing 2 C_q 1.0 corrector_type 0 Tm 365 Om 1 LB 1 H0 -3 ls_eps 1e-10 ls_max_iter 200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth1.txt Ta_file Ta1.txt SoC_file SOC_inputs2.txt Kc_file Kc1.txt C00 2.78 C01 2.78 C02 2.78 C03 2.78
./move_to_folder.sh results/3_9_1year2_noV

time ./z_rothC NB 100 testing 3 Zoutstep 1 C_q 1.0 corrector_type 0 Tm 3650 Om 1 LB 1 H0 -3 ls_eps 1e-10 ls_max_iter 200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth1.txt Ta_file Ta1.txt SoC_file SOC_inputs2.txt Kc_file Kc1.txt C00 2.78 C01 2.78 C02 2.78 C03 2.78
./move_to_folder.sh results/3_9
time ./z_rothC NB 100 testing 3 Zoutstep 1 C_q 0.99 corrector_type 0 Tm 3650 Om 1 LB 1 H0 -3 ls_eps 1e-10 ls_max_iter 200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth1.txt Ta_file Ta1.txt SoC_file SOC_inputs2.txt Kc_file Kc1.txt C00 2.78 C01 2.78 C02 2.78 C03 2.78
./move_to_folder.sh results/3_9_099_0
time ./z_rothC NB 100 testing 3 Zoutstep 1 C_q 0.99 corrector_type 1 Tm 3650 Om 1 LB 1 H0 -3 ls_eps 1e-10 ls_max_iter 200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth1.txt Ta_file Ta1.txt SoC_file SOC_inputs2.txt Kc_file Kc1.txt C00 2.78 C01 2.78 C02 2.78 C03 2.78
./move_to_folder.sh results/3_9_099_1

time ./z_rothC NB 100 testing 3 Zoutstep 1 C_q 1.0 corrector_type 0 Tm 5650 Om 1 LB 1 H0 -3 ls_eps 1e-10 ls_max_iter 200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt Ta_file Ta1.txt SoC_file SOC_inputs2.txt Kc_file Kc1.txt C00 2.78 C01 2.78 C02 2.78 C03 2.78
./move_to_folder.sh results/3_9_restore
