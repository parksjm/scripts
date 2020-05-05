%nproc=24
%chk=H2O_c2v.chk
%mem=7000mb
# tpsstpss/6-31G(d) opt integral=ultrafine scf=tight

H2O_c2v

0 1 
O	0.000000	0.000000	0.121389
H	0.000000	0.766366	-0.485557
H	0.000000	-0.766366	-0.485557

