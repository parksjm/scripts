%nproc=24
%chk=H2O_c2v.chk
%mem=7000mb

# tpsstpss/6-31G(d) freq integral=ultrafine scf=tight

H2O_c2v

0 1 
O	0.000000	0.000000	0.121570
H	0.000000	0.763722	-0.486280
H	0.000000	-0.763722	-0.486280

