%nproc=24
%chk=H2O_c2v.chk
%mem=7000mb

# tpsstpss/6-311+G(d,p) integral=ultrafine scf=tight

H2O_c2v

0 1 
O	0.000000	0.000000	0.121570
H	0.000000	0.763722	-0.486280
H	0.000000	-0.763722	-0.486280

--link1--

%nproc=24
%mem=7000mb
%chk=H2O_c2v.chk
%nosave

# tpsstpss chkbasis integral=ultrafine scf=tight geom=check scrf=(smd,solvent=water)

H2O_c2v

0 1

