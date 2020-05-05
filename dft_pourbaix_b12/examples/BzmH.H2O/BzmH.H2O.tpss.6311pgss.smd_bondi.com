%nproc=24
%chk=BzmH.H2O.chk
%mem=7000mb

# tpsstpss/6-311+G(d,p) integral=ultrafine scf=tight

BzmH.H2O

1 1 
N	1.305136	0.488766	-0.000015
C	0.069981	-0.166250	0.000012
C	-0.921348	0.840727	-0.000006
N	-0.223731	2.054400	0.000015
C	1.100597	1.809647	0.000002
C	-0.258258	-1.528564	0.000027
C	-1.619892	-1.831862	0.000014
C	-2.612811	-0.821346	-0.000011
C	-2.288569	0.536292	-0.000024
H	2.235778	0.003837	-0.000018
H	-0.636135	2.984122	0.000030
H	1.872462	2.567058	0.000006
H	0.509314	-2.297885	0.000045
H	-1.930190	-2.874123	0.000021
H	-3.660615	-1.112372	-0.000015
H	-3.054118	1.308069	-0.000034
O	3.512894	-1.111440	-0.000001
H	4.086265	-1.260463	0.776314
H	4.086052	-1.260753	-0.776419

--link1--

%nproc=24
%mem=7000mb
%chk=BzmH.H2O.chk
%nosave

# tpsstpss chkbasis integral=ultrafine scf=tight geom=check scrf=(smd,solvent=water, read)

BzmH.H2O

1 1

radii=Bondi

