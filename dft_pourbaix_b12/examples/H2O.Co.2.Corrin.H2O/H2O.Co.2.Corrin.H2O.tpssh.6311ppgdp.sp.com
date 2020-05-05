%nproc=24
%chk=H2O.Co.2.Corrin.H2O.chk
%mem=7000mb

# tpssh/6-311++G(d,p) integral=ultrafine scf=tight

H2O.Co.2.Corrin.H2O

1 2 
O	-0.182428	0.206267	-2.306619
H	0.065529	-0.609522	-2.785100
H	0.392992	0.903937	-2.678368
Co	-0.030963	0.000000	-0.000001
N	-1.326193	-1.415207	-0.035824
C	-1.039803	-2.773733	0.034673
C	-2.311344	-3.593028	0.168131
C	-3.426992	-2.556575	-0.064304
C	-2.679085	-1.238636	-0.027661
C	0.224801	-3.323970	0.015984
C	1.405957	-2.546253	-0.096571
N	1.366297	-1.223720	-0.145445
C	2.730415	-0.666532	-0.372121
C	3.671073	-1.802138	0.072943
C	2.827634	-3.070936	-0.211231
C	2.730417	0.666519	0.372125
N	1.366303	1.223714	0.145445
C	1.405968	2.546247	0.096574
C	2.827647	3.070923	0.211240
C	3.671081	1.802121	-0.072935
C	-3.314357	0.000007	0.000004
C	-2.679079	1.238647	0.027666
N	-1.326186	1.415212	0.035822
C	-1.039790	2.773737	-0.034680
C	-2.311328	3.593034	-0.168149
C	-3.426978	2.556591	0.064317
C	0.224816	3.323969	-0.015986
H	2.843799	-0.465782	-1.452107
H	2.843797	0.465768	1.452112
H	-4.402474	0.000010	0.000008
H	0.317579	4.405354	-0.083599
H	0.317559	-4.405356	0.083595
H	3.870695	1.718408	-1.150872
H	4.630925	1.794981	0.454413
H	-4.221190	-2.593578	0.691038
H	3.027942	3.895313	-0.484156
H	-2.368654	4.027642	-1.176420
H	3.005459	3.455700	1.228105
H	-4.221194	2.593596	-0.691005
H	3.870685	-1.718423	1.150880
H	4.630918	-1.795003	-0.454403
H	3.005447	-3.455718	-1.228094
H	3.027923	-3.895325	0.484167
H	-2.332577	4.426171	0.544365
H	-2.332600	-4.426145	-0.544407
H	-3.907606	2.685628	1.045314
H	-2.368666	-4.027664	1.176388
H	-3.907647	-2.685602	-1.045288
O	-0.182445	-0.206261	2.306609
H	0.392972	-0.903928	2.678370
H	0.065503	0.609532	2.785088

