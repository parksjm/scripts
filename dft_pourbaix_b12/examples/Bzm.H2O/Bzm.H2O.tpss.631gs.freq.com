%nproc=24
%chk=Bzm.H2O.chk
%mem=7000mb

# tpsstpss/6-31G(d) freq integral=ultrafine scf=tight

Bzm.H2O

0 1 
N	1.331946	0.807411	-0.039766
C	0.198487	-0.009381	-0.024694
C	-0.979232	0.782362	0.011392
N	-0.514488	2.094025	0.015688
C	0.862093	2.038576	-0.016715
C	0.114514	-1.412135	-0.042212
C	-1.161442	-1.977939	-0.020219
C	-2.327791	-1.176890	0.017998
C	-2.262401	0.219775	0.034007
H	-1.078081	2.937044	0.033854
H	1.465953	2.939032	-0.022070
H	1.028556	-2.002782	-0.077750
H	-1.267242	-3.061702	-0.034703
H	-3.303851	-1.659802	0.033734
H	-3.163931	0.829891	0.061480
O	3.433490	-1.183843	-0.046463
H	3.565987	-1.251026	0.918357
H	2.897108	-0.356175	-0.129998

