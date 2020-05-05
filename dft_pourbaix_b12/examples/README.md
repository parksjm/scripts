
# Reduction potential example
# Computes the proton-coupled reduction potential for diagonal M in Figure 2.

./master.py.v14 redox_saturated.txt E0 -d3bj

# pKa example
# Computes the pKa for the protonation of the DMB tail in base-off cobalamin. 

./master.py.v14 pka_saturated.txt pka -d3bj

# K_on/off example
# Computes the log K value for association/dissociation of the DMB tail in cob(II)alamin

./master.py.v14 keq_saturated.txt keq -d3bj
