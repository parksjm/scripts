#!/usr/bin/python

# from thermochem_g09_scaled import *
from constants import *
import sys
import math
from sys import argv, exit

# This python3 script computes standard reduction potentials, pKas, and
# equilibrium constants for aquacobalamin.
#
# Toward Quantitatively Accurate Calculation of the Redox-Associated Acid–Base 
# and Ligand Binding Equilibria of Aquacobalamin
#
# RC Johnston, J Zhou, JC Smith, JM Parks
# Journal of Physical Chemistry B, 2016, 120, 7307–7318

# Written by Jing Zhou, Ryne C. Johnston and Jerry M. Parks
# Latest revision: Aug. 2019

# Born-Haber cycle
#                                       dG0_gas
# (gas) Ox1    +    [Ox2] + [...] + e-  -----> Red1     +     [Red2] + [...]
#
#        | dG_o->*    | dG_o->*                 | dG_o->*       | dG_o->*
#        v            v                         v               v
#        |            |                         |               |
#        | dG*_solv   | dG*_solv                | dG*_solv      | dG*_solv
#        v            v                         v               v
#                                      dG*_aq
# (Aq)  Ox1    +    [Ox2] + [...] + e-  -----> Red1     +      [Red2] + [...]
#
# Gaussian 09/16 Output File Naming Conventions
#
#
# All Gaussian output file names should include the name of the
# molecule/complex followed by the method and basis set.
#
# Geometry optimization output filenames should end in ".opt.log"
# Vibrational frequency output filenames should end in ".freq.log"
# Continuum solvation output filenames should end in ".<pcm>_<radii>.log",
#   where <pcm> = smd or cpcm and <radii> = bondi or uff
# Single-point energies with a larger basis set should include the name
# of the method and basis set and the filenames should end in ".sp.log"
#
# Examples:
# H2O.Co.3.Corrin.His.bp86.631gs.opt.log
# H2O.Co.3.Corrin.His.bp86.631gs.freq.log
# H2O.Co.3.Corrin.His.bp86.631gs.smd_bondi.log
# H2O.Co.3.Corrin.His.bp86.6311ppgdp.sp.log

# Python script input file format
# The input file should include one reaction per line.
# List all reactant species (oxidized complex, plus any other reactants)
# separated by commas but with no spaces between them, then one or more
# spaces, then list all product species (reduced complex plus any other
# products) separated by commas but with no spaces between them.
#
# Example:
# H2O.Co.3.Corrin.Bzm                              Co.2.Corrin.Bzm,H2O
# H2O.Co.3.Corrin.Bzm                              Co.2.Corrin.Bzm,H2O
# H2O.Co.3.Corrin.Bzm,H2O                          Co.2.Corrin.Bzm,H2O

# Functions
# Get single point energy


def get_dft_b1(b):
    temp_dft = []
    temp_alldft = []
    try:
        dftfile = open(mol + '/' + mol + "." + b + ".sp.log", 'r')
    except IOError:
        print("ERROR:" + mol + '/' + mol + "." + b + ".sp.log does not exist.")
        temp_dft = [0]
        temp_alldft = float(temp_dft[0])
        return temp_alldft
    else:
        for line in dftfile:
            if "SCF Done" in line:
                temp_dft.append(float(line.split()[4]))
        temp_alldft = temp_dft[-1]
        dftfile.close()
        return temp_alldft


# Get thermal corrections
def get_thermal(b, key, vibif):
  temp = []
  temp_all = []
  try:
    zpefile = open(mol + '/' + mol +  "." + b + ".freq.log", 'r')
#    freqfile = mol + '/' + mol +  "." + b + ".freq.log"
  except IOError:
    print("ERROR:" + mol + '/' + mol + "." + b + ".freq.log does not exist.")
    temp = [0]
    temp_all = float(temp[0])
    return temp_all
  else:
    for line in zpefile:
      if key == "Thermal correction to Gibbs Free Energy" and key in line:
        #if vibif == 'y':
        #  temp.append(thermochem(freqfile))
        #else:
          temp.append(float(line.split()[-1]))
      elif key == "SCF Done" and key in line:
        temp.append(float(line.split()[4]))
    temp_all = temp[0]
    zpefile.close()
    return temp_all


# Get energy in solution
def get_solv_E(s, b):
  temp_solv = []
  temp_gas = []
  temp_allgas = []
  temp_allsolv = []
  m = 0
  try:
    solvfile = open (mol + '/' + mol +  "." + b + "."+ s + ".log", 'r')
  except IOError:
    print("ERROR:" + mol + '/' + mol +  "." + b + "." + s + ".log does not exist.")
    temp_gas = [0]
    temp_solv = [0]
    temp_allgas = temp_gas[-1]
    temp_allsolv = temp_solv[-1]
    return temp_allgas, temp_allsolv
  else:
    for line in solvfile:
      if "SCF Done" in line:
        m += 1
    if m != 2:
      print('ERROR:The solvation file does not contain two SCF energies.')
    else:
      n = 0
      solvfile.seek(0)
      for line in solvfile:
        if "SCF Done" in line:
          if n == 0:
            temp_gas.append(float(line.split()[4]))
            n += 1
          else:
            temp_solv.append(float(line.split()[4]))
            
    temp_allgas = temp_gas[-1]
    temp_allsolv = temp_solv[-1]
    solvfile.close()

  return temp_allgas, temp_allsolv


# Get energy in gas phase (low level of theory)
def get_dft_b0(b, key):
  print('#### get_dft_b0 has been called!!! @@@@')
  temp = []
  temp_all = []
  try:
    b0file = open(mol + '/' + mol +  "." + b + ".opt.log", 'r')
  except IOError:
    print("ERROR:" + mol + '/' + mol + "." + b + ".opt.log does not exist.")
    temp = [0]
    temp_all = float(temp[0])
    return temp_all
  else:
    for line in b0file:
      if key == "SCF Done" and key in line:
        temp.append(float(line.split()[4]))
    temp_all = temp[-1]
    b0file.close()
    return temp_all

# Get dispersion corrections
def get_dftd3_b0(b, key, damping):
  temp = []
  temp_all = []
  try:
    b0file = open(mol + '/' + mol + "." + b + ".opt." + damping, 'r')
  except IOError:
    print("ERROR:" + mol + '/' + mol + "." + b + ".opt." + damping + " does not exist.")
    temp = [0]
    temp_all = float(temp[0])
    return temp_all
  else:
    for line in b0file:
      if key == "Edisp" and key in line:
        temp.append(float(line.split()[3]))
    temp_all = temp[-1]
    b0file.close()
    return temp_all

# MAIN

if len(argv) < 3:
   print("Usage: <script.py> <input_file> <E0/pKa/Keq> Optional: -d3bj/-d3zero -l")
   exit()

#Handle user inputs
datatyp = str.lower(argv[2])

for n, arg in enumerate(argv):
  if argv[n] == '-d3bj':
    d3bj = 'y'
  elif argv[n] == '-d3zero':
    d3zero = 'y'
  #elif argv[n] == '-vib':
  #  vib = 'y'
  elif argv[n] == '-l':
    long = 'y'
  else:
    d3bj = 'n'
    d3zero = 'n'
    vib = 'n'
    long = 'n'

# Read oxidized and reduced species from the input file.
temp_list=[]
ox_list = []
red_list = []
exp_list = []

error_list = []
error_list_abs = []
error_list_sq = []

with open(sys.argv[1], 'r') as my_file:
   for line in my_file:
      if list(line)[0] != "#":
        temp_list=line.split()
        ox_list.append(list(temp_list[0].split(',')))
        red_list.append(list(temp_list[1].split(',')))
        # Check to see if exptl E0/pKa values are included in the input file.
        if len(temp_list) == 3:
          exp_list.append(temp_list[2])
        else:
          exp_list.append(9999)

for i in range(len(exp_list)):
   exp_list[i] = float(exp_list[i])

# Select the level of theory, solvent model, etc. Currently hard-coded.

sp = 'tpss.6311ppgdp'
freq = 'tpss.631gs'
solvtheory = 'tpss.6311pgss'
solv = 'smd_bondi'

ssss = [solvtheory, solv]
print('THEORY:  ' + '.'.join(ssss) + '/' + sp + '//' + freq)

# Print header
if datatyp == 'e0':
  typ = 'E (SCE)'
elif datatyp == 'pka':
  typ = 'pKa'
elif datatyp == 'keq' or datatyp == 'logkeq':
  typ = 'logKeq'
  print('Displaying equilibrium constants as logKeq rather than as Keq for easier analysis.')
else:
  print('Unrecognized calculation type. Options are k')

if long == 'n':
  print("%-50s %5s %-50s %9s %9s %9s %9s %9s %9s %9s %9s %9s" %
   ('Reactants', '     ', 'Products', 'dE_g', 'dG_corr', 'dG_gas', 
    'dG_s(red)', 'dG_s(ox)', 'dG0_aq', 'Expt', typ, 'Error'))

for i in range(len(ox_list)):
  b0_gas_ox = []
  b1_gas_ox = []
  b0_pcm_ox = []
  b0_gcorr_ox = []
  b0_disp_ox = []
  H2O_ipso_solv_ox = []

  b0_gas_red = []
  b1_gas_red = []
  b0_pcm_red = []
  b0_gcorr_red = []
  b0_disp_red = []
  H2O_ipso_solv_red = []

  dG_gas = []
  dG_solv_ox = []
  dG_solv_red = []
  ddG_solv = []

  E0_SCE = []
  E0_SHE = []

  # STEP 1: Get energies for oxidized species###     
  for mol in ox_list[i]:
    # 1.1 Get single point energies (b1)
    if mol == "H+":
      b1_gas_ox.append(0.0)
    else:
      b1_gas_ox.append(get_dft_b1(sp))

    # 1.2 Get thermal corrections
    if mol == "H+":
      b0_gcorr_ox.append(G0_proton)
    else:
      key = "Thermal correction to Gibbs Free Energy"
      b0_gcorr_ox.append(get_thermal(freq, key, vib))


# 1.3 Get SCF energies in PCM (b0)
    if mol == "H+":
      b0_gas_ox.append(0.0)
      b0_pcm_ox.append(dGssolv_proton)
    else:
      b0_gas_ox.append(get_solv_E(solv,solvtheory)[0])
      b0_pcm_ox.append(get_solv_E(solv,solvtheory)[1])

    #1.4 Get SCF energies in gas phase (b0) 
#    if mol == "H+":
#       b0_gas_ox.append(0.0)
#    else:
#      key = "SCF Done"
#      b0_gas_ox.append(get_thermal(freq, key, vib))

    #2.5 Get G_liquid to G_star corrections 
    if mol == "H2O" or mol == 'H2O_c2v' or mol == 'H2O_fake':
      H2O_ipso_solv_ox.append(liquid2star)
    else:
      H2O_ipso_solv_ox.append(0.0)

    #2.6 Get dispersion corrections 
    if d3bj == 'y' or d3zero == 'y':
      if mol == "H+":
        b0_disp_red.append(0.0)
      else:
        key = "Edisp"
        if d3bj == 'y':
          damping = 'd3bj'
        elif d3zero == 'y':
          damping = 'd3zero'
        b0_disp_ox.append(get_dftd3_b0(freq, key, damping))

#STEP 2: Get energies for reduced species###
  for mol in red_list[i]:

    # 2.1 Get single-point energy
    if mol == "H+":
      b1_gas_red.append(0.0)
    else:
      b1_gas_red.append(get_dft_b1(sp))

    # 2.2 Get thermal correction
    if mol == "H+":
      b0_gcorr_red.append(G0_proton)
    else:
      key = "Thermal correction to Gibbs Free Energy"
      b0_gcorr_red.append(get_thermal(freq, key, vib))

# 2.3 Get SCF energies in PCM (b0)

    if mol == "H+":
      b0_gas_red.append(0.0)
      b0_pcm_red.append(dGssolv_proton)
    else:
      b0_gas_red.append(get_solv_E(solv,solvtheory)[0])
      b0_pcm_red.append(get_solv_E(solv,solvtheory)[1])

# 2.4 Get SCF energies in the gas phase (b0) 
#    if mol == "H+":
#      b0_gas_red.append(0.0)
#    else:
#      key = "SCF Done"
#      b0_gas_red.append(get_thermal(freq, key, vib))
#    print(b0_gas_red[:1])

# 2.5 Get G_liquid to G_star corrections for water "cluster" of one water 
    if mol == "H2O" or mol == 'H2O_c2v' or mol == 'H2O_fake':
      H2O_ipso_solv_red.append(liquid2star)
    elif mol == "H2O_n6":
      H2O_ipso_solv_red.append(1.0*RT*math.log(55.34/6.0))
    else:
      H2O_ipso_solv_red.append(0.0)

# 2.6 Get dispersion corrections 
    if d3bj == 'y' or d3zero == 'y':
      if mol == "H+":
        b0_disp_red.append(0.0)
      else:
        key = "Edisp"
        if d3bj == 'y':
          damping = 'd3bj'
        elif d3zero == 'y':
          damping = 'd3zero'
        b0_disp_red.append(get_dftd3_b0(freq, key, damping))


  ####STEP 3: Compute E0###
# 3.1 Compute dG_gas (with thermal corrections)
  dE_gas = h2kcal*(sum(b1_gas_red) - sum(b1_gas_ox))
  if datatyp == 'e0':
    dG_thermal = h2kcal*(sum(b0_gcorr_red) - sum(b0_gcorr_ox) - Gs_electron)
  else:
    dG_thermal = h2kcal*(sum(b0_gcorr_red) - sum(b0_gcorr_ox))    
  if d3bj == 'y' or d3zero == 'y':
    dE_disp  = h2kcal*(sum(b0_disp_red) - sum(b0_disp_ox))
    dG_gas   = dE_gas + dG_thermal + dE_disp
  else:
    dG_gas   = dE_gas + dG_thermal          # 1 atm/mol standard state

  dG_gas2star_ox = []
  dG_gas2star_red = []

  if long == 'y':
    print('')
    print("%-20s %27s %-21s %20s" % ('******************', 'REACTION ', str(i+1), '******************'))
    print("%2s %-30s %13s %9s %9s %9s %13s" % ('', 'Species', 'Ggas', 'dG*solv', 'dG0>*', 'dGliq>*', 'Gaq'))

# 3.2.1 Compute dG_solv_ox 
  for j in range(len(b0_gas_ox)):
    Ggas_o = h2kcal*(b1_gas_ox[j] + b0_gcorr_ox[j])
    dGs_solv_o = h2kcal*(b0_pcm_ox[j] - b0_gas_ox[j])
    dG_solv_ox.append(h2kcal*(b0_pcm_ox[j] - b0_gas_ox[j]))
    dG_gas2star_ox.append(gas2star)
    Gaq_o = h2kcal*(b1_gas_ox[j] + b0_gcorr_ox[j]) + gas2star + H2O_ipso_solv_ox[j] + h2kcal*(b0_pcm_ox[j] - b0_gas_ox[j])

    if long == 'y':
      print("%-2s %-30s %13.2f %9.2f %9.2f %9.2f %13.2f" % ('r', ox_list[i][j], Ggas_o, dGs_solv_o, dG_gas2star_ox[j], H2O_ipso_solv_ox[j], Gaq_o))

  if datatyp == 'e0' and long == 'y':
    print("%-2s %-30s %13.2f %9.2f %9.2f %9.2f %13.2f" % ('r', 'e-',  h2kcal*(G0_electron), 0.00, gas2star, H2O_ipso_solv_ox[j], h2kcal*(Gs_electron)))

# 3.2.2 Compute dG_solv_red 
  for k in range(len(b0_gas_red)):
    Ggas_r = h2kcal*(b1_gas_red[k] + b0_gcorr_red[k])
    dGs_solv_r = h2kcal*(b0_pcm_red[k] - b0_gas_red[k])
    dG_solv_red.append(h2kcal*(b0_pcm_red[k] - b0_gas_red[k]))
    dG_gas2star_red.append(gas2star)
    Gaq_r = h2kcal*(b1_gas_red[k] + b0_gcorr_red[k]) + gas2star + H2O_ipso_solv_red[k] + h2kcal*(b0_pcm_red[k] - b0_gas_red[k])

    if long == 'y':
      print("%-2s %-30s %13.2f %9.2f %9.2f %9.2f %13.2f" % ('p', red_list[i][k], Ggas_r, dGs_solv_r, dG_gas2star_red[k], H2O_ipso_solv_red[k], Gaq_r))

# 3.3 Compute aqueous free energy of reduction
  gas2star_corr = sum(dG_gas2star_red) - sum(dG_gas2star_ox)
  liquid2star_corr = sum(H2O_ipso_solv_red) - sum(H2O_ipso_solv_ox)
  dG0_aq = dG_gas + gas2star_corr + sum(dG_solv_red) - sum(dG_solv_ox) + liquid2star_corr

# 3.4 Compute E0
  E0_SCE = -1.0*dG0_aq / faraday - SCE 
  pKa = dG0_aq / (RT * math.log(10.0)) 
  Keq = math.exp(-1.0 * dG0_aq / RT)

  if exp_list[i] != 9999.0:
    if datatyp == 'e0':
      error = E0_SCE - exp_list[i]
    elif datatyp == 'pka':
      error = pKa - exp_list[i]
    else:
      error = math.log10(Keq) - math.log10(exp_list[i])

  error_list.append(float(error))
  error_list_abs.append(abs(float(error)))
  error_list_sq.append(float(error*error))

  if datatyp == 'e0':
    finaldatum = E0_SCE
  elif datatyp == 'pka':
    finaldatum = pKa
  else:
    finaldatum = math.log10(Keq)

# 3.5 Print output
  reactants = ' + '.join(ox_list[i])
  products  = ' + '.join(red_list[i])
  arrow     = '---> '

  if long == 'y':
    print("-------------------------------------------------------------------------------------------")
    print("%-2s %-30s %13.2f %9.2f %9.2f %9.2f %13.2f" % ('', 'Reaction energies', dG_gas, sum(dG_solv_red) - sum(dG_solv_ox), gas2star_corr, liquid2star_corr, dG0_aq))
    print("%-2s %-74s %13.2f" % ('', typ, finaldatum))
    if exp_list[i] != 9999.0:
      print("%-2s %-74s %13.2f" % ('', 'Error', error))
  else:
    if exp_list[i] != 9999.0:
      if datatyp == 'logkeq':
        print("%-52s %5s %-50s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f" % (reactants, '---> ', products, dE_gas, dG_thermal, dG_gas, sum(dG_solv_red), sum(dG_solv_ox), dG0_aq, math.log10(exp_list[i]), finaldatum, error ))
      else:
        print("%-52s %5s %-50s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f" % (reactants, '---> ', products, dE_gas, dG_thermal, dG_gas, sum(dG_solv_red), sum(dG_solv_ox), dG0_aq, exp_list[i], finaldatum, error ))
    else:
      print("%-52s %5s %-50s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f" % (reactants, '---> ', products, dE_gas, dG_thermal, dG_gas, sum(dG_solv_red), sum(dG_solv_ox), dG0_aq, finaldatum ))

# Error analysis
# print('error_list', error_list)
# print('error_list_abs', error_list_abs)
# print('error_list_sq', error_list_sq)

print("%-189s %9.2f" % ('MSE', sum(error_list)/float(len(error_list))))
print("%-189s %9.2f" % ('MUE', sum(error_list_abs)/float(len(error_list))))
print("%-189s %9.2f" % ('RMSE', math.sqrt(sum(error_list_sq)/float(len(error_list)))))
min_err = min(err for err in error_list)
max_err = max(err for err in error_list)
if (abs(min_err) > abs(max_err)):
  print("%-189s %9.2f" % ('MAXE', min(err for err in error_list)))
else:
  print("%-189s %9.2f" % ('MAXE', max(err for err in error_list)))
