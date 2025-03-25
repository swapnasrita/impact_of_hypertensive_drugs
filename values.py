# The values here are updated from the male human model. 

N=200
# physical constants
RT    = 2.57
RTosm = 1.93e4
F     = 96.5
visc  = 6.4e-6
# number of variables for H-K-ATPase
Natp = 14

# Vwbar = Molar volume of water  (cm3/mmole)
# Formula weight of H20 = 2*1.00797 + 15.9994 = 18.0153
# There is 1 mole H20/18.0153 gm of H20.
# Density of water at 37 C = 0.99335 gm/cm3 (CRC Handbook)
# Vwbar = (gm/millimole)/(gm/cm3)=(cm3/millimole)

Vwbar = 18.0153/1000/0.99335

# normalization constants
Cref = 1e-3  # 1 mM = 0.01 mmol/cm3
Vref = 1e-4  # in cm3/cm2 epith
Pfref = 1.0   # in cm/s
href = 1e-5  # in cm/s
EPref= 1e-3  # in volts

# acid-base pair balance constants
pKHCO3 = 3.57
pKHPO4 = 6.80
pKNH3 = 9.15
pKbuf = 7.5
pKHCO2 = 3.76

# transporter parameters for NKCC2 F isoform
poppnkccF = 3.9280e4
popnkccF = 3.578e5
pnkccpF = 1e4
pnkccppF = 1.098e3
pnmccpF = 2.0e3
pnmccppF = 219.6
bn2F = 58.93
bc2F = 13.12
bk2F = 9.149
bm2F = 9.149

# transporter parameters for NKCC2 A isoform
poppnkccA = 7.535e4
popnkccA = 2.594e5
pnkccpA = 1e4
pnkccppA = 2.904e3
pnmccpA = 2.0e3
pnmccppA = 580.8
bn2A = 118.8
bc2A = 0.08834
bk2A = 1.8710e4
bm2A = 1.8710e4

# transporter parameters for NKCC2 B isoform
poppnkccB = 2.517e5
popnkccB = 2.596e5
pnkccpB = 1e4
pnkccppB = 9.695e3
pnmccpB = 2.0e3
pnmccppB = 1.939e3
bn2B = 275.0
bc2B = 0.08157
bk2B = 5.577e3
bm2B = 5.577e3

# transporter parameters for KCC (KCC4 isoform)
poppkcc = 3.928e4
popkcc = 3.577e5
pkccp = 1.0e4
pkccpp = 1.098e3
pmccp = 2.0e3
pmccpp = 219.6
bckcc = 21.08
bkkcc = 1.45
bmkcc = 1.45

# transporter parameters for NCC
popncc = 4.295e6
poppncc = 0.1e6
pnpncc = 7692.0
pnppncc = 179.0
dKnncc = 0.293
dKcncc = 112.7
dKncncc = 0.565
dKcnncc = 1.47e-3

# transporter parameters for Pendrin
Pclepd = 10000.0
Pcliclepd = 1.239
Pbieclepd = 10.76
Poheclepd = 0.262
dKclpd = 3.01
dKbipd = 5.94
dKohpd = 1.38e-6

# transporter parameters for AE1
Pbp = 1247.0
Pbpp = 135.0
Pcp = 562.0
Pcpp = 61.0
dKbp = 198.0
dKbpp = 198.0
dKcp = 50.0
dKcpp = 50.0

# parameters
pHplasma = 7.4
phpap = 7.3
pKHCO3 = 3.57
pKHPO4 = 6.8
pKNH3 = 9.15
pKHCO2 = 3.76