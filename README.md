# supplementary information for "Modelling Molecular Emitters in Organic LightEmitting Diodes with the Quantum MechanicalBespoke Force Field"

The following folders contain:

## FF_and_pdb
 * ptz_s0.xml - parameterised QUBE FF for the s0 state of PTZ-DBTO2
 * ptz_s1.xml - parameterised QUBE FF for the s1 state of PTZ-DBTO2
 * CBP.xml - parameterised QUBE FF for the s0 state of CBP
 * PYD2.xml - parameterised QUBE FF for the s0 state of PYD2
 * ptz_cbp_box.pdb - a 5x5x5 nm^3 box of CBP molecules with a PTZ-DBTO2 at the centre
 * ptz_pyd2_box.pdb - a 5x5x5 nm^3 box of PYD molecules with a PTZ-DBTO2 at the centre
 
## sample_input_scripts
 * OpenMM_sampling_s0.py - a script for running s0 state sampling using OpenMM
 * OpenMM_sampling_s1.py - a script for running s1 state sampling using OpenMM
 * tda.inp - a Gaussian 09 input for calculating the lowest 20 singlet and triplet excitations of PTZ-DBTO2 (TDDFT+TDA + backgroud point charges)


