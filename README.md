# Supplementary Information for "Modelling Molecular Emitters in Organic LightEmitting Diodes with the Quantum MechanicalBespoke Force Field"

QUBE FF are generated with QUBEKit, exapmle configuration files are given in **QUBEKit_configs**. To run QUBEKit in command line, simply type:
   
    # PTZ s0 state
    QUBEKit -i PTZ.pdb -config ptz-s0.ini
    # or for s1 state
    QUBEKit -i PTZ.pdb -config ptz-s1.ini

For more details on QUBEKit and troubleshooting, see https://github.com/qubekit/QUBEKit.

The following folders contain:

## FF_and_pdb
 * ptz_s0.xml - parameterised QUBE FF for the s0 state of PTZ-DBTO2
 * ptz_s1.xml - parameterised QUBE FF for the s1 state of PTZ-DBTO2
 * CBP.xml - parameterised QUBE FF for the s0 state of CBP
 * PYD2.xml - parameterised QUBE FF for the s0 state of PYD2
 * ptz_cbp_box.pdb - a 5x5x5 nm^3 box of CBP molecules with a PTZ-DBTO2 at the centre
 * ptz_pyd2_box.pdb - a 5x5x5 nm^3 box of PYD molecules with a PTZ-DBTO2 at the centre
 
## sample_input_scripts
 * OpenMM_sampling_s0.py - a python script for running s0 state sampling using OpenMM
 * OpenMM_sampling_s1.py - a python script for running s1 state sampling using OpenMM
 * tda.inp - a Gaussian 09 input file for calculating the lowest 20 singlet and triplet excitations of PTZ-DBTO2 (TDDFT+TDA) with backgroud point charges

## trajectory
 * Raw trajectories of PTZ-DBTO2 S<sub>0</sub> and S<sub>1</sub> states samplings in different environments.
 * To access the data, see https://doi.org/10.5281/zenodo.4507905.
 * Each *.dcd file corresponds to a 100 ns MD simulation, which contains 500 frames. 
