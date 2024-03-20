# Second harmonic generation with the FDTD technique
This code solves Maxwell's equations using the finite-difference time-domain (FDTD) method. It incorporates a classical Lorentz model modified by including a quadratic term in the polarization to generate the second harmonic.

The code was developed as part of a master's thesis project, the full text of which can be found at the following link: 
https://gredos.usal.es/bitstream/handle/10366/151749/tfm_isabel_rodriguez_perez.pdf?sequence=1&isAllowed=y

## Dispersion and SHG 
This code simulates both linear and nonlinear pulse propagation in vacuum and dielectrics. Dispersion in the dielectric is introduced through the differential equation provided by the Lorentz model. Depending on the region of interest, the pulse experiences either positive or negative chirping as it propagates through the dispersive medium due to Group Velocity Dispersion (GVD). The difference between phase velocity and group velocity is also highlighted. Nonlinear propagation in the medium, under suitable conditions, leads to the generation of the second harmonic and higher-order harmonics.

## Birrefringence
After achieving the second harmonic, a birefringent material is implemented, which may be in phase or phase mismatch, generating a second harmonic in the direction perpendicular to the incident. In this section, we include code to simulate a nonlinear birefringent dispersive medium with an incident pulse with the electric field oriented along the Z direction propagating in the X direction. Due to phase adjustment through birefringence, an electric field oriented along the Y axis and a magnetic field oriented in the Z direction are generated, both perpendicular to the propagation direction in X.

## Bibliography
Rose M Joseph and Allen Taflove. FDTD maxwell’s equations models for nonlinear electrodynamics and optics. IEEE Transactions on Antennas and Propagation, 45(3):364–374, 1997.
John B Schneider. Understanding the finite-difference time-domain method. School of electrical engineering and computer science Washington State University, 28, 2010.
