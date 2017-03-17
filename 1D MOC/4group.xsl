ORNL 47-Group Library, Group 47
 4 5
 2.0E+07 1.0E5 1.0E3 1.0
!
!Comments can appear after the first 3 lines and between macro/micro blocks
!
!In the second line, the first number is number of groups and the other is
!number of cross section sets. 
!
!In the third line the energy group bounds are made up.
!
! Template for XS set
! XSMACRO <name> <scattering order>
! xs_absorption xs_nufission xs_fission chi ! repeat for each group
! Scattering matrix for scatters from column to row
! Repeat scattering matrix for each scattering order (starting with 0)
!

!Moderator
XSMACRO Moderator 0
  6.17E-04  0.00E+00  0.00E+00  0.00E+00
  2.28E-03  0.00E+00  0.00E+00  0.00E+00
  2.07E-02  0.00E+00  0.00E+00  0.00E+00
  3.72E-02  0.00E+00  0.00E+00  0.00E+00
  4.40E-01  0.00E+00  0.00E+00  0.00E+00
  1.31E-01  6.61E-01  7.14E-05  0.00E+00
  5.55E-05  4.99E-01  1.35E+00  1.32E-01
  1.05E-06  1.26E-02  5.99E-01  2.48E+00


!Guide tube
XSMACRO GuideTube 0
  5.87E-04  0.00E+00  0.00E+00  0.00E+00
  1.47E-03  0.00E+00  0.00E+00  0.00E+00
  1.26E-02  0.00E+00  0.00E+00  0.00E+00
  2.32E-02  0.00E+00  0.00E+00  0.00E+00
  3.66E-01  0.00E+00  0.00E+00  0.00E+00
  5.30E-02  3.55E-01  3.73E-05  0.00E+00
  2.22E-05  2.04E-01  6.24E-01  4.98E-02
  4.21E-07  5.14E-03  2.63E-01  1.10E+00


!Gap -- Actually just guide tube
XSMACRO Gap 0
  5.87E-04  0.00E+00  0.00E+00  0.00E+00
  1.47E-03  0.00E+00  0.00E+00  0.00E+00
  1.26E-02  0.00E+00  0.00E+00  0.00E+00
  2.32E-02  0.00E+00  0.00E+00  0.00E+00
  3.66E-01  0.00E+00  0.00E+00  0.00E+00
  5.30E-02  3.55E-01  3.73E-05  0.00E+00
  2.22E-05  2.04E-01  6.24E-01  4.98E-02
  4.21E-07  5.14E-03  2.63E-01  1.10E+00


! Material 4: Fuel
XSMACRO FUEL 0
  1.17E-02  2.21E-02  8.03E-03  1.00E+00
  1.23E-01  6.09E-02  2.50E-02  3.30E-04
  1.41E-01  2.45E-01  1.01E-01  0.00E+00
  2.83E-01  5.26E-01  2.16E-01  0.00E+00
  4.94E-01  0.00E+00  0.00E+00  0.00E+00
  1.64E-03  9.06E-01  1.25E-04  0.00E+00
  0.00E+00  5.57E-03  5.49E-01  8.55E-03
  0.00E+00  0.00E+00  1.68E-02  2.73E-01


!Control Rod
XSMACRO CRod 0
  1.01E-02  0.00E+00  0.00E+00  0.00E+00
  4.82E-01  0.00E+00  0.00E+00  0.00E+00
  1.63E+00  0.00E+00  0.00E+00  0.00E+00
  1.18E+00  0.00E+00  0.00E+00  0.00E+00
  6.86E-01  0.00E+00  0.00E+00  0.00E+00
  7.84E-04  1.37E+00  6.56E-05  0.00E+00
  0.00E+00  1.46E-03  4.15E-01  3.53E-03
  0.00E+00  0.00E+00  4.75E-03  6.59E-01