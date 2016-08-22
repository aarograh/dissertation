ORNL 47-Group Library, Group 47
 4 1
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

! Material 4: Fuel
XSMACRO FUEL 0
8.0248E-03 2.005998E-02 7.21206E-03 1.0
3.7174E-03 2.027303E-03 8.19301E-04 0.0
1.1126E-01 2.020901E-01 8.30348E-02 0.0
2.8278E-01 5.257105E-01 2.16004E-01 0.0
1.27537E-01 0.00000E+00 0.00000E+00 0.00000E+00
4.23780E-02 3.24456E-01 0.00000E+00 0.00000E+00
0.00000E+00 1.00000E-02 2.65802E-01 0.0
0.00000E+00 1.00000E-03 1.68090E-02 2.73080E-01