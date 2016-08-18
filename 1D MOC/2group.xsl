ORNL 47-Group Library, Group 47
 2 5
 2.0E+07 1.0
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

! Material 1: Moderator
XSMACRO MOD 0
1.394922836 0.0 0.0 0.0
3.0 0.0 0.0 0.0
0.5 0.0
0.25 0.5

! Material 2: Clad
XSMACRO CLAD 0
0.099908964 0.0 0.0 0.0
0.2 0.0 0.0 0.0
0.5 0.0
0.05 0.5

! Material 3: Gap
XSMACRO GAP 0
3.44E-05 0.0 0.0 0.0
3.44E-05 0.0 0.0 0.0
1.0E-05 0.0
0.0 1.0E-05

! Material 4: Fuel
XSMACRO FUEL 0
0.389061075 0.3 0.125 0.95
0.6 0.9 0.375 0.05
0.8 0.0
0.5 1.2

! Material 5: Control rod
XSMACRO CONTROL 0
1.13440383 0.0 0.0 0.0
18.0 0.0 0.0 0.0
1.0 0.0
0.5 1.0