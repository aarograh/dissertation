ORNL 47-Group Library, Group 47
 1 5
 2.0E+07
!
!Comments can appear after the first 3 lines and between macro/micro blocks
!
!In the second line, the first number is number of groups and the other is
!number of cross section sets. 
!
!In the third line the energy group bounds are made up.
!
! Template for XS set.  Keyword is "xsmacro" in ALL CAPS
! <keyword> <name> <scattering order>
! xs_absorption xs_nufission xs_fission chi ! repeat for each group
! Scattering matrix for scatters from column to row
! Repeat scattering matrix for each scattering order (starting with 0)
!

! Material 1: Moderator
XSMACRO MOD 0
4.394922836 0.0 0.0 0.0
0.0

! Material 2: Clad
XSMACRO CLAD 0
0.299908964 0.0 0.0 0.0
0.0

! Material 3: Gap
XSMACRO GAP 0
3.44E-05 0.0 0.0 0.0
0.0

! Material 4: Fuel
XSMACRO FUEL 0
1.489061075 0.0 0.0 0.0
0.0

! Material 5: Control rod
XSMACRO CONTROL 0
19.13440383 0.0 0.0 0.0
0.0