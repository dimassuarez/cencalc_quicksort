
!*****************************************************************************************
MODULE parameters                                                                !
!*****************************************************************************************
        implicit none
        integer, parameter :: DP = 8 
        integer, parameter :: K4B=selected_int_kind(9)  !  this is kind=4
        integer, parameter :: K4C=selected_int_kind(10) !  this is kind=4 
        real(DP),parameter :: R=8.3144621d0/4.184d0     !Gas Constant in cal/molK
        real(DP),parameter :: Pi=3.14159265359d0
END module parameters
!-----------------------------------------------------------------------------------------

