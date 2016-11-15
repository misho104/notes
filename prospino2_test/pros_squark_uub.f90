program main
  use xx_kinds
  use xx_prospino_subroutine
  implicit none

  integer                              :: inlo,isq_ng_in,icoll_in,i_error_in,ipart1_in,ipart2_in,isquark1_in,isquark2_in
  real(kind=double)                    :: energy_in
  logical                              :: lfinal
  character(len=2)                     :: ss, sb

  inlo = 0          ! specify LO only[0] or complete NLO (slower)[1]
  isq_ng_in = 1     ! specify degenerate [0] or free [1] squark masses
  icoll_in = 1      ! collider : tevatron[0], lhc[1]
  energy_in = 13000 ! collider energy in GeV
  i_error_in = 0    ! with central scale [0] or scale variation [1]
  ipart1_in = 1
  ipart2_in = 1

  ss = 'ss'
  sb = 'sb'
!---------------------------------------------------------------------------!
!  for LO with light-squark flavor in the final state                       !
!  isquark1_in     =  -5,-4,-3,-2,-1,+1,+2,+3,+4,+5                         !
!                    (bL cL sL dL uL uR dR sR cR bR) in CteQ ordering       !
!  isquark1_in     = 0 sum over light-flavor squarks throughout             !
!                      (the squark mass in the data files is then averaged) !
!  flavors in initial state: only light-flavor partons, no bottoms          !
!                            bottom partons only for Higgs channels         !
!  flavors in final state: light-flavor quarks summed over five flavors     !
!---------------------------------------------------------------------------!
  
  call PROSPINO_OPEN_CLOSE(0)
  call PROSPINO(inlo, isq_ng_in, icoll_in, energy_in, i_error_in, sb, ipart1_in, ipart2_in, 4, 4)
  call PROSPINO_OPEN_CLOSE(1)
end program main
