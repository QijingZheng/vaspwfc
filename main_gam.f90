program vaspwfc
  use prec
  use info
  use lattice
  use gvector_gam
  use wavecar
  use constants
  use vaspwfc_gam

  implicit none

  type(system) :: MySys
  integer :: i, j
  integer :: Nb, Nk
  ! type(psi), allocatable, dimension(:) :: BRAS, KETS
  type(psi) :: kets

  character(len=128) :: buf

  MySys%WAVECAR = "WAVECAR"

  do i = 1, 3
    MySys%NGPTAR(i) = -1
  end do

  write(*,*) "Reset grid size (NGX, NGY, NGZ)? (yes/no)"
  read(*,*) buf
  buf = ADJUSTL(buf)
  if ((buf(1:1) == 'y') .or. (buf(1:1) == 'Y')) then
    write(*,*) "Enter NGX, NGY, NGZ:"
    read(*,*) MySys%NGPTAR(1), MySys%NGPTAR(2), MySys%NGPTAR(3)
    ! write(*,*) MySys%NGPTAR(1), MySys%NGPTAR(2), MySys%NGPTAR(3)
  end if

  ! initialize the system
  call INITIALIZE("WAVECAR", MySys)
  call WRITEINFO(MySys)

  write(*,*) 'Enter ISPIN, IKPTS, IBANDS'
  read(*,*), kets%ispin, kets%ikpts, kets%iband

  call wavefft(kets, MySys)

  call freemem(MySys)

end program vaspwfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize(WAVEFILE, MySys)
  use info
  use lattice
  use gvector_gam
  use wavecar
  implicit none

  type(system), intent(inout) :: MySys
  character(len=*) :: WAVEFILE

  ! read ISPIN, NKPTS, NBANDS, ENCUT, LATTICE
  call sysinfo(MySys, WAVEFILE)
  ! get reciprocal space basis from real one
  call lattic(MySys%MyLatt)
  ! get NGPTAR, and G-vector index
  call GCUTS(MySys)
  ! count number of plane waves for each k-point and compare the value with the
  ! one read from WAVECAR. Also, the index of each G-vector are produced.
  call CNT_PLW_NUM(MySys)

end subroutine

subroutine writeinfo(MySys)
  use info
  type(system), intent(in) :: MySys

  integer :: i, j

  write(*,'(A)') '************************************************************'
  write(*, '(A, I4)') 'ISPIN = ', MySys%ISPIN
  write(*, '(A, I4)') 'NKPTS = ', MySys%NKPTS
  write(*, '(A, I4)') 'NBAND = ', MySys%NBANDS
  write(*, '(A, F8.2, A)') 'ENCUT = ', MySys%ENCUT, ' eV'
  write(*, '(A, 3I5)') 'NGPTAR: ', (MySys%NGPTAR(i), i=1,3)
  write(*, '(A)') 'Real Space Basis:'
  write(*, '(3F10.4)') ((MySys%MyLatt%A(i,j), i=1,3), j=1,3)
  write(*, '(A)') 'Reciprocal Space Basis:'
  write(*, '(3F10.4)') ((MySys%MyLatt%B(i,j), i=1,3), j=1,3)
  write(*,'(A)') '************************************************************'
end subroutine
