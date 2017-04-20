module vaspwfc
  use prec
  use info
  use constants
  use lattice
  use wavecar

  use, intrinsic :: iso_c_binding 
  implicit none

  include 'fftw3.f03'

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate FT of single KS wavefunction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wavefft(ket, MySys)
    implicit none

    type(psi), intent(in) :: ket
    type(system), intent(in) :: MySys

    integer :: i, j, k, nplw
    complex(kind=qs), allocatable, dimension(:) :: cket

    type(C_PTR) :: plan
    complex(C_DOUBLE_COMPLEX), dimension(MySys%NGPTAR(1), MySys%NGPTAR(2), MySys%NGPTAR(3)/2+1) :: wfc_k
    real(C_DOUBLE), dimension(MySys%NGPTAR(1), MySys%NGPTAR(2), MySys%NGPTAR(3)) :: wfc_r
    
    call openwav(trim(MySys%WAVECAR))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(cket(MySys%NPLWS(ket%ikpts)))
    ! load the G coefficients
    call LOADWAVE(cket, ket, MySys)

    ! initialize the input array to zero
    wfc_k = (0., 0.)
    ! from 1d coefficients to 3d grid
    do nplw = 1, MySys%NPLWS(ket%ikpts)
      i = MySys%grid_index(1,nplw,ket%ikpts)
      j = MySys%grid_index(2,nplw,ket%ikpts)
      k = MySys%grid_index(3,nplw,ket%ikpts)
      wfc_k(i,j,k) = cket(nplw)
    end do
    ! add those 
    do j = 1, MySys%NGPTAR(2)
      do i = 1, MySys%NGPTAR(1)
        if (MySys%LPCTY(j)>0) cycle
        if (MySys%LPCTY(j)==0 .AND. MySys%LPCTX(i)>0) cycle

        wfc_k(i,j,k) = conjg(wfc_k(MySys%NGPTAR(1) + 1 - i, MySys%NGPTAR(2) + 1 - j, 1))
      end do
    end do

    plan = fftw_plan_dft_c2r_3d(MySys%NGPTAR(3), MySys%NGPTAR(2), MySys%NGPTAR(1), wfc_k, wfc_r, FFTW_BACKWARD)

    call fftw_execute_dft_c2r(plan, wfc_k, wfc_r)
    call fftw_destroy_plan(plan)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call outwfc(wfc_r, ket, MySys)

    call closewav()

    deallocate(cket)
  end subroutine

  subroutine outwfc(wfc_r, ket, MySys)
    implicit none

    real(C_DOUBLE), intent(in), dimension(:,:,:) :: wfc_r
    type(system), intent(in) :: MySys
    type(psi), intent(in) :: ket

    character(len=256) :: outF
    character(len=40) :: fileformat
    integer :: i, j, k, ierr, nwritten

    fileformat = '(1(1X,E17.11))'

    write(outF, '(A,I0.4,A,I0.4,A,I0.4)') 'wfc_', ket%ispin, '-', ket%ikpts, '-', ket%iband
    open(unit=23, file=trim(outF), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) 'Error creating wfc files'
    end if

    write(23,'(A,3I5)') '# ', MySys%NGPTAR(1), MySys%NGPTAR(2), MySys%NGPTAR(3)
    nwritten = 0
    do k = 1, MySys%NGPTAR(3)
      do j = 1, MySys%NGPTAR(2)
        do i = 1, MySys%NGPTAR(1)
          nwritten=nwritten+1
          if ( mod(nwritten,10)==0 ) then
             write(23, fileformat) real(wfc_r(i,j,k))
          else
             write(23, fileformat, advance='no') real(wfc_r(i,j,k))
          endif
        end do
      end do
    end do

    ! write(23,'(/,A,3I5)') '# ', MySys%NGPTAR(1), MySys%NGPTAR(2), MySys%NGPTAR(3)
    ! nwritten = 0
    ! do k = 1, MySys%NGPTAR(3)
    !   do j = 1, MySys%NGPTAR(2)
    !     do i = 1, MySys%NGPTAR(1)
    !       nwritten=nwritten+1
    !       if ( mod(nwritten,10)==0 ) then
    !          write(23, fileformat) aimag(wfc_r(i,j,k))
    !       else
    !          write(23, fileformat, advance='no') aimag(wfc_r(i,j,k))
    !       endif
    !     end do
    !   end do
    ! end do

    close(23)
  end subroutine outwfc
end module vaspwfc

