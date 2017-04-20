module gvector
  use prec
  use lattice
  use wavecar
  use constants

  implicit none
  
  contains

    subroutine GCUTS(MySys)
      use info
      implicit none
      type(system), intent(inout) :: MySys

      real(kind=q) :: CUTOF(3)
      integer :: i, ng

      ! logical :: fftchk_furth

      do i = 1, 3
        ! Gcut^2 / 2 = ENCUT
        CUTOF(i) =SQRT(MySys%ENCUT / RYTOEV)/(TPI/(MySys%MyLatt%ANORM(i)/AUTOA))
        ! twice Gcut incorporate all the points within the sphere
        ! MySys%NGPTAR(i) = 2 * NINT(CUTOF(i)) + 1
        ng = 2 * (2 * NINT(CUTOF(i)) + 1)
        ! vasp use even number in NGX, NGY, NGZ
        if (mod(ng, 2) /= 0) then
          ng = ng + 1
        end if
        MySys%NGPTAR(i) = max(ng, MySys%NGPTAR(i))

        ! choose a valid value for FFT
        do
          if (fftchk_furth(MySys%NGPTAR(i))) then
            exit
          else
            MySys%NGPTAR(i) = MySys%NGPTAR(i) + 1
          end if
        end do
      end do

      allocate(MySys%LPCTX(MySys%NGPTAR(1)))
      allocate(MySys%LPCTY(MySys%NGPTAR(2)))
      allocate(MySys%LPCTZ(MySys%NGPTAR(3)))

      ! G-vector index stored in LPCT[XYZ]
      call GINDEX(MySys%NGPTAR(1), MySys%NGPTAR(2), MySys%NGPTAR(3), & 
                  MySys%LPCTX, MySys%LPCTY, MySys%LPCTZ)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! G-vector indexs within Gcut
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine GINDEX(NGX, NGY, NGZ, LPCTX, LPCTY, LPCTZ)
      implicit none

      integer, intent(in) :: NGX, NGY, NGZ
      integer, intent(inout) :: LPCTX(NGX), LPCTY(NGY), LPCTZ(NGZ)
      integer i

      do i=1, NGX
        if (i > (NGX/2+1)) then
          LPCTX(i) = i-1-NGX
        else
          LPCTX(i) = i - 1
        end if
      end do

      do i=1, NGY
        if (i > (NGY/2+1)) then
          LPCTY(i) = i-1-NGY
        else
          LPCTY(i) = i - 1
        end if
      end do

      do i=1, NGZ
        if (i > (NGZ/2+1)) then
          LPCTZ(i) = i-1-NGZ
        else
          LPCTZ(i) = i - 1
        end if
      end do

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! count number of plane waves of specific k-point
    ! parallel gamma version
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine CNT_PLW_NUM(MySys)
      use info
      implicit none

      type(system), intent(inout) :: MySys
      integer :: i, j, k, n, itmp
      integer :: NGX, NGY, NGZ
      real(q) :: G1, G2, G3, GX, GY, GZ, ENERG

      NGX = MySys%NGPTAR(1)
      NGY = MySys%NGPTAR(2)
      NGZ = MySys%NGPTAR(3)

      allocate(MySys%GX(MySys%MAXPLWS, MySys%NKPTS), &
               MySys%GY(MySys%MAXPLWS, MySys%NKPTS), &
               MySys%GZ(MySys%MAXPLWS, MySys%NKPTS), &
               MySys%grid_index(3, MySys%MAXPLWS, MySys%NKPTS))

      do n = 1, MySys%NKPTS
        itmp = 0
        do k = 1, NGZ
          do j = 1, NGY
            do i = 1, NGX
               ! exclude some component of parallel gamma version (wNGZhalf)
               ! (C(G) = C*(-G))
               if (MySys%LPCTZ(k)<0) cycle
               if (MySys%LPCTZ(k)==0 .AND. MySys%LPCTY(j)<0) cycle
               if (MySys%LPCTZ(k)==0 .AND. MySys%LPCTY(j)==0 .AND. MySys%LPCTX(i)<0) cycle
               
               G1=(MySys%LPCTX(i)+MySys%VKPTS(1,n))
               G2=(MySys%LPCTY(j)+MySys%VKPTS(2,n))
               G3=(MySys%LPCTZ(k)+MySys%VKPTS(3,n))

               GX= (G1*MySys%MyLatt%B(1,1)+G2*MySys%MyLatt%B(1,2)+G3*MySys%MyLatt%B(1,3))*TPI
               GY= (G1*MySys%MyLatt%B(2,1)+G2*MySys%MyLatt%B(2,2)+G3*MySys%MyLatt%B(2,3))*TPI
               GZ= (G1*MySys%MyLatt%B(3,1)+G2*MySys%MyLatt%B(3,2)+G3*MySys%MyLatt%B(3,3))*TPI

               ENERG =HSQDTM*((GX**2)+(GY**2)+(GZ**2))

               if (ENERG < MySys%ENCUT) then
                 itmp = itmp + 1
                 MySys%GX(itmp,n) = MySys%LPCTX(i)
                 MySys%GY(itmp,n) = MySys%LPCTY(j)
                 MySys%GZ(itmp,n) = MySys%LPCTZ(k)

                 MySys%grid_index(1,itmp,n) = i
                 MySys%grid_index(2,itmp,n) = j
                 MySys%grid_index(3,itmp,n) = k
               end if
            end do
          end do
        end do
        if (MySys%NPLWS(n) /= itmp) then
          write(*,*) "Error! Plane waves number not equal!"
        end if
      end do

    end subroutine

    ! check whether the NIN is valid for FFT
    logical function fftchk_furth(nin)
      use prec
      implicit none

      integer, intent(in) :: nin

      integer, parameter :: nfact = 4
      integer, parameter, dimension(nfact) :: ifact = (/2, 3, 5, 7/)
      integer, dimension(nfact) :: ncount
      integer :: i, n, next

      n=nin
      do i=1,nfact
        ncount(i) = 0
        do
          next = n / ifact(i)
          if (next*ifact(i) == n) then
            n = next
            ncount(i) = ncount(i) + 1
          else
            exit
          end if
        end do
      end do

      if (n==1 .and. (ncount(1)/=0)) then
        fftchk_furth=.true.
      else
        fftchk_furth=.false.
      endif

      return
    end function fftchk_furth
    

end module gvector
