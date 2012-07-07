module math

contains

!===============================================================================
! WIGNER samples a Wigner distribution of energy level spacings. Note that this
! scheme is C50 in the Monte Carlo Sampler from Los Alamos (LA-9721-MS).
!===============================================================================

  function wigner(D_avg) result (D)

    real(8), intent(in) :: D_avg ! average level spacing
    real(8)             :: D     ! sampled level spacing

    real(8) :: c

    c = -4.*D_avg*D_avg/PI * log(prn())
    D = sqrt(c)

  end function wigner

!===============================================================================
! CHI_SQUARED samples a chi-squared distribution with n degrees of freedom. The
! distribution of resonance widths in the unresolved region is given by a
! chi-squared distribution. For the special case of n=1, this is a Porter-Thomas
! distribution. For cases with n odd, rule C64 is used whereas for cases with n
! even, rule C45 is used.
!===============================================================================

  function chi_squared(n, G_avg) result(G)

    integer, intent(in)           :: n     ! number of degrees of freedom
    real(8), intent(in), optional :: G_avg ! average resonance width

    integer :: i       ! loop index
    real(8) :: G       ! sampled random variable (or resonance width)
    real(8) :: x, y, c ! dummy variables
    real(8) :: r1, r2  ! psuedorandom numbers

    select case (mod(n,2))
    case (0)
       ! Even number of degrees of freedom can be sampled via rule C45. We can
       ! sample x as -2/n*log(product(r_i, i = 1 to n/2))
       x = ONE
       do i = 1, n/2
          x = x * prn()
       end do
       x = -2./n * log(x)

    case (1)
       ! Odd number of degrees of freedom can be sampled via rule C64. We can
       ! sample x as -2/n*(log(r)*cos^2(pi/2*r) + log(product(r_i, i = 1 to
       ! floor(n/2)))

       ! Note that we take advantage of integer division on n/2
       y = ONE
       do i = 1, n/2
          y = y * prn()
       end do

       r1 = prn()
       r2 = prn()
       c = cos(PI/2.*r2)
       x = -2./n * (log(y) + log(r1)*c*c)
    end select

    ! If sampling a chi-squared distribution for a resonance width and the
    ! average resonance width has been given, return the sampled resonance
    ! width.
    if (present(G_avg)) then
       G = x * G_avg
    else
       G = x
    end if

  end function chi_squared

end module math
