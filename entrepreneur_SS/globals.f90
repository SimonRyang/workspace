module globals

  use toolbox
  use sorting

  implicit none

  ! number of years the household lives
  integer, parameter :: JJ = 16

  ! number of years a household can be a worker (+1)
  integer, parameter :: JR = 10

  ! number of years a household can be an entrepreneur (+1)
  integer, parameter :: JE = 13

  ! number of permanent skill classes
  integer, parameter :: NS = 3

  ! number of transitory shock process values (worker)
  integer, parameter :: NW = 5

  ! number of transitory shock process values (entrepreneur)
  integer, parameter :: NE = 5

  ! number of points on the asset grid (-1)
  integer, parameter :: NA = 24

  ! number of points on the pension claim grid (-1)
  integer, parameter :: NP = 5

  ! number of occupations (-1)
  integer, parameter :: NO = 1

  ! household parameters
  real*8 :: gamma, sigma, mu_b, beta, l_bar, swc

  ! production parameters
  real*8 :: alpha, delta, nu

  ! numerical parameters
  real*8 :: a_l, a_u, a_grow
  real*8 :: p_l, p_u
  real*8 :: damp, tol
  integer :: itermax

  ! macroeconomic variables
  real*8:: gy, by
  real*8 :: r, w, inc_pen, inc_tax, pinv
  real*8 :: KK, KC, KE, AA, LC, HH
  real*8 :: YY, YC, YE, CC, CCX, II, GG, NEX
  real*8 :: BB, BF, BQ
  real*8 :: TAc, TAr, TAw, TAy

  ! government variables
  real*8 :: tauc, taur, tauy, taup
  real*8 :: kappa, sscc, lambda, phi, mu
  real*8 :: pen(0:NP, JJ), PP, PE, PRE, PC, BP = 0d0

  ! progressive income tax
  real*8, parameter :: t1 = 0.14d0, t2 = 0.24d0, t3 = 0.45d0
  real*8 :: r1, r2, r3, b1, b2

  ! cohort aggregate variables
  real*8 :: pop_w(NS), pop_e(NS), pop_r(NS), pop_re(NS)
  real*8 :: bqs(NS), beq(NS, JJ)
  real*8 :: vv_coh(JJ) = 0d0
  real*8 :: c_coh(0:1, JJ), a_coh(0:1, JJ), k_coh(JJ)
  real*8 :: inc_coh(0:1, JJ), o_coh(0:1, 0:1, JJ), os_coh(0:1, 0:1, NS, JJ), flc_coh(JJ)

  ! the shock process
  real*8 :: dist_skill(NS)
  real*8 :: eta(NW, NS), dist_eta(NW, NS), pi_eta(NW, NW, NS)
  real*8 :: theta(NE, NS), dist_theta(NE, NS), pi_theta(NE, NE, NS)

  ! demographic and other model parameters
  real*8 :: eff(JJ, NS), rpop(NS, JJ), pop(JJ), psi(NS, JJ+1), Gama(JJ), n_p

  ! individual variables
  real*8 :: a(0:NA), p(0:NP)
  real*8, allocatable :: aplus(:, :, :, :, :, :, :)
  real*8, allocatable :: pplus(:, :, :, :, :, :, :)
  real*8, allocatable :: c(:, :, :, :, :, :, :)
  real*8, allocatable :: cx(:, :, :, :, :, :, :)
  real*8, allocatable :: l(:, :, :, :, :, :, :)
  real*8, allocatable :: k(:, :, :, :, :, :, :)
  real*8, allocatable :: oplus(:, :, :, :, :, :, :)
  real*8, allocatable :: pencon(:, :, :, :, :, :, :)
  real*8, allocatable :: inctax(:, :, :, :, :, :, :)
  real*8, allocatable :: captax(:, :, :, :, :, :, :)
  real*8, allocatable :: VV(:, :, :, :, :, :, :)
  real*8, allocatable :: EV(:, :, :, :, :, :, :)
  real*8, allocatable :: m(:, :, :, :, :, :, :)

  ! numerical variables
  integer :: io_com, ia_com, ip_com, iw_com, ie_com, is_com, ij_com
  real*8 :: c_com, cx_com, l_com, k_com, pplus_com, oplus_com, pencon_com, inctax_com, captax_com, DIFF

  ! statistical variables
  logical :: gini_on = .false.
  real*8 :: gini_w, percentiles_w(6), gini_i, percentiles_i(6)

  ! switches
  logical :: smopec = .false.   ! .true. = economcy is smopec
  logical :: ent = .true.     ! .true. = endogenous decision to become an entrepreneur

  !$omp threadprivate(io_com, ia_com, ip_com, iw_com, ie_com, is_com, ij_com)
  !$omp threadprivate(c_com, cx_com, l_com, k_com, pplus_com, oplus_com, pencon_com, inctax_com, captax_com)

contains


  !##############################################################################
  ! FUNCTION valuefunc_w
  !
  ! Determines the value function of a worker
  !##############################################################################
  function valuefunc_w(xy)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: xy(:)
    real*8 :: valuefunc_w

    !##### OTHER VARIABLES ####################################################
    real*8 :: a_plus, wage, c_help, valuefunc_help

    ! tomorrow's assets
    a_plus = xy(1)

    ! today's labor
    l_com  = xy(2)

    ! calculate the wage rate and next periods pension claims
    wage = w*eff(ij_com, is_com)*eta(iw_com, is_com)
    pplus_com = (p(ip_com)*dble(ij_com-1) + mu*(lambda &
           + (1d0-lambda)*min(wage*l_com/inc_pen, sscc)))/dble(ij_com)

    ! worker do not invest
    k_com = 0d0

    ! calculate contribution to pension system
    pencon_com = taup*min(wage*l_com, sscc*inc_pen)

    ! calculate income tax
    inctax_com = tarif(max(wage*l_com - 0.08d0*wage*l_com - 0.04d0*inc_tax - pencon_com, 0d0) + pen(ip_com, ij_com))

    ! calculate capital gains tax
    captax_com = taur*1.055d0*max(r*a(ia_com) - 0.08d0*r*a(ia_com) - 2d0*0.0267d0*inc_tax, 0d0)

    ! calculate consumption
    c_com = ((1d0+r)*a(ia_com) + wage*l_com + beq(is_com, ij_com) + pen(ip_com, ij_com) &
             - pencon_com - inctax_com - captax_com - a_plus)*pinv
    c_help = ((1d0+r)*a(ia_com) + wage*l_com + beq(is_com, ij_com) + pen(ip_com, ij_com) &
              - pencon_com - inctax_com - captax_com - a_plus - swc)*pinv

    ! calculate tomorrow's part of the value function and occupational decision
    valuefunc_w = 0d0
    valuefunc_help = 0d0
    cx_com = 0d0
    oplus_com = 0d0

    if (ij_com < JJ) then

      ! interpolate next period's value function as a worker/retiree
      valuefunc_w = util(c_com, l_com) + beta*psi(is_com, ij_com+1)*interpolate_EV(0, a_plus, pplus_com, iw_com, ie_com, is_com, ij_com+1)**(1d0-gamma)/(1d0-gamma)

      ! interpolate next period's value function as an entrepreneur
      if (ij_com < JR-1) then

        valuefunc_help = util(c_help, l_com) + beta*psi(is_com, ij_com+1)*interpolate_EV(1, a_plus, pplus_com, iw_com, ie_com, is_com, ij_com+1)**(1d0-gamma)/(1d0-gamma)

        ! set next period's occupational decision
        if (valuefunc_help > valuefunc_w .and. ent) then
          valuefunc_w = valuefunc_help
          c_com = c_help
          cx_com = swc
          oplus_com = 1d0
        endif

      endif

    endif

    ! add today's part and discount
    valuefunc_w = -valuefunc_w

  end function


  !##############################################################################
  ! FUNCTION valuefunc_e
  !
  ! Determines the value function of an entrepreneur
  !##############################################################################
  function valuefunc_e(xy)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: xy(:)
    real*8 :: valuefunc_e

    !##### OTHER VARIABLES ####################################################
    real*8 :: a_plus, profit, c_help, valuefunc_help

    ! tomorrow's assets
    a_plus = xy(1)

    ! today's investment
    k_com = xy(2)

    ! no investment without any assets
    if (ia_com == 0) k_com = 0d0

    ! today's fixed labor
    l_com = l_bar

    ! entrepreneur's profit
    profit = theta(ie_com, is_com)*(k_com**alpha*(eff(ij_com, is_com)*eta(iw_com, is_com)*l_bar)**(1d0-alpha))**nu &
             - delta*k_com - r*max(k_com-a(ia_com), 0d0)

    ! calculate contribution to pension system
    if (ij_com < JR) then
        pencon_com = phi*taup*min(profit, sscc*inc_pen)
    else
        pencon_com = 0d0
    endif

    ! calculate income tax
    inctax_com = tarif(max(profit - 0.08d0*profit - 0.04d0*inc_tax - pencon_com, 0d0) + pen(ip_com, ij_com))

    ! calcualte capital gains tax
    captax_com = taur*1.055d0*max(r*max(a(ia_com)-k_com, 0d0) - 0.08d0*r*a(ia_com) - 2d0*0.0267d0*inc_tax, 0d0)

    ! calculate consumption
    c_com =  (a(ia_com) + r*max(a(ia_com)-k_com, 0d0) + profit + beq(is_com, ij_com) + pen(ip_com, ij_com)  &
           - captax_com - inctax_com - pencon_com - a_plus - swc)*pinv
    c_help =  (a(ia_com) + r*max(a(ia_com)-k_com, 0d0) + profit + beq(is_com, ij_com) + pen(ip_com, ij_com)  &
           - captax_com - inctax_com - pencon_com - a_plus)*pinv

    ! calculate next periods pension claims
    if (ij_com < JR) then
      pplus_com = (p(ip_com)*dble(ij_com-1) + mu*phi*(lambda &
             + (1d0-lambda)*min(profit/inc_pen, sscc)))/dble(ij_com)
    else
      pplus_com = p(ip_com)
    endif

    ! calculate tomorrow's part of the value function and occupational decision
    valuefunc_e = 0d0
    valuefunc_help = 0d0
    cx_com = swc
    oplus_com = 0d0

    if (ij_com < JJ) then

      ! interpolate next period's value function as a worker/retiree
      valuefunc_e = util(c_com, l_com) + beta*psi(is_com, ij_com+1)*interpolate_EV(0, a_plus, pplus_com, iw_com, ie_com, is_com, ij_com+1)**(1d0-gamma)/(1d0-gamma)

      ! interpolate next period's value function as an entrepreneur
      if (ij_com < JE-1) then

        valuefunc_help = util(c_help, l_com) + beta*psi(is_com, ij_com+1)*interpolate_EV(1, a_plus, pplus_com, iw_com, ie_com, is_com, ij_com+1)**(1d0-gamma)/(1d0-gamma)

        ! set next period's occupational decision
        if (valuefunc_help > valuefunc_e .and. ent) then
          valuefunc_e = valuefunc_help
          c_com = c_help
          cx_com = 0d0
          oplus_com = 1d0
        endif

      endif

    endif

    ! add today's part and discount
    if (ij_com >= JR) then
      valuefunc_e = -(valuefunc_e + (1d0-psi(is_com, ij_com+1))*mu_b*a_plus**(1d0-gamma))
    else
      valuefunc_e = -valuefunc_e
    endif


  end function


  !##############################################################################
  ! FUNCTION valuefunc_r
  !
  ! Determines the value function of a retired individual
  !##############################################################################
  function valuefunc_r(xy)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: xy
    real*8 :: valuefunc_r

    !##### OTHER VARIABLES ####################################################
    real*8 :: a_plus, varphi
    integer :: ial, iar

    ! tomorrow's assets
    a_plus = xy

    ! today's labor
    l_com  = 0d0

    ! calculate the wage rate and next periods pension claims
    pplus_com = p(ip_com)

    ! retirees do not invest
    k_com = 0d0

    ! calculate contribution to pension system
    pencon_com = 0d0

    ! calculate income tax
    inctax_com = tarif(pen(ip_com, ij_com))

    ! calculate capital gains tax
    captax_com = taur*1.055d0*max(r*a(ia_com) - 0.08d0*r*a(ia_com) - 2d0*0.0267d0*inc_tax, 0d0)

    ! calculate consumption
    c_com = ((1d0+r)*a(ia_com) + beq(is_com, ij_com) + pen(ip_com, ij_com) &
         - inctax_com - captax_com - a_plus)*pinv

    ! calculate tomorrow's part of the value function and occupational decision
    valuefunc_r = 0d0
    oplus_com = 0d0

    if (ij_com < JJ) then

      ! interpolate next period's value function as a worker/retiree
      valuefunc_r = interpolate_EV(0, a_plus, pplus_com, iw_com, ie_com, is_com, ij_com+1)**(1d0-gamma)/(1d0-gamma)

    endif

    ! add today's part and discount
    valuefunc_r = -(util(c_com, l_com) + beta*psi(is_com, ij_com+1)*valuefunc_r &
                    + (1d0-psi(is_com, ij_com+1))*mu_b*a_plus**(1d0-gamma))

  end function



  !##############################################################################
  ! FUNCTION profent
  !
  ! Computes profit of an entrepreneur
  !##############################################################################
  function profent(k, ij, ia, is, iw, ie)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: k
    integer, intent(in) :: ij, ia, is, iw, ie
    real*8 :: profent

    ! compute profit
    profent = theta(ie, is)*(k**alpha*(eff(ij, is)*eta(iw, is)*l_bar)**(1d0-alpha))**nu - delta*k - r*max(k-a(ia), 0d0)

  end function


  !##############################################################################
  ! FUNCTION tarif
  !
  ! Calculates income tax
  !##############################################################################
  function tarif(zn)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: zn
    real*8 :: tarif

    !##### OTHER VARIABLES ####################################################
    real*8 :: tr2, tr3, znr

    ! determine parameters of tarif
    tr2 =     (r2-r1)*t1 + (r2-r1)**2d0*b1*0.5d0
    tr3 = tr2 + (r3-r2)*t2 + (r3-r2)**2d0*b2*0.5d0

    ! apply german income tax tarif
    if (zn <= r1) then
      tarif= 0.0d0
    elseif (zn <= r2) then
      znr = zn - r1
      tarif = znr*t1 + znr**2d0*b1*0.5d0
    elseif (zn <= r3) then
      znr = zn - r2
      tarif = tr2 + znr*t2 + znr**2d0*b2*0.5d0
    else
      znr = zn - r3
      tarif = tr3 + znr*t3
    endif

  end function


  !##############################################################################
  ! FUNCTION util
  !
  ! Determines instantaneous utility
  !##############################################################################
  function util(cons, lab)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: cons, lab
    real*8 :: util

    ! determine utility
    if (cons <= 0d0 .or. lab < 0d0 .or. lab >= 1d0) then
      if(cons <= 0d0)util=-1d18*(1d0+abs(cons))
      if(lab < 0d0)util=-1d18*(1d0+abs(lab))
      if(lab >= 1d0)util=-1d18*lab
    else
      util = (cons**sigma*(1d0-lab)**(1d0-sigma))**(1d0-gamma)/&
          (1d0-gamma)
    endif

  end function

  !##############################################################################
  ! FUNCTION margu
  !
  ! Determines marginal utility
  !##############################################################################
  function margu(cons, lab)

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: cons, lab
    real*8 :: margu

    ! determine marginal utility
    margu = sigma*(cons**sigma*(1d0-lab)**(1d0-sigma))**(1d0-gamma)/((1d0+tauc)*cons)

  end function


  !##############################################################################
  ! FUNCTION interpolate_EV
  !
  ! Interpolates the expected valuefunction
  !##############################################################################
  function interpolate_EV(io, a_plus, p_plus, iw, ie, is, ij)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: a_plus, p_plus
    integer, intent(in) :: io, iw, ie, is, ij
    real*8 :: interpolate_EV

    !##### OTHER VARIABLES ####################################################
    real*8 :: varphi, varpsi
    integer :: ial, iar, ipl, ipr

    ! interpolate value function
    call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)
    call linint_Equi(p_plus, p_l, p_u, NP, ipl, ipr, varpsi)

    interpolate_EV = (varphi*varpsi*EV(io, ial, ipl, iw, ie, is, ij) &
                      + varphi*(1d0-varpsi)*EV(io, ial, ipr, iw, ie, is, ij) &
                      + (1d0-varphi)*varpsi*EV(io, iar, ipl, iw, ie, is, ij) &
                      + (1d0-varphi)*(1d0-varpsi)*EV(io, iar, ipr, iw, ie, is, ij))

  end function


  !##############################################################################
  ! FUNCTION check_grid
  !
  ! Checks for the maximum gridpoint used
  !##############################################################################
  subroutine check_grid(iamax)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer :: iamax(JJ)

    !##### OTHER VARIABLES ####################################################
    integer :: ia, ip, iw, ie, is, ij

    iamax = 0

    do ij = 1, JJ

      ! check for the maximum asset grid point used at a certain age
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ia = 0, NA
                  if (m(0, ia, ip, iw, ie, is, ij) > 0d0) iamax(ij)=ia
                  if (m(1, ia, ip, iw, ie, is, ij) > 0d0) iamax(ij)=ia
              enddo ! ix
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is

    enddo ! ij


  end subroutine


  !##############################################################################
  ! FUNCTION sigma5
  !
  ! Converts sigma of an one-year process into sigma of an five-year process
  !##############################################################################
  function sigma5(rho, sigma)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: rho, sigma
    real*8 :: sigma5

    !##### OTHER VARIABLES ####################################################
    real*8 :: Gamma5
    integer :: i, j

    Gamma5 = 0d0
    do i = 0, 4
      do j = 0, 4
        Gamma5 = Gamma5 + rho**abs(j-i)
      enddo
    enddo
    Gamma5 = Gamma5/25d0

    sigma5 = Gamma5*sigma*(1d0-rho**10d0)/(1d0-rho**2d0)

  end function

  !##############################################################################
  ! FUNCTION gini_wealth
  !
  ! Calculates GINI coefficient of wealth
  !##############################################################################
  subroutine gini_wealth()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: ia, ip, iw, ie, is, ij, ic, ICMAX
    real*8, allocatable :: xs(:), ys(:), xcum(:), ycum(:)

    if(allocated(xs))deallocate(xs)
    if(allocated(ys))deallocate(ys)
    if(allocated(xcum))deallocate(xcum)
    if(allocated(ycum))deallocate(ycum)

    allocate(xs(2*NS*NE*NW*(NP+1)*(NA+1)*JJ))
    allocate(ys(2*NS*NE*NW*(NP+1)*(NA+1)*JJ))
    allocate(xcum(2*NS*NE*NW*(NP+1)*(NA+1)*JJ))
    allocate(ycum(2*NS*NE*NW*(NP+1)*(NA+1)*JJ))

    ic = 1

    do ij = 1, JJ
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ia = 0, NA

                if(m(0, ia, ip, iw, ie, is, ij) > 0d0) then
                  ys(ic) = m(0, ia, ip, iw, ie, is, ij)
                  xs(ic) = a(ia)
                  ic = ic + 1
                endif
                if(m(1, ia, ip, iw, ie, is, ij) > 0d0) then
                  ys(ic) = m(1, ia, ip, iw, ie, is, ij)
                  xs(ic) = a(ia)
                  ic = ic + 1
                endif

              enddo ! ia
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is
    enddo ! ij

    ! get array size and normalize ys
    ICMAX = ic - 1
    ys(1:ICMAX) = ys(1:ICMAX)/sum(ys(1:ICMAX))

    ! sort array
    call quick_sort(xs(1:ICMAX), ys(1:ICMAX))

    ! calculate cumulative distributions
    xcum(1) = 0d0
    ycum(1) = 0d0
    do ic = 2, ICMAX+1
      xcum(ic) = xcum(ic-1) + xs(ic-1)*ys(ic-1)
      ycum(ic) = ycum(ic-1) + ys(ic-1)
    enddo

    gini_w = 0d0
    do ic = 2, ICMAX+1
      gini_w = gini_w + ys(ic-1)*(xcum(ic-1)+xcum(ic))
    enddo
    gini_w = 1d0-gini_w/xcum(ICMAX+1)

    percentiles_w = 0d0
    do ic = ICMAX+1, 2, -1
      if (1d0-ycum(ic) > 0.01 .and. percentiles_w(1) <= 0d0) then
        percentiles_w(1) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.05 .and. percentiles_w(2) <= 0d0) then
        percentiles_w(2) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.10 .and. percentiles_w(3) <= 0d0) then
        percentiles_w(3) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.20 .and. percentiles_w(4) <= 0d0) then
        percentiles_w(4) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.40 .and. percentiles_w(5) <= 0d0) then
        percentiles_w(5) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.60 .and. percentiles_w(6) <= 0d0) then
        percentiles_w(6) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
    enddo

  end subroutine


  !##############################################################################
  ! FUNCTION gini_income
  !
  ! Calculates GINI coefficient of income
  !##############################################################################
  subroutine gini_income()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: is, ie, iw, ip, ia, ij, ic, ICMAX
    real*8, allocatable :: xs(:), ys(:), xcum(:), ycum(:)

    if(allocated(xs))deallocate(xs)
    if(allocated(ys))deallocate(ys)
    if(allocated(xcum))deallocate(xcum)
    if(allocated(ycum))deallocate(ycum)

    allocate(xs(2*NS*NE*NW*(NP+1)*(NA+1)*JJ))
    allocate(ys(2*NS*NE*NW*(NP+1)*(NA+1)*JJ))
    allocate(xcum(2*NS*NE*NW*(NP+1)*(NA+1)*JJ))
    allocate(ycum(2*NS*NE*NW*(NP+1)*(NA+1)*JJ))

    ic = 1

    do ij = 1, JJ
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ia = 0, NA

                if(m(0, ia, ip, iw, ie, is, ij) > 0d0) then
                  ys(ic) = m(0, ia, ip, iw, ie, is, ij)
                  xs(ic) = a(ia)*r + pen(ip, ij) + eff(ij, is)*eta(iw, is)*l(0, ia, ip, iw, ie, is, ij)*w
                  ic = ic + 1
                endif
                if(m(1, ia, ip, iw, ie, is, ij) > 0d0) then
                  ys(ic) = m(1, ia, ip, iw, ie, is, ij)
                  xs(ic) = max(a(ia)-k(1, ia, ip, iw, ie, is, ij), 0d0)*r + pen(ip, ij) &
                           + profent(k(1, ia, ip, iw, ie, is, ij), ij, ia, is, iw, ie)
                  ic = ic + 1
                endif

              enddo ! ia
            enddo ! ip
          enddo ! iw
        enddo ! iw
      enddo ! is
    enddo ! ij


    ! get array size and normalize ys
    ICMAX = ic - 1
    ys(1:ICMAX) = ys(1:ICMAX)/sum(ys(1:ICMAX))

    ! sort array
    call quick_sort(xs(1:ICMAX), ys(1:ICMAX))

    ! calculate cumulative distributions
    xcum(1) = 0d0
    ycum(1) = 0d0
    do ic = 2, ICMAX+1
      xcum(ic) = xcum(ic-1) + xs(ic-1)*ys(ic-1)
      ycum(ic) = ycum(ic-1) + ys(ic-1)
    enddo

    gini_i = 0d0
    do ic = 2, ICMAX+1
      gini_i = gini_i + ys(ic-1)*(xcum(ic-1)+xcum(ic))
    enddo
    gini_i = 1d0-gini_i/xcum(ICMAX+1)

    percentiles_i = 0d0
    do ic = ICMAX+1, 2, -1
      if (1d0-ycum(ic) > 0.01d0 .and. percentiles_i(1) <= 0d0) then
        percentiles_i(1) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.05d0 .and. percentiles_i(2) <= 0d0) then
        percentiles_i(2) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.10d0 .and. percentiles_i(3) <= 0d0) then
        percentiles_i(3) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.20d0 .and. percentiles_i(4) <= 0d0) then
        percentiles_i(4) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.40d0 .and. percentiles_i(5) <= 0d0) then
        percentiles_i(5) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.60d0 .and. percentiles_i(6) <= 0d0) then
        percentiles_i(6) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
    enddo

  end subroutine


end module
