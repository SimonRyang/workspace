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

  ! number of transition periods
  integer, parameter :: TT = 2

  ! number of permanent skill classes
  integer, parameter :: NS = 3

  ! number of transitory shock process values (worker)
  integer, parameter :: NW = 5

  ! number of transitory shock process values (entrepreneur)
  integer, parameter :: NE = 5

  ! number of points on the asset grid (-1)
  integer, parameter :: NA = 512

  ! number of points on the annuitized asset grid (-1)
  integer, parameter :: NX = 1

  ! number of points on the pension claim grid (-1)
  integer, parameter :: NP = 6

  ! number of occupations (-1)
  integer, parameter :: NO = 1

  ! household parameters
  real*8 :: gamma, sigma, beta, l_bar, phi1, phi2

  ! production parameters
  real*8 :: alpha, delta, nu, suc

  ! numerical parameters
  real*8 :: a_l, a_u, a_grow
  real*8 :: x_l, x_u, x_grow
  real*8 :: p_l, p_u
  real*8 :: damp, tol
  integer :: itermax

  ! counter variables
  integer :: iter

  ! macroeconomic variables
  real*8:: gy, by
  real*8 :: r(0:TT), w(0:TT), inc_bar(0:TT), psix(NS, JJ, 0:TT), pinv(0:TT)
  real*8 :: KK(0:TT), KC(0:TT), KE(0:TT), AA(0:TT), AX(0:TT), XB(0:TT), LC(0:TT), HH(0:TT)
  real*8 :: YY(0:TT), YC(0:TT), YE(0:TT), CC(0:TT), II(0:TT), GG(0:TT), NEX(0:TT)
  real*8 :: BB(0:TT), BF(0:TT), BQ(0:TT)
  real*8 :: TAc(0:TT), TAr(0:TT), TAw(0:TT), TAy(0:TT)

  ! government variables
  real*8 :: tauc(0:TT), taur(0:TT), tauy(0:TT), taup(0:TT)
  real*8 :: kappa(0:TT), sscc(0:TT), lambda(0:TT), phi(0:TT), mu(0:TT)
  real*8 :: pen(0:NP, JJ, 0:TT), PP(0:TT), PE(0:TT), PRE(0:TT), PC(0:TT), BP(0:TT) = 0d0

  ! LSRA variables
  real*8 :: BA(0:TT) = 0d0, SV(0:TT) = 0d0, lsra_comp, lsra_all, Vstar
  logical :: lsra_on

  ! progressive income tax
  real*8, parameter :: t1 = 0.14d0, t2 = 0.24d0, t3 = 0.45d0
  real*8 :: r1, r2, r3, b1, b2

  ! cohort aggregate variables
  real*8 :: pop_w(NS, 0:TT), pop_e(NS, 0:TT), pop_r(NS, 0:TT), pop_re(NS, 0:TT)
  real*8 :: bqs(NS, 0:TT)
  real*8 :: x_coh(JJ, 0:TT), bx_coh(JJ, 0:TT), os_coh(0:1, 0:1, NS, JJ, 0:TT)
  real*8 :: vv_coh(JJ, 0:TT) = 0d0

  ! the shock process
  real*8 :: dist_skill(NS)
  real*8 :: eta(NW, NS), dist_eta(NW, NS), pi_eta(NW, NW, NS)
  real*8 :: theta(NE, NS), dist_theta(NE, NS), pi_theta(NE, NE, NS)

  ! demographic and other model parameters
  real*8 :: eff(JJ, NS), rpop(NS, JJ, 0:TT), pop(JJ, 0:TT), psi(NS, JJ+1), beq(NS, JJ, 0:TT), Gama(JJ), n_p

  ! individual variables
  real*8 :: a(0:NA), x(0:NX), p(0:NP)
  real*8, allocatable :: aplus(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: xplus(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: pplus(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: c(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: l(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: k(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: mx(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: oplus(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: pencon(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: inctax(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: captax(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: VV(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: EV(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: m(:, :, :, :, :, :, :, :, :)
  real*8, allocatable :: v(:, :, :, :, :, :, :, :, :)

  ! numerical variables
  integer :: io_com, ia_com, ix_com, ip_com, iw_com, ie_com, is_com, ij_com, it_com
  real*8 :: c_com, l_com, k_com, mx_com, xplus_com, pplus_com, oplus_com, pencon_com, inctax_com, captax_com, DIFF(0:TT)

  ! statistical variables
  logical :: gini_on = .false.
  real*8 :: gini_w(0:TT), percentiles_w(6, 0:TT), gini_i(0:TT), percentiles_i(6, 0:TT)

  ! switches
  logical :: smopec = .false.   ! .true. = economcy is smopec
  logical :: ann = .false.       ! .true. = wealth can be annuitized
  logical :: ent = .true.       ! .true. = endogenous decision to become an entrepreneur
  logical :: labor = .true.     ! .true. = endogenous labor decision of worker
  logical :: pen_debt = .false. ! .true. = pension system can run into debts

  !$omp threadprivate(io_com, ia_com, ix_com, ip_com, iw_com, ie_com, is_com, ij_com, it_com)
  !$omp threadprivate(c_com, l_com, k_com, mx_com, xplus_com, pplus_com, oplus_com, pencon_com, inctax_com, captax_com)

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
    real*8 :: a_plus, wage, v_ind, valuefunc_help
    integer :: itp

    ! tomorrow's assets
    a_plus = xy(1)

    ! today's labor
    if (labor) then
      l_com  = xy(2)
    else
      l_com = l(0, ia_com, ix_com, ip_com, iw_com, ie_com, is_com, ij_com, 0)
    endif

    ! today's investment in annuitized assets
    mx_com = 0d0
    if (ij_com == JR-1 .and. ann) mx_com = xy(3)

    xplus_com = 0d0
    ! calculate tommorrow's annuitized asset stock
    if (ij_com == JR-1 .and. ann) then
      xplus_com = mx_com
    endif

    ! get tomorrow's year
    itp = year(it_com, ij_com, ij_com+1)

    ! get lsra transfer payment
    v_ind = v(0, ia_com, ix_com, ip_com, iw_com, ie_com, is_com, ij_com, it_com)

    ! calculate the wage rate and next periods pension claims
    wage = w(it_com)*eff(ij_com, is_com)*eta(iw_com, is_com)
    pplus_com = (p(ip_com)*dble(ij_com-1) + mu(it_com)*(lambda(it_com) &
           + (1d0-lambda(it_com))*min(wage*l_com/inc_bar(it_com), sscc(it_com))))/dble(ij_com)

    ! worker do not invest
    k_com = 0d0

    ! calculate contribution to pension system
    pencon_com = taup(it_com)*min(wage*l_com, sscc(it_com)*inc_bar(it_com))

    ! calculate income tax
    inctax_com = tarif(max(wage*l_com - 0.00d0*wage*l_com - 0.04d0*inc_bar(0) - pencon_com, 0d0) + pen(ip_com, ij_com, it_com))

    ! calculate capital gains tax
    captax_com = taur(it_com)*1.055d0*max(r(it_com)*a(ia_com)-2d0*0.0267d0*inc_bar(0), 0d0)

    ! calculate consumption
    c_com = ((1d0+r(it_com))*a(ia_com) + wage*l_com + beq(is_com, ij_com, it_com) + pen(ip_com, ij_com, it_com) + v_ind &
         - pencon_com - inctax_com - captax_com - mx_com - a_plus)*pinv(it_com)

    ! calculate tomorrow's part of the value function and occupational decision
    valuefunc_w = 0d0
    valuefunc_help = 0d0
    oplus_com = 0d0

    if (ij_com < JJ) then

      ! interpolate next period's value function as a worker/retiree
      valuefunc_w = interpolate_EV(0, a_plus, xplus_com, pplus_com, iw_com, ie_com, is_com, ij_com+1, itp)

      ! interpolate next period's value function as an entrepreneur
      if (ij_com < JR-1) then

        valuefunc_help = interpolate_EV(1, a_plus, xplus_com, pplus_com, iw_com, ie_com, is_com, ij_com+1, itp) - suc

        ! set next period's occupational decision
        if (valuefunc_help > valuefunc_w .and. ent) then
          valuefunc_w = valuefunc_help
          oplus_com = 1d0
        endif

      endif

    endif

    ! add today's part and discount
    valuefunc_w = -(util(c_com, l_com) + beta*psi(is_com, ij_com+1)*valuefunc_w + (1d0-psi(is_com, ij_com+1))*phi1*(1d0+a_plus*phi2)**(1d0-gamma))

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
    real*8 :: a_plus, p_hat, profit, v_ind, valuefunc_help
    real*8 :: temp1, temp2
    integer :: ij, itp, iij, itj

    ! tomorrow's assets
    a_plus = xy(1)

    ! today's investment
    k_com = xy(2)

    ! today's investment in annuitized assets
    mx_com = 0d0
    if (ij_com == JR-1 .and. ann) mx_com = xy(3)

    ! calculate annuities
    p_hat = 0d0
    if (ij_com >= JR) then
      temp1 = 0d0
      do ij = ij_com, JJ
        temp2 = 1d0
        do iij = ij_com+1, ij
          itj = year(it_com, ij_com, iij)
          temp2 = temp2*(1d0+r(itj))*psix(is_com, iij, itj)
        enddo
        temp1 = temp1 + 1d0/temp2
      enddo
      p_hat = x(ix_com)*(1d0+r(it_com))*psix(is_com, ij_com, it_com)/temp1
    endif

    xplus_com = 0d0
    ! calculate tommorrow's annuitized asset stock
    if (ij_com == JR-1 .and. ann) then
      xplus_com = mx_com
    elseif (ij_com >= JR) then
      xplus_com = x(ix_com)*(1d0+r(it_com))*psix(is_com, ij_com, it_com) - p_hat
    endif

    ! no investment without any assets
    if (ia_com == 0) k_com = 0d0

    ! today's fixed labor
    l_com = l_bar

    ! get tomorrows year
    itp = year(it_com, ij_com, ij_com+1)

    ! get lsra transfer payment
    v_ind = v(1, ia_com, ix_com, ip_com, iw_com, ie_com, is_com, ij_com, it_com)

    ! entrepreneur's profit
    profit = theta(ie_com, is_com)*(k_com**alpha*(eff(ij_com, is_com)*l_bar)**(1d0-alpha))**nu - delta*k_com - r(it_com)*max(k_com-a(ia_com), 0d0)

    ! calculate contribution to pension system
    if (ij_com < JR) then
        pencon_com = phi(it_com)*taup(it_com)*min(profit, sscc(it_com)*inc_bar(it_com))
    else
        pencon_com = 0d0
    endif

    ! calculate income tax
    inctax_com = tarif(max(profit - 0.00d0*profit - 0.04d0*inc_bar(0) - pencon_com, 0d0) + pen(ip_com, ij_com, it_com))

    ! calcualte capital gains tax
    captax_com = taur(it_com)*1.055d0*max(r(it_com)*max(a(ia_com)-k_com, 0d0) - 2d0*0.0267d0*inc_bar(0), 0d0)

    ! calculate consumption
    c_com =  (a(ia_com) + r(it_com)*max(a(ia_com)-k_com, 0d0) + profit + beq(is_com, ij_com, it_com) + pen(ip_com, ij_com, it_com) + p_hat + v_ind  &
           - captax_com - inctax_com - pencon_com - mx_com - a_plus)*pinv(it_com)

    ! calculate next periods pension claims
    if (ij_com < JR) then
      pplus_com = (p(ip_com)*dble(ij_com-1) + mu(it_com)*phi(it_com)*(lambda(it_com) &
             + (1d0-lambda(it_com))*min(profit/inc_bar(it_com), sscc(it_com))))/dble(ij_com)
    else
      pplus_com = p(ip_com)
    endif

    ! calculate tomorrow's part of the value function and occupational decision
    valuefunc_e = 0d0
    valuefunc_help = 0d0
    oplus_com = 0d0

    if (ij_com < JJ) then

      ! interpolate next period's value function as a worker/retiree
      valuefunc_e = interpolate_EV(0, a_plus, xplus_com, pplus_com, iw_com, ie_com, is_com, ij_com+1, itp) - suc

      ! interpolate next period's value function as an entrepreneur
      if (ij_com < JE-1) then

        valuefunc_help = interpolate_EV(1, a_plus, xplus_com, pplus_com, iw_com, ie_com, is_com, ij_com+1, itp)

        ! set next period's occupational decision
        if (valuefunc_help > valuefunc_e .and. ent) then
          valuefunc_e = valuefunc_help
          oplus_com = 1d0
        endif

      endif

    endif

    ! add today's part and discount
    valuefunc_e = -(util(c_com, l_com) + beta*psi(is_com, ij_com+1)*valuefunc_e + (1d0-psi(is_com, ij_com+1))*phi1*(1d0+a_plus*phi2)**(1d0-gamma))

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
    real*8 :: a_plus, p_hat, v_ind
    real*8 :: temp1, temp2
    integer :: ij, iij, itj, itp

    ! tomorrow's assets
    a_plus = xy

    ! today's labor
    l_com  = 0d0

    ! today's investment in annuitized assets
    mx_com = 0d0

    ! get tomorrow's year
    itp = year(it_com, ij_com, ij_com+1)

    ! get lsra transfer payment
    v_ind = v(0, ia_com, ix_com, ip_com, iw_com, ie_com, is_com, ij_com, it_com)

    ! calculate the wage rate and next periods pension claims
    pplus_com = p(ip_com)

    ! retirees do not invest
    k_com = 0d0

    ! calculate annuities
    p_hat = 0d0
    temp1 = 0d0
    do ij = ij_com, JJ
      temp2 = 1d0
      do iij = ij_com+1, ij
        itj = year(it_com, ij_com, iij)
        temp2 = temp2*(1d0+r(itj))*psix(is_com, iij, itj)
      enddo
      temp1 = temp1 + 1d0/temp2
    enddo
    p_hat = x(ix_com)*(1d0+r(it_com))*psix(is_com, ij_com, it_com)/temp1

    ! calculate tomorrow's annuitized asset stock
    xplus_com = x(ix_com)*(1d0+r(it_com))*psix(is_com, ij_com, it_com) - p_hat

    ! calculate contribution to pension system
    pencon_com = 0d0

    ! calculate income tax
    inctax_com = tarif(pen(ip_com, ij_com, it_com))

    ! calculate capital gains tax
    captax_com = taur(it_com)*1.055d0*max(r(it_com)*a(ia_com)-2d0*0.0267d0*inc_bar(0), 0d0)

    ! calculate consumption
    c_com = ((1d0+r(it_com))*a(ia_com) + beq(is_com, ij_com, it_com) + pen(ip_com, ij_com, it_com) + p_hat + v_ind &
         - inctax_com - captax_com - a_plus)*pinv(it_com)

    ! calculate tomorrow's part of the value function and occupational decision
    valuefunc_r = 0d0
    oplus_com = 0d0

    if (ij_com < JJ) then

      ! interpolate next period's value function as a retiree
      valuefunc_r = interpolate_EV(0, a_plus, xplus_com, pplus_com, iw_com, ie_com, is_com, ij_com+1, itp)

    endif

    ! add today's part and discount
    valuefunc_r = -(util(c_com, l_com) + beta*psi(is_com, ij_com+1)*valuefunc_r + (1d0-psi(is_com, ij_com+1))*phi1*(1d0+a_plus*phi2)**(1d0-gamma))

  end function



  !##############################################################################
  ! FUNCTION profent
  !
  ! Computes profit of an entrepreneur
  !##############################################################################
  function profent(k, ij, ia, is, ie, it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: k
    integer, intent(in) :: ij, ia, is, ie, it
    real*8 :: profent

    ! compute profit
    profent = theta(ie, is)*(k**alpha*(eff(ij, is)*l_bar)**(1d0-alpha))**nu - delta*k - r(it)*max(k-a(ia), 0d0)

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
  function margu(cons, lab, it)

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: cons, lab
    integer, intent(in) :: it
    real*8 :: margu

    ! determine marginal utility
    margu = sigma*(cons**sigma*(1d0-lab)**(1d0-sigma))**(1d0-gamma)/((1d0+tauc(it))*cons)

  end function


  !##############################################################################
  ! FUNCTION interpolate_EV
  !
  ! Interpolates the expected valuefunction
  !##############################################################################
  function interpolate_EV(io, a_plus, x_plus, p_plus, iw, ie, is, ij, it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: a_plus, x_plus, p_plus
    integer, intent(in) :: io, iw, ie, is, ij, it
    real*8 :: interpolate_EV

    !##### OTHER VARIABLES ####################################################
    real*8 :: varphi, varchi, varpsi
    integer :: ial, iar, ixl, ixr, ipl, ipr

    ! interpolate value function
    call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)
    call linint_Grow(x_plus, x_l, x_u, x_grow, NX, ixl, ixr, varchi)
    call linint_Equi(p_plus, p_l, p_u, NP, ipl, ipr, varpsi)

    interpolate_EV = (varphi*varchi*varpsi*EV(io, ial, ixl, ipl, iw, ie, is, ij, it) &
                      + varphi*varchi*(1d0-varpsi)*EV(io, ial, ixl, ipr, iw, ie, is, ij, it) &
                      + varphi*(1d0-varchi)*varpsi*EV(io, ial, ixr, ipl, iw, ie, is, ij, it) &
                      + varphi*(1d0-varchi)*(1d0-varpsi)*EV(io, ial, ixr, ipr, iw, ie, is, ij, it) &
                      + (1d0-varphi)*varchi*varpsi*EV(io, iar, ixl, ipl, iw, ie, is, ij, it) &
                      + (1d0-varphi)*varchi*(1d0-varpsi)*EV(io, iar, ixl, ipr, iw, ie, is, ij, it) &
                      + (1d0-varphi)*(1d0-varchi)*varpsi*EV(io, iar, ixr, ipl, iw, ie, is, ij, it) &
                      + (1d0-varphi)*(1d0-varchi)*(1d0-varpsi)*EV(io, iar, ixr, ipr, iw, ie, is, ij, it)) &
                     **(1d0-gamma)/(1d0-gamma)

    end function


  !##############################################################################
  ! FUNCTION year
  !
  ! Calculates year at which age ij agent is ijj
  !##############################################################################
  function year(it, ij, ijj)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it, ij, ijj
    integer :: year

    ! calculate year
    year = it + ijj - ij

    if(it == 0 .or. year <= 0)year = 0
    if(it == TT .or. year >= TT)year = TT

  end function


  !##############################################################################
  ! FUNCTION check_grid
  !
  ! Checks for the maximum gridpoint used
  !##############################################################################
  subroutine check_grid(iamax, ixmax, it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it
    integer :: iamax(JJ), ixmax(JJ)

    !##### OTHER VARIABLES ####################################################
    integer :: ia, ix, ip, iw, ie, is, ij

    iamax = 0
    ixmax = 0

    do ij = 1, JJ

      ! check for the maximum asset grid point used at a certain age
      do ia = NA, 0, -1
        if (sum(m(:, ia, :, :, :, :, :, ij, it)) > 0d0) then
          iamax(ij) = ia
          exit
        endif
      enddo ! ia

      ! check for the maximum annuitie grid point used at a certain age
      if (ann) then
        do ix = NX, 0, -1
          if (sum(m(:, :, ix, :, :, :, :, ij, it)) > 0d0) then
            ixmax(ij) = ix
            exit
          endif
        enddo ! ix
      endif

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
  ! SUBROUTINE gini_wealth
  !
  ! Calculates GINI coefficient of wealth
  !##############################################################################
  subroutine gini_wealth(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES ####################################################
    integer :: is, ie, iw, ip, ix, ia, ij, ic, ICMAX
    real*8, allocatable :: xs(:), ys(:), xcum(:), ycum(:)

    if(allocated(xs))deallocate(xs)
    if(allocated(ys))deallocate(ys)
    if(allocated(xcum))deallocate(xcum)
    if(allocated(ycum))deallocate(ycum)

    allocate(xs(2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ))
    allocate(ys(2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ))
    allocate(xcum(2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ))
    allocate(ycum(2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ))

    ic = 1

    do ij = 1, JJ
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ia = 0, NA

                  if(m(0, ia, ip, ix, iw, ie, is, ij, it) > 0d0) then
                    ys(ic) = m(0, ia, ip, ix, iw, ie, is, ij, it)
                    xs(ic) = a(ia)
                    ic = ic + 1
                  endif
                  if(m(1, ia, ip, ix, iw, ie, is, ij, it) > 0d0) then
                    ys(ic) = m(1, ia, ip, ix, iw, ie, is, ij, it)
                    xs(ic) = a(ia)
                    ic = ic + 1
                  endif

                enddo ! ia
              enddo ! ix
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

    gini_w(it) = 0d0
    do ic = 2, ICMAX+1
      gini_w(it) = gini_w(it) + ys(ic-1)*(xcum(ic-1)+xcum(ic))
    enddo
    gini_w(it) = 1d0-gini_w(it)/xcum(ICMAX+1)

    percentiles_w = 0d0
    do ic = ICMAX, 1, -1
      if (1d0-ycum(ic) > 0.01 .and. percentiles_w(1, it) <= 0d0) then
        percentiles_w(1, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.05 .and. percentiles_w(2, it) <= 0d0) then
        percentiles_w(2, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.10 .and. percentiles_w(3, it) <= 0d0) then
        percentiles_w(3, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.20 .and. percentiles_w(4, it) <= 0d0) then
        percentiles_w(4, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.40 .and. percentiles_w(5, it) <= 0d0) then
        percentiles_w(5, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.60 .and. percentiles_w(6, it) <= 0d0) then
        percentiles_w(6, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
    enddo

  end subroutine


  !##############################################################################
  ! SUBROUTINE gini_income
  !
  ! Calculates GINI coefficient of income
  !##############################################################################
  subroutine gini_income(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES ####################################################
    integer :: is, ie, iw, ip, ix, ia, ij, ic, ICMAX
    real*8, allocatable :: xs(:), ys(:), xcum(:), ycum(:)

    if(allocated(xs))deallocate(xs)
    if(allocated(ys))deallocate(ys)
    if(allocated(xcum))deallocate(xcum)
    if(allocated(ycum))deallocate(ycum)

    allocate(xs(2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ))
    allocate(ys(2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ))
    allocate(xcum(2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ))
    allocate(ycum(2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ))

    ic = 1

    do ij = 1, JJ
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ia = 0, NA

                  if(m(0, ia, ip, ix, iw, ie, is, ij, it) > 0d0) then
                    ys(ic) = m(0, ia, ip, ix, iw, ie, is, ij, it)
                    xs(ic) = a(ia)*r(it) + pen(ip, ij, it) + eff(ij, is)*eta(iw, is)*l(0, ia, ip, ix, iw, ie, is, ij, it)*w(it)
                    ic = ic + 1
                  endif
                  if(m(1, ia, ip, ix, iw, ie, is, ij, it) > 0d0) then
                    ys(ic) = m(1, ia, ip, ix, iw, ie, is, ij, it)
                    xs(ic) = max(a(ia)-k(1, ia, ip, ix, iw, ie, is, ij, it), 0d0)*r(it) + pen(ip, ij, it) + profent(k(1, ia, ip, ix, iw, ie, is, ij, it), ij, ia, is, ie, it)
                    ic = ic + 1
                  endif

                enddo ! ia
              enddo ! ix
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

    gini_i(it) = 0d0
    do ic = 2, ICMAX+1
      gini_i(it) = gini_i(it) + ys(ic-1)*(xcum(ic-1)+xcum(ic))
    enddo
    gini_i(it) = 1d0-gini_i(it)/xcum(ICMAX+1)

    percentiles_i = 0d0
    do ic = ICMAX, 1, -1
      if (1d0-ycum(ic) > 0.01d0 .and. percentiles_i(1, it) <= 0d0) then
        percentiles_i(1, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.05d0 .and. percentiles_i(2, it) <= 0d0) then
        percentiles_i(2, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.10d0 .and. percentiles_i(3, it) <= 0d0) then
        percentiles_i(3, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.20d0 .and. percentiles_i(4, it) <= 0d0) then
        percentiles_i(4, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.40d0 .and. percentiles_i(5, it) <= 0d0) then
        percentiles_i(5, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
      if (1d0-ycum(ic) > 0.60d0 .and. percentiles_i(6, it) <= 0d0) then
        percentiles_i(6, it) = (xcum(ICMAX+1)-xcum(ic-1))/xcum(ICMAX+1)
      endif
    enddo

  end subroutine


  !##############################################################################
  ! FUNCTION gini
  !
  ! Calculates GINI coefficient of x with mass y
  !##############################################################################
  function gini(x, y)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: x(:), y(:)

    !##### OTHER VARIABLES ####################################################
    integer :: ii, ic, ICMAX
    real*8, allocatable :: xs(:), ys(:), xcum(:), ycum(:)
    real*8 :: gini

    if(allocated(xs))deallocate(xs)
    if(allocated(ys))deallocate(ys)
    if(allocated(xcum))deallocate(xcum)
    if(allocated(ycum))deallocate(ycum)

    allocate(xs(size(x)))
    allocate(ys(size(y)))
    allocate(xcum(size(x)))
    allocate(ycum(size(y)))

    ic = 1

    do ii = 1, size(x)
      if(y(ii) > 0d0) then
        ys(ic) = y(ii)
        xs(ic) = x(ii)
        ic = ic + 1
      endif
    enddo ! ii

    ! get array size and normalize ys
    ICMAX = ic - 1
    ys(1:ICMAX) = ys(1:ICMAX)/sum(ys(1:ICMAX))

    ! sort array
    call quick_sort(xs(1:ICMAX), ys(1:ICMAX))

    ! calculate cumulative distributions
    xcum = 0d0
    ycum = 0d0
    xcum(1) = xs(1)*ys(1)
    ycum(1) = ys(1)
    do ic = 2, ICMAX
      xcum(ic) = xcum(ic-1) + xs(ic)*ys(ic)
      ycum(ic) = ycum(ic-1) + ys(ic)
    enddo

    gini = ys(1)*xcum(1)
    do ic = 2, ICMAX
      gini = gini + ys(ic-1)*(xcum(ic-1)+xcum(ic))
    enddo
    gini = 1d0 - gini/xcum(ICMAX)

  end function


  !##############################################################################
  ! FUNCTION gini
  !
  ! Calculates percentiles p of x with mass y
  !##############################################################################
  function percentiles(x, y, p)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: x(:), y(:), p(:)

    !##### OTHER VARIABLES ####################################################
    integer :: ii, ic, ICMAX
    real*8, allocatable :: xs(:), ys(:), xcum(:), ycum(:)
    real*8 :: percentiles(size(p))

    if(allocated(xs))deallocate(xs)
    if(allocated(ys))deallocate(ys)
    if(allocated(xcum))deallocate(xcum)
    if(allocated(ycum))deallocate(ycum)

    allocate(xs(size(x)))
    allocate(ys(size(y)))
    allocate(xcum(size(x)))
    allocate(ycum(size(y)))

    ic = 1

    do ii = 1, size(x)
      if(y(ii) > 0d0) then
        ys(ic) = y(ii)
        xs(ic) = x(ii)
        ic = ic + 1
      endif
    enddo ! ii

    ! get array size and normalize ys
    ICMAX = ic - 1
    ys(1:ICMAX) = ys(1:ICMAX)/sum(ys(1:ICMAX))

    ! sort array
    call quick_sort(xs(1:ICMAX), ys(1:ICMAX))

    ! calculate cumulative distributions
    xcum = 0d0
    ycum = 0d0
    xcum(1) = xs(1)*ys(1)
    ycum(1) = ys(1)
    do ic = 2, ICMAX
      xcum(ic) = xcum(ic-1) + xs(ic)*ys(ic)
      ycum(ic) = ycum(ic-1) + ys(ic)
    enddo

    percentiles = 0d0
    do ic = ICMAX, 1, -1
      do ii = 1, size(p)

        if (1d0-ycum(ic) >= p(ii) .and. percentiles(ii) <= 0d0) then
          percentiles(ii) = (xcum(ICMAX)-xcum(ic-1))/xcum(ICMAX)*100d0
          !percentiles(ii) = xs(ic)
        endif

      end do
    enddo

  end function


end module
