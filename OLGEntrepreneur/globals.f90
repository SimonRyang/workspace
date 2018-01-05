module globals

  use toolbox

  implicit none

  ! number of years the household lives
  integer, parameter :: JJ = 16

  ! number of years the household retires
  integer, parameter :: JR = 10

  ! number of permanent skill classes
  integer, parameter :: NS = 3

  ! number of productivity (eta) shocks
  integer, parameter :: NW = 5

  ! number of entrepreneurial ability (theta) shocks
  integer, parameter :: NE = 2

  ! number of points on the total asset grid
  integer, parameter :: NQ = 16

  ! number of points on the liquid asset grid
  integer, parameter :: NA = 16

  ! number of points on the capital grid
  integer, parameter :: NK = 0

  ! number of points on the annuity asset grid
  integer, parameter :: NX = 0

  ! number of points on the pension claim grid
  integer, parameter :: NP = 4

  ! demographic parameters
  real*8, parameter :: n_p = (1d0+0.007d0)**5-1d0 !
  real*8, parameter :: dist_skill(NS) = (/0.1520d0, 0.5547d0, 0.2933d0/) !

  ! macroeconomic parameters
  real*8, parameter :: gy = 0.19d0 !
  real*8, parameter :: by = 0.60d0/5d0 !

  ! government parameters
  real*8, parameter :: tauy = 0.150d0 !
  real*8, parameter :: taur = 0.250d0 !

  ! household preference parameters
  real*8, parameter :: gamma = 0.5d0 !
  real*8, parameter :: egam  = 1d0 - 1d0/gamma !
  real*8, parameter :: sigma = 0.320d0
  real*8, parameter :: beta  = 0.98d0**5
  real*8, parameter :: mu_b  = 0.10d0

  ! maximum investment in annuities
  real*8, parameter :: mx_max = 0.10d0

  ! capital parameters
  real*8, parameter :: delta_k = 1d0 - (1d0-0.06d0)**5d0 !
  real*8, parameter :: xi      = 1d0/3d0

  ! production parameters
  real*8, parameter :: Omega = 1.0d0
  real*8, parameter :: k_min = 0.2d0
  real*8, parameter :: phi_k = 0.25d0
  real*8, parameter :: alpha = 0.36d0 !
  real*8, parameter :: nu    = 0.88d0

  ! size of the total asset grid
  real*8, parameter :: Q_l    = 0d0
  real*8, parameter :: Q_u    = 12d0
  real*8, parameter :: Q_grow = 0.05d0

  ! size of the liquid asset grid
  real*8, parameter :: a_l    = Q_l
  real*8, parameter :: a_u    = Q_u
  real*8, parameter :: a_grow = Q_grow

  ! size of the capital grid
  real*8, parameter :: k_l    = k_min
  real*8, parameter :: k_u    = 0.25d0*Q_u/(1d0-xi)
  real*8, parameter :: k_grow = Q_grow

  ! size of the annuity asset grid
  real*8, parameter :: x_l    = Q_l
  real*8, parameter :: x_u    = 4d0*Q_u
  real*8, parameter :: x_grow = Q_grow

  ! size of the pension claim grid
  real*8, parameter :: p_l    = 0d0
  real*8, parameter :: p_u    = 2d0

  ! pension fraction of average income
  real*8, parameter :: kappa = 0.55d0

  ! numerical parameters
  integer, parameter :: itermax = 200
  real*8, parameter :: sig  = 1d-6
  real*8, parameter :: damp = 0.45d0

  ! measure time
  integer :: time

  ! discretized shocks
  real*8 :: dist_eta(NW, NS), pi_eta(NW, NW, NS), eta(NW, NS), dist_theta(NE, NS), pi_theta(NE, NE, NS), theta(NE, NS)

  ! demographic and other model parameters
  real*8 :: eff(NS, JJ)
  real*8 :: pen(0:NP, JJ), ann(0:NX, NS, JJ), ans(0:NX, NS, JJ)
  real*8 :: psi(NS, JJ+1), rpop(NS, JJ)
  real*8 :: bqs(NS), beq(NS, JJ), Gama(JJ)

  ! government variables
  real*8 :: mu, phi, lambda
  real*8 :: tauc, taup

  ! progressive income tax
  real*8, parameter :: t1 = 0.14d0, t2 = 0.24d0, t3 = 0.45d0
  real*8 :: r1, r2, r3, b1, b2

  ! macroeconomic variables
  real*8 :: r, w
  real*8 :: ybar, pinv
  real*8 :: AA, AX, BQ, PBEN, PCON
  real*8 :: KK, KC, KE, LC, BB
  real*8 :: YY, YC, YE, CC, II, TC, GG
  real*8 :: TAc, TAr, TAw, TAy

  ! different grids to discretize the state space
  real*8 :: Q(0:NQ), a(0:NA), k(0:NK), x(0:NX), p(0:NP)

  ! variables to store the policy functions
  real*8 :: Q_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: a_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), k_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: x_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), p_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: inctax(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), captax(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: penben(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), pencon(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: c(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), l(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)

  ! variables to store the value function
  real*8 :: V(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), EV(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)

  ! variables for temporary policy and value functions
  real*8 :: Q_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: a_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), k_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: x_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), p_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: inctax_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), captax_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: penben_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), pencon_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: c_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), l_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: V_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)

  ! variables to store the savings choice decisions
  real*8 :: omega_x_t(0:1, 0:NQ, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), omega_k_t(0:1, 0:NQ, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)
  real*8 :: S(0:1, 0:NQ, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)

  ! weights for the different gridpoints on the discretized state space
  real*8 :: m(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ), m_Q(0:NQ, 0:NK, 0:NX, 0:NP, NW, NE, NS, JJ)

  ! numerical variables
  integer :: iter
  integer :: ij_com, iq_com, ia_com, ix_com, ip_com, ik_com, iw_com, ie_com, is_com
  integer :: iq_p_com, ia_p_com, ip_p_com, io_p_com
  integer :: iqmax(JJ), iamax(JJ), ixmax(JJ), ikmax(JJ)
  real*8 :: cons_com, lab_com, x_plus_com, p_plus_com
  real*8 :: inctax_com, captax_com, pencon_com, aas_com
  real*8 :: DIFF

  !$omp threadprivate(ij_com, iq_com, ia_com, ix_com, ip_com, ik_com, iw_com, ie_com, is_com)
  !$omp threadprivate(iq_p_com, ia_p_com, ip_p_com, io_p_com)
  !$omp threadprivate(cons_com, lab_com, x_plus_com, p_plus_com)
  !$omp threadprivate(inctax_com, captax_com, pencon_com, aas_com)

contains


  !#############################################################################
  ! FUNCTION solve_worker
  !
  ! solves the future worker's decision of
  ! how much wealth to invest into annuities
  !#############################################################################
  subroutine solve_worker(iq_p, ik, ix, ip_p, iw, ie, is, ij)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: iq_p, ik, ix, ip_p, iw, ie, is, ij

    !##### OTHER VARIABLES #####################################################
    real*8 :: x_in, fret

    ! set up communication variables
    iq_p_com = iq_p; ik_com = ik; ix_com = ix; ip_p_com = ip_p
    iw_com = iw; ie_com = ie; is_com = is; ij_com = ij

    if (Q(iq_p) > 0d0) then

      ! get best initial guess from future period
      x_in = max(omega_x_t(0, iq_p, ik, ix, ip_p, iw, ie, is, ij), 1d-4)

      ! solve the household problem using fminsearch
      call fminsearch(x_in, fret, 0d0, 1d0, inv_w)

      ! wealth share for annuities and firm capital
      omega_x_t(0, iq_p, ik, ix, ip_p, iw, ie, is, ij) = x_in
      omega_k_t(0, iq_p, ik, ix, ip_p, iw, ie, is, ij) = 0d0
      S(0, iq_p, ik, ix, ip_p, iw, ie, is, ij) = -fret

    else

      ! wealth share for annuities and firm capital
      omega_x_t(0, iq_p, ik, ix, ip_p, iw, ie, is, ij) = 0d0
      omega_k_t(0, iq_p, ik, ix, ip_p, iw, ie, is, ij) = 0d0
      S(0, iq_p, ik, ix, ip_p, iw, ie, is, ij) = -inv_w(0d0)

    endif

  end subroutine


  !#############################################################################
  ! FUNCTION solve_entrepreneur
  !
  ! solve the future entrepreneur's decision
  ! of how much wealth to invest into annuities and into firm capital
  !#############################################################################
  subroutine solve_entrepreneur(iq_p, ik, ix, ip_p, iw, ie, is, ij)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: iq_p, ik, ix, ip_p, iw, ie, is, ij

    !##### OTHER VARIABLES #####################################################
    real*8 :: x_in(2), fret

    ! set up communication variables
    iq_p_com = iq_p; ik_com = ik; ix_com = ix; ip_p_com = ip_p
    iw_com = iw; ie_com = ie; is_com = is; ij_com = ij

    if (Q(iq_p) > (1d0-xi)*k_min + tr(k(ik), k_min)) then

      ! get best initial guess from future period
      x_in(1) = max(omega_x_t(1, iq_p, ik, ix, ip_p, iw, ie, is, ij), 1d-4)
      x_in(2) = max(omega_k_t(1, iq_p, ik, ix, ip_p, iw, ie, is, ij), 1d-4)

      ! solve the household problem using fminsearch
      call fminsearch(x_in, fret, (/0d0, 0d0/), (/1d0, 1d0/), inv_e)

      ! wealth share for annuities and firm capital
      omega_x_t(1, iq_p, ik, ix, ip_p, iw, ie, is, ij) = x_in(1)
      omega_k_t(1, iq_p, ik, ix, ip_p, iw, ie, is, ij) = x_in(2)
      S(1, iq_p, ik, ix, ip_p, iw, ie, is, ij) = -fret

    else

      ! wealth share for annuities and firm capital
      omega_x_t(1, iq_p, ik, ix, ip_p, iw, ie, is, ij) = 1d0
      omega_k_t(1, iq_p, ik, ix, ip_p, iw, ie, is, ij) = 1d0
      S(1, iq_p, ik, ix, ip_p, iw, ie, is, ij) = 1d-13**egam/egam !-inv_e((/1d0, 1d0/))

    endif

  end subroutine


  !#############################################################################
  ! FUNCTION solve_retiree
  !
  ! solve the future retiree's decision
  !#############################################################################
  subroutine solve_retiree(iq_p, ik, ix, ip_p, iw, ie, is, ij)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: iq_p, ik, ix, ip_p, iw, ie, is, ij

    !##### OTHER VARIABLES #####################################################
    integer :: ial, iar
    real*8 :: a_p, EV_temp, S_temp, varphi_a

    a_p = Q(iq_p)

    ! derive interpolation weights
    call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

    ! restrict values to grid just in case
    ial = min(ial, NA)
    iar = min(iar, NA)
    varphi_a = max(min(varphi_a, 1d0),0d0)

    ! calculate future part of the value function
    EV_temp = (varphi_a      *(egam*EV(ial, 0, ix, ip_p, iw, ie, is, ij+1))**(1d0/egam) + &
               (1d0-varphi_a)*(egam*EV(iar, 0, ix, ip_p, iw, ie, is, ij+1))**(1d0/egam))**egam/egam

    ! calculate bequest part of the value function
    S_temp = (1d0-psi(is, ij+1))*mu_b*max(a_p, 1d-13)**egam/egam

    ! wealth share for annuities and firm capital
    omega_x_t(:, iq_p, ik, ix, ip_p, iw, ie, is, ij) = 0d0
    omega_k_t(:, iq_p, ik, ix, ip_p, iw, ie, is, ij) = 0d0
    S(:, iq_p, ik, ix, ip_p, iw, ie, is, ij) = psi(is, ij+1)*beta*EV_temp + S_temp

  end subroutine


  !#############################################################################
  ! FUNCTION solve_consumption
  !
  ! solves the household's consumption-savings-working decision
  !#############################################################################
  subroutine solve_consumption(io_p, ia, ik, ix, ip, iw, ie, is, ij)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: io_p, ia, ik, ix, ip, iw, ie, is, ij

    !##### OTHER VARIABLES #####################################################
    real*8 :: x_in(2), fret, x_p, mx, k_p, varphi_q, varphi_p
    integer :: iql, iqr, ipl, ipr

    ! set up communication variables
    io_p_com = io_p; ia_com = ia; ik_com = ik; ix_com = ix; ip_com = ip
    iw_com = iw; ie_com = ie; is_com = is; ij_com = ij

    ! get best initial guess from future period
    x_in(1) = max(Q_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij), 1d-4)
    x_in(2) = max(l_t(io_p, ia, ik, ix, ip, iw, ie, is, ij), 1d-4)

    ! solve the household problem using fminsearch
    if (ij < JR) then
      call fminsearch(x_in, fret, (/Q_l, 0d0/), (/Q_u, 0.99d0/), cons_o)
    else
      call fminsearch(x_in(1), fret, Q_l, Q_u, cons_r)
    endif

    ! derive interpolation weights
    call linint_Grow(x_in(1), Q_l, Q_u, Q_grow, NQ, iql, iqr, varphi_q)
    call linint_Equi(p_plus_com, p_l, p_u, NP, ipl, ipr, varphi_p)

    ! restrict values to grid just in case
    iql = min(iql, NQ)
    iqr = min(iqr, NQ)
    varphi_q = max(min(varphi_q, 1d0),0d0)

    ! restrict values to grid just in case
    ipl = min(ipl, NP)
    ipr = min(ipr, NP)
    varphi_p = max(min(varphi_p, 1d0),0d0)

    ! determine future investment
    k_p = 0d0

    if (io_p == 1) then
     if (varphi_q <= varphi_p) then
       k_p = ((1d0-xi)*k_min + (varphi_q           *omega_k_t(io_p, iql, ik, ix, ipl, iw, ie, is, ij) +  &
                                (varphi_p-varphi_q)*omega_k_t(io_p, iqr, ik, ix, ipl, iw, ie, is, ij) +  &
                                (1d0-varphi_p)     *omega_k_t(io_p, iqr, ik, ix, ipr, iw, ie, is, ij))*(x_in(1)-(1d0-xi)*k_min))/(1d0-xi)
     else
       k_p = ((1d0-xi)*k_min + (varphi_p           *omega_k_t(io_p, iql, ik, ix, ipl, iw, ie, is, ij) +  &
                                (varphi_q-varphi_p)*omega_k_t(io_p, iql, ik, ix, ipr, iw, ie, is, ij) +  &
                                (1d0-varphi_q)     *omega_k_t(io_p, iqr, ik, ix, ipr, iw, ie, is, ij))*(x_in(1)-(1d0-xi)*k_min))/(1d0-xi)
     endif
    endif

    ! determine future annuity stock
    x_p = 0d0
    mx = 0d0

    if (ij < JR) then
      if (varphi_q <= varphi_p) then
        mx = min((varphi_q           *omega_x_t(io_p, iql, ik, ix, ipl, iw, ie, is, ij) +  &
                  (varphi_p-varphi_q)*omega_x_t(io_p, iqr, ik, ix, ipl, iw, ie, is, ij) +  &
                  (1d0-varphi_p)     *omega_x_t(io_p, iqr, ik, ix, ipr, iw, ie, is, ij))*x_in(1), x_in(1) - (1d0-xi)*k_p - tr(k(ik), k_p))
      else
        mx = min((varphi_p           *omega_x_t(io_p, iql, ik, ix, ipl, iw, ie, is, ij) +  &
                  (varphi_q-varphi_p)*omega_x_t(io_p, iql, ik, ix, ipr, iw, ie, is, ij) +  &
                  (1d0-varphi_q)     *omega_x_t(io_p, iqr, ik, ix, ipr, iw, ie, is, ij))*x_in(1), x_in(1) - (1d0-xi)*k_p - tr(k(ik), k_p))
      endif
      x_p = (1d0+r)/psi(is, ij)*x(ix)+ mx
    else
      x_p = x(ix)
    endif

    ! copy decisions
    Q_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = x_in(1)
    a_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = x_in(1) - (1d0-xi)*k_p - mx - tr(k(ik), k_p)
    k_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = k_p
    x_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = x_p
    p_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = p_plus_com
    inctax_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = inctax_com
    captax_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = captax_com
    penben_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = pen(ip, ij)
    pencon_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = pencon_com
    c_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) =  (aas_com - x_in(1))*pinv
    l_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = lab_com
    V_t(io_p, ia, ik, ix, ip, iw, ie, is, ij) = -fret

  end subroutine


  !#############################################################################
  ! FUNCTION inv_w
  !
  ! calculates the future worker's future part of the value function
  !#############################################################################
  function inv_w(x_in)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    real*8, intent(in) :: x_in
    real*8 :: inv_w

    !##### OTHER VARIABLES #####################################################
    real*8 :: a_p, a_temp, x_p, EV_temp, S_temp, omega_x, varphi_a, varphi_x
    integer :: ial, iar, ixl, ixr

    ! store annuity share
    omega_x  = x_in

    ! determine future liquid wealth and future annuity asset stock
    x_p = (1d0+r)/psi(is_com, ij_com)*x(ix_com) + min(omega_x*Q(iq_p_com), mx_max)
    a_temp = (1d0-omega_x)*Q(iq_p_com)
    a_p = max(a_temp, 0d0)

    ! derive interpolation weights
    call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
    if (NX > 0) then
      call linint_Grow(x_p, x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)
    else
      ixl = 0; ixr = 0; varphi_x = 1d0
    endif

    ! restrict values to grid just in case
    ial = min(ial, NA)
    iar = min(iar, NA)
    varphi_a = max(min(varphi_a, 1d0),0d0)

    ! restrict values to grid just in case
    ixl = min(ixl, NX)
    ixr = min(ixr, NX)
    varphi_x = max(min(varphi_x, 1d0),0d0)

    ! calculate future part of the value function
    if (varphi_a <= varphi_x) then
      EV_temp = (varphi_a           *(egam*EV(ial, 0, ixl, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
                 (varphi_x-varphi_a)*(egam*EV(iar, 0, ixl, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
                 (1d0-varphi_x)     *(egam*EV(iar, 0, ixr, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam))**egam/egam
    else
      EV_temp = (varphi_x           *(egam*EV(ial, 0, ixl, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
                 (varphi_a-varphi_x)*(egam*EV(ial, 0, ixr, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
                 (1d0-varphi_a)     *(egam*EV(iar, 0, ixr, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam))**egam/egam
    endif

    ! calculate bequest part of the value function
    S_temp = 0d0 !(1d0-psi(is_com, ij_com+1))*mu_b*max(a_p, 1d-13)**egam/egam

    ! calculate future part of the value function
    if (a_temp < 0d0) then
      inv_w = -1d-13**egam/egam*(1d0+abs(a_temp))
    else
      inv_w = - (psi(is_com, ij_com+1)*beta*EV_temp + S_temp)
    endif

  end function


  !#############################################################################
  ! FUNCTION inv_e
  !
  ! calculates the future entrepreneur's future part of the value function
  !#############################################################################
  function inv_e(x_in)

  implicit none

  !##### INPUT/OUTPUT VARIABLES ##############################################
  real*8, intent(in) :: x_in(:)
  real*8 :: inv_e

  !##### OTHER VARIABLES ######################################################
  real*8 :: a_p, x_p, k_p, EV_temp, S_temp, omega_x, omega_k, varphi_a, varphi_k, varphi_x, a_temp
  integer :: ial, iar, ikl, ikr, ixl, ixr

    ! store annuity and firm capital share
    omega_x  = x_in(1)
    omega_k  = x_in(2)

    ! determine future liquid wealth, future firm capital and future annuity asset stock
    x_p = (1d0+r)/psi(is_com, ij_com)*x(ix_com) + min(omega_x*Q(iq_p_com), mx_max)
    k_p = ((1d0-xi)*k_min + omega_k*(Q(iq_p_com) - (1d0-xi)*k_min))/(1d0-xi)
    a_temp = Q(iq_p_com) - omega_x*Q(iq_p_com) - (1d0-xi)*k_p - tr(k(ik_com), k_p)
    a_p = max(a_temp, 0d0)

    ! derive interpolation weights
    call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
    if (NX > 0) then
      call linint_Grow(x_p, x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)
    else
      ixl = 0; ixr = 0; varphi_x = 1d0
    endif
    if (NK > 0) then
      call linint_Grow(k_p, k_l, k_u, k_grow, NK-1, ikl, ikr, varphi_k)
      ikl = min(ikl+1, NK)
      ikr = min(ikr+1, NK)
      varphi_k = max(min(varphi_k, 1d0), 0d0)
    else
      ikl = 0; ikr = 0; varphi_k = 1d0
    endif

    ! restrict values to grid just in case
    ial = min(ial, NA)
    iar = min(iar, NA)
    varphi_a = max(min(varphi_a, 1d0),0d0)

    ! restrict values to grid just in case
    ixl = min(ixl, NX)
    ixr = min(ixr, NX)
    varphi_x = max(min(varphi_x, 1d0),0d0)

    ! calculate future part of the value function
    EV_temp = (varphi_a*varphi_x*varphi_k                  *(egam*EV(ial, ikl, ixl, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
               varphi_a*varphi_k*(1d0-varphi_x)            *(egam*EV(ial, ikl, ixr, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
               varphi_a*(1d0-varphi_k)*varphi_x            *(egam*EV(ial, ikr, ixl, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
               varphi_a*(1d0-varphi_k)*(1d0-varphi_x)      *(egam*EV(ial, ikr, ixr, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
               (1d0-varphi_a)*varphi_k*varphi_x            *(egam*EV(iar, ikl, ixl, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
               (1d0-varphi_a)*varphi_k*(1d0-varphi_x)      *(egam*EV(iar, ikl, ixr, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
               (1d0-varphi_a)*(1d0-varphi_k)*varphi_x      *(egam*EV(iar, ikr, ixl, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam) + &
               (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*(egam*EV(iar, ikr, ixr, ip_p_com, iw_com, ie_com, is_com, ij_com+1))**(1d0/egam))**egam/egam

     ! calculate bequest part of the value function
    S_temp = 0d0 !(1d0-psi(is_com, ij_com+1))*mu_b*max(a_p + (1d0-xi)*k_p, 1d-13)**egam/egam

    ! calculate future part of the value function
    if (a_temp < 0d0) then
      inv_e = -1d-13**egam/egam*(1d0+abs(a_temp))
    else
      inv_e = - (psi(is_com, ij_com+1)*beta*EV_temp + S_temp)
    endif

  end function


  !##############################################################################
  ! FUNCTION cons_o
  !
  ! calculates the value function of a worker/entrpreneur
  !##############################################################################
  function cons_o(x_in)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ###############################################
    real*8, intent(in) :: x_in(:)
    real*8 :: cons_o

    !##### OTHER VARIABLES ######################################################
    real*8 :: Q_plus, ind_o, income, tomorrow, varphi_q, varphi_p
    integer :: iql, iqr, ipl, ipr

    ! define tomorrow's assets
    Q_plus  = x_in(1)

    ! define labor supply
    lab_com = max(x_in(2), 0d0)

    ! compute current occupation
    ind_o = abs(dble(ik_com > 0))

    ! calculate current income
    income = (1d0-ind_o)*w*eff(is_com, ij_com)*eta(iw_com, is_com)*lab_com + &
             ind_o*theta(ie_com, is_com)*(k(ik_com)**alpha*(eff(is_com, ij_com)*lab_com)**(1d0-alpha))**nu

    ! calculate income tax
    inctax_com = tarif(income)

    ! calculate capital tax
    captax_com = taur*r*max(a(ia_com)-xi*k(ik_com), 0d0)

    ! calculate pension contribution
    pencon_com = (1d0-(1d0-phi)*ind_o)*min(income, 2d0*ybar)

    ! available assets
    aas_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + (1d0-delta_k)*k(ik_com) + income + beq(is_com, ij_com) &
               - inctax_com - captax_com - taup*pencon_com

    ! calculate consumption
    cons_com = (aas_com - Q_plus)*pinv

    ! calculate future earning points
    p_plus_com = (p(ip_com)*dble(ij_com-1) + (1d0-(1d0-phi)*ind_o)*mu*(lambda + (1d0-lambda)*min(income/ybar, 2d0)))/dble(ij_com)

    ! derive interpolation weights
    call linint_Grow(Q_plus, Q_l, Q_u, Q_grow, NQ, iql, iqr, varphi_q)
    call linint_Equi(p_plus_com, p_l, p_u, NP, ipl, ipr, varphi_p)

    ! restrict values to grid just in case
    iql = min(iql, NQ)
    iqr = min(iqr, NQ)
    varphi_q = max(min(varphi_q, 1d0),0d0)

    ! restrict values to grid just in case
    ipl = min(ipl, NP)
    ipr = min(ipr, NP)
    varphi_p = max(min(varphi_p, 1d0),0d0)

    ! get next period value function
    tomorrow = 0d0

    if(varphi_q <= varphi_p) then
      tomorrow = (varphi_q           *(egam*S(io_p_com, iql, ik_com, ix_com, ipl, iw_com, ie_com, is_com, ij_com))**(1d0/egam) +  &
                  (varphi_p-varphi_q)*(egam*S(io_p_com, iqr, ik_com, ix_com, ipl, iw_com, ie_com, is_com, ij_com))**(1d0/egam) +  &
                  (1d0-varphi_p)     *(egam*S(io_p_com, iqr, ik_com, ix_com, ipr, iw_com, ie_com, is_com, ij_com))**(1d0/egam))**egam/egam
    else
      tomorrow = (varphi_p           *(egam*S(io_p_com, iql, ik_com, ix_com, ipl, iw_com, ie_com, is_com, ij_com))**(1d0/egam) +  &
                  (varphi_q-varphi_p)*(egam*S(io_p_com, iql, ik_com, ix_com, ipr, iw_com, ie_com, is_com, ij_com))**(1d0/egam) +  &
                  (1d0-varphi_q)     *(egam*S(io_p_com, iqr, ik_com, ix_com, ipr, iw_com, ie_com, is_com, ij_com))**(1d0/egam))**egam/egam
     endif

    ! calculate today's value function
    if(cons_com <= 0d0)then
      cons_o = -1d-13**egam/egam*(1d0+abs(cons_com))
    elseif(lab_com < 0d0) then
      cons_o = -1d-13**egam/egam*(1d0+abs(lab_com))
    elseif(lab_com >= 1d0) then
      cons_o = -1d-13**egam/egam*lab_com
    else
      cons_o = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + tomorrow)
    endif

  end function


  !##############################################################################
  ! FUNCTION cons_r
  !
  ! calculates the value function of a retiree
  !##############################################################################
  function cons_r(x_in)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ###############################################
    real*8, intent(in) :: x_in
    real*8 :: cons_r

    !##### OTHER VARIABLES ######################################################
    real*8 :: Q_plus, tomorrow, varphi_q
    integer :: iql, iqr

    ! define tomorrow's assets
    Q_plus  = x_in

    ! define labor supply
    lab_com = 0d0

    ! pension contribution
    pencon_com = 0d0

    ! calculate income tax
    inctax_com = tarif(pen(ip_com, ij_com) + ann(ix_com, is_com, ij_com))

    ! calculate capital tax
    captax_com = taur*r*a(ia_com)

    ! available assets
    aas_com = (1d0+r)*a(ia_com) + pen(ip_com, ij_com) + ann(ix_com, is_com, ij_com) &
              - inctax_com - captax_com

    ! calculate consumption
    cons_com = (aas_com - Q_plus)*pinv

    ! define future earning points
    p_plus_com = p(ip_com)

    ! derive interpolation weights
    call linint_Grow(Q_plus, Q_l, Q_u, Q_grow, NQ, iql, iqr, varphi_q)

    ! restrict values to grid just in case
    iql = min(iql, NQ)
    iqr = min(iqr, NQ)
    varphi_q = max(min(varphi_q, 1d0),0d0)

    ! get next period value function
    tomorrow = 0d0
    if (ij_com < JJ .or. mu_b /= 0d0) then
      tomorrow = (varphi_q      *(egam*S(io_p_com, iql, ik_com, ix_com, ip_com, iw_com, ie_com, is_com, ij_com))**(1d0/egam) +  &
                  (1d0-varphi_q)*(egam*S(io_p_com, iqr, ik_com, ix_com, ip_com, iw_com, ie_com, is_com, ij_com))**(1d0/egam))**egam/egam
    endif

    ! calculate today's value function
    if(cons_com <= 0d0)then
       cons_r = -1d-13**egam/egam*(1d0+abs(cons_com))
    else
       cons_r = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + tomorrow)
    endif

  end function


  !##############################################################################
  ! FUNCTION tr
  !
  ! calculates the transaction costs
  !##############################################################################
  function tr(k, k_p)

   implicit none

   !##### INPUT/OUTPUT VARIABLES ################################################
   real*8, intent(in) :: k, k_p
   real*8 :: tr

   tr = 0d0
   if (k <= 0d0) then
      tr = phi_k*(k_p - k)
   endif

  end function


  !##############################################################################
  ! FUNCTION tarif
  !
  ! Calculates income tax
  !##############################################################################
  function tarif(zn)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ###############################################
    real*8, intent(in) :: zn
    real*8 :: tarif

    !##### OTHER VARIABLES ######################################################
    real*8 :: tr2, tr3, znr

    ! determine parameters of tarif
    tr2 =       (r2-r1)*t1 + (r2-r1)**2d0*b1*0.5d0
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
  ! FUNCTION check_grid
  !
  ! Checks for the maximum gridpoints used
  !##############################################################################
  subroutine check_grid(iqmax, iamax, ikmax, ixmax)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ###############################################
    integer :: iqmax(JJ), iamax(JJ), ikmax(JJ), ixmax(JJ)

    !##### OTHER VARIABLES ######################################################
    integer :: iq, ia, ik, ix, ij

    iqmax = 0
    iamax = 0
    ikmax = 0
    ixmax = 0

    do ij = 1, JJ

      ! check for the maximum total asset grid point used at a certain age
      do iq = NQ, 0, -1
        if (sum(m_Q(iq, :, :, :, :, :, :, ij)) > 1d-10) then
          iqmax(ij) = iq
          exit
        endif
      enddo ! iq

      ! check for the maximum liquid asset grid point used at a certain age
      do ia = NA, 0, -1
        if (sum(m(ia, :, :, :, :, :, :, ij)) > 1d-10) then
          iamax(ij) = ia
          exit
        endif
      enddo ! ia

      ! check for the maximum investment grid point used at a certain age
      do ik = NK, 0, -1
        if (sum(m(:, ik, :, :, :, :, :, ij)) > 1d-10) then
          ikmax(ij) = ik
          exit
        endif
      enddo ! ik

      ! check for the maximum annuity grid point used at a certain age
      do ix = NX, 0, -1
        if (sum(m(:, :, ix, :, :, :, :, ij)) > 1d-10) then
          ixmax(ij) = ix
          exit
        endif
      enddo ! ix

    enddo ! ij


  end subroutine


  !##############################################################################
  ! FUNCTION sigma5
  !
  ! converts sigma of an one-year process into sigma of an five-year process
  !##############################################################################
  function sigma5(rho, sigma)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ###############################################
    real*8, intent(in) :: rho, sigma
    real*8 :: sigma5

    !##### OTHER VARIABLES ######################################################
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
  ! SUBROUTINE tick
  !
  ! sets the starting time
  !##############################################################################
  subroutine tick(t)

    !##### INPUT/OUTPUT VARIABLES ###############################################
    integer, intent(OUT) :: t

    ! save time
    call system_clock(t)

  end subroutine


  !##############################################################################
  ! SUBROUTINE tick
  !
  ! returns time in seconds from now to time described by t
  !##############################################################################
  subroutine tock(t)

    !##### INPUT/OUTPUT VARIABLES ###############################################
    integer, intent(in) :: t
    real*8 :: calctime

    !##### OTHER VARIABLES ######################################################
    integer :: now, clock_rate

    ! get time
    call system_clock(now, clock_rate)

    ! calculate time in seconds from now to time described by t
    calctime = real(now - t)/real(clock_rate)

    ! transform into hours and minutes
    if (calctime < 60) then
        write(*,'(/a,f7.3,a)')'Time elapsed: ', calctime, ' s'
    elseif (calctime < 3600) then
        write(*,'(/a,i3,a,f7.3,a)')'Time elapsed: ', floor(calctime/60d0), ' min ', mod(calctime, 60d0), ' s'
    else
        write(*,'(/a,i3,a,i3,a,f7.3,a)')'Time elapsed: ', floor(calctime/3600d0), ' h ', &
                                        floor(mod(calctime, 3600d0)/60d0), ' min ', mod(calctime, 60d0), ' s'
    endif

  end subroutine tock


end module
