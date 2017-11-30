module globals

    use toolbox

    implicit none

    ! number of years the household retires
    integer, parameter :: JR = 10

    ! number of years the household lives
    integer, parameter :: JJ = 16

    ! number of productivity (eta) shocks
    integer, parameter :: NW = 3

    ! number of entrepreneurial ability (theta) shocks
    integer, parameter :: NE = 3

    ! number of points on the asset grid
    integer, parameter :: NQ = 12

    ! number of points on the liquid asset grid
    integer, parameter :: NA = 12

    ! number of points on the capital grid
    integer, parameter :: NK = 12

    ! number of points on the annuity asset grid
    integer, parameter :: NX = 1

    ! number of points on the pension claim grid
    integer, parameter :: NP = 4

    ! demographic parameters
    real*8, parameter :: n_p = 0d0! (1d0+0.005d0)**5-1d0

    ! household preference parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: sigma = 0.3d0
    real*8, parameter :: beta = 0.96d0**5
    real*8, parameter :: mu_b = 0d0 !0.25d0

    ! maximum investment in annuities
    real*8, parameter :: mx_max = 0.10d0

    ! capital parameters
    real*8, parameter :: delta_k = 0d0 !0.06d0
    real*8, parameter :: xi = 1d0/3d0

    ! production parameters
    real*8, parameter :: k_min = 0.2d0
    real*8, parameter :: phi_k = 0d0!0.4d0
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: nu = 0.88d0

    ! size of the asset grid
    real*8, parameter :: Q_l    = 0d0
    real*8, parameter :: Q_u    = 10d0
    real*8, parameter :: Q_grow = 0.10d0

    ! size of the liquid asset grid
    real*8, parameter :: a_l    = Q_l
    real*8, parameter :: a_u    = Q_u
    real*8, parameter :: a_grow = Q_grow

    ! size of the capital grid
    real*8, parameter :: k_l = k_min
    real*8, parameter :: k_u = 0.5d0*Q_u/(1d0-xi)
    real*8, parameter :: k_grow = Q_grow

    ! size of the annuity grid
    real*8, parameter :: x_l    = Q_l
    real*8, parameter :: x_u    = Q_u
    real*8, parameter :: x_grow = Q_grow

    ! size of the pension claim grid
    real*8, parameter :: p_l    = 0d0
    real*8, parameter :: p_u    = 2d0

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.45d0

    ! numerical parameters
    integer, parameter :: itermax = 20
    real*8, parameter :: sig = 1d-6
    real*8, parameter :: damp = 0.50d0

    ! measure time
    integer :: time

    ! discretized shocks
    real*8 :: dist_eta(NW), pi_eta(NW, NW), eta(NW), dist_theta(NE), pi_theta(NE, NE), theta(NE)

    ! demographic and other model parameters
    real*8 :: eff(JJ), pen(0:NP, JJ), ann(0:NX, JJ), psi(JJ+1), rpop(0:JJ+1), workpop, b(JJ)

    ! government variables
    real*8 :: lambda, phi, mu
    real*8 :: taup

    ! macroeconomic variables
    real*8 :: r, w
    real*8 :: ybar
    real*8 :: AA, BQ, PBEN, PCON
    real*8 :: YY, YC, YE, CC, II, KK, KC, KE, LC

    ! production variables
    real*8 :: Omega

    ! cohort aggregate variables
    real*8 :: c_coh(0:1, JJ), y_coh(0:1, JJ), l_coh(0:1, JJ), o_coh(JJ)
    real*8 :: a_coh(0:1, JJ), x_coh(0:1, JJ), k_coh(JJ), penb_coh(JJ), penc_coh(JJ)

    ! different grids to discretize the state space
    real*8 :: Q(0:NQ), a(0:NA), k(0:NK), x(0:NX), p(0:NP)

    ! variables to store the policy functions
    real*8 :: Q_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: a_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), x_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: p_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), k_plus(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: penb(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), penc(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: c(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), l(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)

    ! variables for temporary policy and value functions
    real*8 :: Q_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: a_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), x_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: p_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), k_plus_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: penb_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), penc_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: c_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), l_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: V_t(0:1, 0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)

    ! variables to store the portfolio choice decisions
    real*8 :: omega_x_t(0:1, 0:NQ, 0:NK, 0:NX, 0:NP, NW, NE, JJ), omega_k_t(0:1, 0:NQ, 0:NK, 0:NX, 0:NP, NW, NE, JJ)
    real*8 :: S(0:1, 0:NQ, 0:NK, 0:NX, 0:NP, NW, NE, JJ)

    ! variables to store the value functions
    real*8 :: V(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ), EV(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)

    ! weights for the different gridpoints on the discretized state space
    real*8 :: m_Q(0:NQ, 0:NK, 0:NX, 0:NP, NW, NE, JJ), m(0:NA, 0:NK, 0:NX, 0:NP, NW, NE, JJ)

    ! numerical variables
    integer :: ij_com, iq_com, ia_com, ix_com, ip_com, ik_com, iw_com, ie_com, ia_p_com, iq_p_com, ip_p_com, io_p_com, iter
    integer :: iqmax(JJ), iamax(JJ), ixmax(JJ), ikmax(JJ)
    real*8 :: cons_com, lab_com, x_plus_com, p_plus_com, penc_com
    real*8 :: DIFF

    !$omp threadprivate(ij_com, iq_com, ia_com, ix_com, ip_com, ik_com, iw_com, ie_com, ia_p_com, ip_p_com, iq_p_com, io_p_com)
    !$omp threadprivate(cons_com, lab_com, x_plus_com, p_plus_com, penc_com)


  contains

    ! solve the household's decision of how much wealth to invest into capital
    subroutine solve_worker(iq_p, ik, ix, ip_p, iw, ie, ij)

      implicit none

      integer, intent(in) :: iq_p, ik, ix, ip_p, iw, ie, ij
      real*8 :: x_in, fret

      ! set up communication variables
      iq_p_com = iq_p; ik_com = ik; ix_com = ix; ip_p_com = ip_p; iw_com = iw; ie_com = ie; ij_com = ij

      if (Q(iq_p) > 0d0) then

         ! get best guess for the root of foc_real
         x_in = max(omega_x_t(0, iq_p, ik, ix, ip_p, iw, ie, ij), 1d-2)

         ! solve the household problem using fminsearch
         call fminsearch(x_in, fret, 0d0, 1d0, inv_w)

         ! portfolio share for capital
         omega_x_t(0, iq_p, ik, ix, ip_p, iw, ie, ij) = x_in
         omega_k_t(0, iq_p, ik, ix, ip_p, iw, ie, ij) = 0d0
         S(0, iq_p, ik, ix, ip_p, iw, ie, ij) = -fret

      else

        omega_x_t(0, iq_p, ik, ix, ip_p, iw, ie, ij) = 0d0
        omega_k_t(0, iq_p, ik, ix, ip_p, iw, ie, ij) = 0d0
        S(0, iq_p, ik, ix, ip_p, iw, ie, ij) = -inv_w(0d0)

      endif

    end subroutine

    ! solve the household's decision of how much wealth to invest into firm capital
    subroutine solve_entrepreneur(iq_p, ik, ix, ip_p, iw, ie, ij)

        implicit none

        integer, intent(in) :: iq_p, ik, ix, ip_p, iw, ie, ij
        real*8 :: x_in(2), fret

        ! set up communication variables
        iq_p_com = iq_p; ik_com = ik; ix_com = ix; ip_p_com = ip_p; iw_com = iw; ie_com = ie; ij_com = ij

        if (Q(iq_p) > (1d0-xi)*k_min) then

         ! get best guess for the root of foc_real
          x_in(1) = max(omega_x_t(1, iq_p, ik, ix, ip_p, iw, ie, ij), 1d-2)
          x_in(2) = max(omega_k_t(1, iq_p, ik, ix, ip_p, iw, ie, ij), 1d-2)

         ! solve the household problem using fminsearch
         call fminsearch(x_in, fret, (/0d0, 0d0/), (/1d0, 1d0/), inv_e)

         ! portfolio share for capital
         omega_x_t(1, iq_p, ik, ix, ip_p, iw, ie, ij) = x_in(1)
         omega_k_t(1, iq_p, ik, ix, ip_p, iw, ie, ij) = x_in(2)
         S(1, iq_p, ik, ix, ip_p, iw, ie, ij) = -fret

      else

        omega_x_t(1, iq_p, ik, ix, ip_p, iw, ie, ij) = 0d0
        omega_k_t(1, iq_p, ik, ix, ip_p, iw, ie, ij) = 0d0
        S(1, iq_p, ik, ix, ip_p, iw, ie, ij) = -inv_e((/1d0, 0d0/))

      endif

    end subroutine

    ! solve the household's decision of how much wealth to invest into capital
    subroutine solve_retiree(iq_p, ik, ix, ip_p, iw, ie, ij)

      implicit none

      integer, intent(in) :: iq_p, ik, ix, ip_p, iw, ie, ij
      integer :: ial, iar, ixl, ixr
      real*8 :: a_p, x_p, EV_temp, S_temp, varphi_a, varphi_x

      a_p = Q(iq_p)
      x_p = x(ix)

     ! calculate linear interpolation for future assets
     call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
     call linint_Grow(x_p, x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)

     ! restrict values to grid just in case
     ial = min(ial, NA)
     iar = min(iar, NA)
     varphi_a = max(min(varphi_a, 1d0),0d0)

     ! restrict values to grid just in case
     ixl = min(ixl, NX)
     ixr = min(ixr, NX)
     varphi_x = max(min(varphi_x, 1d0),0d0)

     ! calculate future part of the value function
     S_temp = (1d0-psi(ij+1))*mu_b*max(Q(iq_p), 1d-13)**egam/egam

    if (varphi_a <= varphi_x) then
       EV_temp = (varphi_a             *(egam*EV(ial, 0, ixl, ip_p, iw, ie, ij+1))**(1d0/egam) + &
                  (varphi_x - varphi_a)*(egam*EV(iar, 0, ixl, ip_p, iw, ie, ij+1))**(1d0/egam) + &
                  (1d0-varphi_x)       *(egam*EV(iar, 0, ixr, ip_p, iw, ie, ij+1))**(1d0/egam))**egam/egam
      else
        EV_temp = (varphi_x             *(egam*EV(ial, 0, ixl, ip_p, iw, ie, ij+1))**(1d0/egam) + &
                   (varphi_a - varphi_x)*(egam*EV(ial, 0, ixr, ip_p, iw, ie, ij+1))**(1d0/egam) + &
                   (1d0-varphi_a)       *(egam*EV(iar, 0, ixr, ip_p, iw, ie, ij+1))**(1d0/egam))**egam/egam
      endif

     omega_x_t(:, iq_p, ik, ix, ip_p, iw, ie, ij) = 0d0
     omega_k_t(:, iq_p, ik, ix, ip_p, iw, ie, ij) = 0d0
     S(:, iq_p, ik, ix, ip_p, iw, ie, ij) = psi(ij+1)*EV_temp + S_temp

    end subroutine


  ! solve the household's consumption-savings decision
  subroutine solve_consumption(io_p, ia, ik, ix, ip, iw, ie, ij)

      implicit none

      integer, intent(in) :: io_p, ia, ik, ix, ip, iw, ie, ij
      real*8 :: x_in(2), fret, x_p, mx, k_p, varphi_q, varphi_p
      integer :: iql, iqr, ipl, ipr

      ! set up communication variables
      io_p_com = io_p; ia_com = ia; ik_com = ik; ix_com = ix; ip_com = ip; iw_com = iw; ie_com = ie; ij_com = ij

      ! get best initial guess from future period
      x_in(1) = max(Q_plus_t(io_p, ia, ik, ix, ip, iw, ie, ij), 1d-2)
      x_in(2) = max(l_t(io_p, ia, ik, ix, ip, iw, ie, ij), 0.33d0)

      ! solve the household problem using fminsearch
      if (ij < JR) then
        call fminsearch(x_in, fret, (/Q_l, 0d0/), (/Q_u, 0.99d0/), cons_o)
      else
        call fminsearch(x_in(1), fret, Q_l, Q_u, cons_r)
      endif

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

       ! get next period's capital size
       if (varphi_q <= varphi_p) then
         k_p = ((1d0-xi)*k_min + (varphi_q             *omega_k_t(io_p, iql, ik, ix, ipl, iw, ie, ij) +  &
                                  (varphi_p-varphi_q)  *omega_k_t(io_p, iqr, ik, ix, ipl, iw, ie, ij) +  &
                                  (1d0-varphi_p)       *omega_k_t(io_p, iqr, ik, ix, ipr, iw, ie, ij))*(x_in(1)-(1d0-xi)*k_min))/(1d0-xi)
       else
         k_p = ((1d0-xi)*k_min + (varphi_p             *omega_k_t(io_p, iql, ik, ix, ipl, iw, ie, ij) +  &
                                  (varphi_q-varphi_p)  *omega_k_t(io_p, iql, ik, ix, ipr, iw, ie, ij) +  &
                                  (1d0-varphi_q)       *omega_k_t(io_p, iqr, ik, ix, ipr, iw, ie, ij))*(x_in(1)-(1d0-xi)*k_min))/(1d0-xi)
       endif

      endif

      ! determine future annuity stock
      x_p = 0d0
      mx = 0d0

      if (ij < JR) then

        if (varphi_q <= varphi_p) then
          mx = (varphi_q           *omega_x_t(io_p, iql, ik, ix, ipl, iw, ie, ij) +  &
                (varphi_p-varphi_q)*omega_x_t(io_p, iqr, ik, ix, ipl, iw, ie, ij) +  &
                (1d0-varphi_p)     *omega_x_t(io_p, iqr, ik, ix, ipr, iw, ie, ij))*x_in(1)
        else
          mx = (varphi_p           *omega_x_t(io_p, iql, ik, ix, ipl, iw, ie, ij) +  &
                (varphi_q-varphi_p)*omega_x_t(io_p, iql, ik, ix, ipr, iw, ie, ij) +  &
                (1d0-varphi_q)     *omega_x_t(io_p, iqr, ik, ix, ipr, iw, ie, ij))*x_in(1)
        endif

        x_p = (1d0+r)/psi(ij)*x(ix) + mx

      else

        x_p = x(ix)

      endif

      if (mx>1d-12) write(*,*) mx, x_p, io_p, ia, ik, ix, ip, iw, ie, ij
      if (mx>1d-12) write(*,*) omega_x_t(io_p, iql, ik, ix, ipl, iw, ie, ij), omega_x_t(io_p, iql, ik, ix, ipr, iw, ie, ij), omega_x_t(io_p, iqr, ik, ix, ipr, iw, ie, ij)
      if (mx>1d-12) write(*,*) omega_k_t(io_p, iql, ik, ix, ipl, iw, ie, ij), omega_k_t(io_p, iql, ik, ix, ipr, iw, ie, ij), omega_k_t(io_p, iqr, ik, ix, ipr, iw, ie, ij)

      ! copy decisions
      Q_plus_t(io_p, ia, ik, ix, ip, iw, ie, ij) = x_in(1)
      a_plus_t(io_p, ia, ik, ix, ip, iw, ie, ij) = x_in(1) - (1d0-xi)*k_p - mx
      k_plus_t(io_p, ia, ik, ix, ip, iw, ie, ij) = k_p
      x_plus_t(io_p, ia, ik, ix, ip, iw, ie, ij) = x_p
      p_plus_t(io_p, ia, ik, ix, ip, iw, ie, ij) = p_plus_com
      penb_t(io_p, ia, ik, ix, ip, iw, ie, ij) = pen(ip, ij)
      penc_t(io_p, ia, ik, ix, ip, iw, ie, ij) = penc_com
      c_t(io_p, ia, ik, ix, ip, iw, ie, ij) =  (1d0+r)*a(ia) + pen(ip, ij) + ann(ix, ij) - Q_plus_t(io_p, ia, ik, ix, ip, iw, ie, ij)
      l_t(io_p, ia, ik, ix, ip, iw, ie, ij) = lab_com
      V_t(io_p, ia, ik, ix, ip, iw, ie, ij) = -fret

      !if (ij>JR) write(*,*) Q_plus_t(io_p, ia, ik, ix, ip, iw, ie, ij) + c_t(io_p, ia, ik, ix, ip, iw, ie, ij) - (1d0+r)*a(ia) - pen(ip, ij) - ann(ix, ij)

  end subroutine

  ! the first order condition with respect to next period real estate
  function inv_w(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: inv_w, a_p, x_p, EV_temp, S_temp, omega_x, varphi_a, varphi_x, a_temp
      integer :: ial, iar, ixl, ixr

      ! store real estate share
      omega_x  = x_in

      ! determine future liquid wealth and future downpayment
      x_p = (1d0+r)/psi(ij_com)*x(ix_com) + min(omega_x*Q(iq_p_com), mx_max)
      a_temp = (1d0-omega_x)*Q(iq_p_com)
      a_p = max(a_temp, 0d0)

      ! derive interpolation weights
      call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
      call linint_Grow(x_p, x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)

      ! restrict values to grid just in case
      ial = min(ial, NA)
      iar = min(iar, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)

      ! restrict values to grid just in case
      ixl = min(ixl, NX)
      ixr = min(ixr, NX)
      varphi_x = max(min(varphi_x, 1d0),0d0)

      S_temp = 0d0 !(1d0-psi(ij_com+1))*mu_b*max((1d0-omega_x)*Q(iq_p_com), 1d-13)**egam/egam

      ! get optimal investment strategy
      if (varphi_a <= varphi_x) then
        EV_temp = (varphi_a            *(egam*EV(ial, 0, ixl, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                   (varphi_x-varphi_a) *(egam*EV(iar, 0, ixl, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                   (1d0-varphi_x)      *(egam*EV(iar, 0, ixr, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam))**egam/egam
      else
        EV_temp = (varphi_x            *(egam*EV(ial, 0, ixl, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                   (varphi_a-varphi_x) *(egam*EV(ial, 0, ixr, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                   (1d0-varphi_a)      *(egam*EV(iar, 0, ixr, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam))**egam/egam
      endif

      if (a_temp < 0d0) then
        inv_w = -1d-13**egam/egam*(1d0+abs(a_temp))
      else
        inv_w = - (psi(ij_com+1)*EV_temp + S_temp)
      endif

  end function

  ! the first order condition with respect to next period real estate
  function inv_e(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in(:)

      ! variable declarations
      real*8 :: inv_e, a_p, x_p, k_p, EV_temp, S_temp, omega_x, omega_k, varphi_a, varphi_k, varphi_x, a_temp
      integer :: ial, iar, ikl, ikr, ixl, ixr

      ! store real estate share
      omega_x  = x_in(1)
      omega_k  = x_in(2)

      ! determine future liquid wealth and future downpayment
      x_p = (1d0+r)/psi(ij_com)*x(ix_com) + min(omega_x*Q(iq_p_com), mx_max)
      k_p = ((1d0-xi)*k_min + omega_k*(Q(iq_p_com)-(1d0-xi)*k_min))/(1d0-xi)
      a_temp = Q(iq_p_com) - omega_x*Q(iq_p_com) - (1d0-xi)*k_p - tr(k(ik_com), k_p)
      a_p = max(a_temp, 0d0)

      ! derive interpolation weights
      call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
      call linint_Grow(x_p, x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)
      call linint_Grow(k_p, k_l, k_u, k_grow, NK-1, ikl, ikr, varphi_k)

      ! restrict values to grid just in case
      ial = min(ial, NA)
      iar = min(iar, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)

      ! restrict values to grid just in case
      ixl = min(ixl, NX)
      ixr = min(ixr, NX)
      varphi_x = max(min(varphi_x, 1d0),0d0)

      ! restrict values to grid just in case
      ikl = min(ikl+1, NK)
      ikr = min(ikr+1, NK)
      varphi_k = max(min(varphi_k, 1d0), 0d0)

      S_temp = 0d0 !(1d0-psi(ij_com+1))*mu_b*max((1d0-omega_x)*Q(iq_p_com), 1d-13)**egam/egam

      ! get optimal investment strategy
      EV_temp = (varphi_a*varphi_x*varphi_k                  *(egam*EV(ial, ikl, ixl, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                 varphi_a*varphi_k*(1d0-varphi_x)            *(egam*EV(ial, ikl, ixr, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                 varphi_a*(1d0-varphi_k)*varphi_x            *(egam*EV(ial, ikr, ixl, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                 varphi_a*(1d0-varphi_k)*(1d0-varphi_x)      *(egam*EV(ial, ikr, ixr, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                 (1d0-varphi_a)*varphi_k*varphi_x            *(egam*EV(iar, ikl, ixl, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                 (1d0-varphi_a)*varphi_k*(1d0-varphi_x)      *(egam*EV(iar, ikl, ixr, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                 (1d0-varphi_a)*(1d0-varphi_k)*varphi_x      *(egam*EV(iar, ikr, ixl, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam) + &
                 (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*(egam*EV(iar, ikr, ixr, ip_p_com, iw_com, ie_com, ij_com+1))**(1d0/egam))**egam/egam

      if (a_temp < 0d0) then
        inv_e = -1d-13**egam/egam*(1d0+abs(a_temp))
      else
        inv_e = - (psi(ij_com+1)*EV_temp + S_temp)
      endif

  end function

    ! the first order condition regarding consumption
    function cons_o(x_in)

        implicit none

        ! input variable
        real*8, intent(in) :: x_in(:)

        ! variable declarations
        real*8 :: cons_o, Q_plus, ind_o, income, tomorrow, varphi_q, varphi_p
        integer :: iql, iqr, ipl, ipr

        ! define tomorrow's assets
        Q_plus  = x_in(1)

        ! define labor supply
        lab_com = max(x_in(2), 0d0)

        ! compute current occupation
        ind_o = abs(dble(ik_com > 0))

        ! calculate current income
        income = (1d0-ind_o)*w*eff(ij_com)*eta(iw_com)*lab_com + &
                 ind_o*theta(ie_com)*(k(ik_com)**alpha*(eff(ij_com)*lab_com)**(1d0-alpha))**nu

        ! pension contribution
        penc_com = (1d0-(1d0-phi)*ind_o)*min(income, p_u)

        ! calculate consumption-savings
        cons_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + (1d0-delta_k)*k(ik_com) + income + b(ij_com) &
                   - taup*penc_com - Q_plus

        ! calculate future earning points
        p_plus_com = (p(ip_com)*dble(ij_com-1) + (1d0-(1d0-phi)*ind_o)*mu*(lambda+(1d0-lambda)*min(income, p_u)))/dble(ij_com)

        ! calculate linear interpolation for future part of value function
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
        if (ij_com < JJ .or. mu_b /= 0d0) then

          if(varphi_q <= varphi_p) then
            tomorrow = (varphi_q           *(egam*S(io_p_com, iql, ik_com, ix_com, ipl, iw_com, ie_com, ij_com))**(1d0/egam) +  &
                        (varphi_p-varphi_q)*(egam*S(io_p_com, iqr, ik_com, ix_com, ipl, iw_com, ie_com, ij_com))**(1d0/egam) +  &
                        (1d0-varphi_p)     *(egam*S(io_p_com, iqr, ik_com, ix_com, ipr, iw_com, ie_com, ij_com))**(1d0/egam))**egam/egam
          else
            tomorrow = (varphi_p           *(egam*S(io_p_com, iql, ik_com, ix_com, ipl, iw_com, ie_com, ij_com))**(1d0/egam) +  &
                        (varphi_q-varphi_p)*(egam*S(io_p_com, iql, ik_com, ix_com, ipr, iw_com, ie_com, ij_com))**(1d0/egam) +  &
                        (1d0-varphi_q)     *(egam*S(io_p_com, iqr, ik_com, ix_com, ipr, iw_com, ie_com, ij_com))**(1d0/egam))**egam/egam
           endif

        endif

        ! calculate today's value function
        if(cons_com <= 0d0)then
           cons_o = -1d-13**egam/egam*(1d0+abs(cons_com))
        elseif(lab_com < 0d0) then
          cons_o = -1d-13**egam/egam*(1d0+abs(lab_com))
        elseif(lab_com >= 1d0) then
          cons_o = -1d-13**egam/egam*lab_com
        else
           cons_o = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + beta*tomorrow)
        endif

    end function


    ! the first order condition regarding consumption
    function cons_r(x_in)

        implicit none

        ! input variable
        real*8, intent(in) :: x_in

        ! variable declarations
        real*8 :: cons_r, Q_plus, tomorrow, varphi_q
        integer :: iql, iqr

        ! define tomorrow's assets
        Q_plus  = x_in

        ! define labor supply
        lab_com = 0d0

        ! pension contribution
        penc_com = 0d0

        ! calculate consumption
        cons_com = (1d0+r)*a(ia_com) + pen(ip_com, ij_com) + ann(ix_com, ij_com) &
                    - Q_plus

        ! define future earning points
        p_plus_com = p(ip_com)

        ! calculate linear interpolation for future part of value function
        call linint_Grow(Q_plus, Q_l, Q_u, Q_grow, NQ, iql, iqr, varphi_q)

        ! restrict values to grid just in case
        iql = min(iql, NQ)
        iqr = min(iqr, NQ)
        varphi_q = max(min(varphi_q, 1d0),0d0)

        ! get next period value function
        tomorrow = 0d0
        if (ij_com < JJ .or. mu_b /= 0d0) then
            tomorrow = (varphi_q           *(egam*S(io_p_com, iql, ik_com, ix_com, ip_com, iw_com, ie_com, ij_com))**(1d0/egam) +  &
                        (1d0-varphi_q)     *(egam*S(io_p_com, iqr, ik_com, ix_com, ip_com, iw_com, ie_com, ij_com))**(1d0/egam))**egam/egam
        endif

        ! calculate today's value function
        if(cons_com <= 0d0)then
           cons_r = -1d-13**egam/egam*(1d0+abs(cons_com))
        else
           cons_r = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + beta*tomorrow)
        endif

    end function

    ! function that defines adjustment cost
    function tr(k, k_p)

       implicit none

       ! input variable
       real*8, intent(in) :: k, k_p
       real*8 :: tr

       tr = 0d0
       if (k <= 0d0) then

          tr = phi_k*(k_p - k)

       endif

    end function

    !##############################################################################
    ! FUNCTION check_grid
    !
    ! Checks for the maximum gridpoints used
    !##############################################################################
    subroutine check_grid(iqmax, iamax, ikmax, ixmax)

      implicit none

      !##### INPUT/OUTPUT VARIABLES #############################################
      integer :: iqmax(JJ), iamax(JJ), ikmax(JJ), ixmax(JJ)

      !##### OTHER VARIABLES ####################################################
      integer :: iq, ia, ik, ix, ij

      iqmax = 0
      iamax = 0
      ikmax = 0
      ixmax = 0

      do ij = 1, JJ

        ! check for the maximum total asset grid point used at a certain age
        do iq = NQ, 0, -1
          if (sum(m_Q(iq, :, :, :, :, :, ij)) > 1d-10) then
            iqmax(ij) = iq
            exit
          endif
        enddo ! iq

        ! check for the maximum liquid asset grid point used at a certain age
        do ia = NA, 0, -1
          if (sum(m(ia, :, :, :, :, :, ij)) > 1d-10) then
            iamax(ij) = ia
            exit
          endif
        enddo ! ia

        ! check for the maximum investment grid point used at a certain age
        do ik = NK, 0, -1
          if (sum(m(:, ik, :, :, :, :, ij)) > 1d-10) then
            ikmax(ij) = ik
            exit
          endif
        enddo ! ik

        ! check for the maximum annuity grid point used at a certain age
        do ix = NX, 0, -1
          if (sum(m(:, :, ix, :, :, :, ij)) > 1d-10) then
            ixmax(ij) = ix
            exit
          endif
        enddo ! ix

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

    subroutine tick(t)

      integer, intent(OUT) :: t

      call system_clock(t)

    end subroutine tick

    ! returns time in seconds from now to time described by t
    subroutine tock(t)

      integer, intent(in) :: t
      real*8 :: calctime
      integer :: now, clock_rate

      call system_clock(now, clock_rate)

      calctime = real(now - t)/real(clock_rate)

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
