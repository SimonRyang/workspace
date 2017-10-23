module globals

    use toolbox

    implicit none

    ! number of years the household retires
    integer, parameter :: JR = 10

    ! number of years the household lives
    integer, parameter :: JJ = 16

    ! number of productivity (eta) shocks
    integer, parameter :: NW = 5

    ! number of entrepreneurial ability (theta) shocks
    integer, parameter :: NE = 5

    ! number of points on the asset grid
    integer, parameter :: NQ = 32

    ! number of points on the liquid asset grid
    integer, parameter :: NA = 32

    ! number of points on the annuity asset grid
    integer, parameter :: NX = 32

    ! number of points on the capital grid
    integer, parameter :: NK = 32

    ! number of points on the pension claim grid
    integer, parameter :: NP = 4

    ! number of occupation states
    integer, parameter :: NO = 1

    ! household preference parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: sigma = 0.3d0
    real*8, parameter :: beta = 0.99d0**5d0
    real*8, parameter :: mu_b = 0.15d0
    ! convert variables into per period values

    ! risk free rate and risk premium
    real*8, parameter :: r  = 0.1d0

    ! capital parameters
    real*8, parameter :: delta_k = 0.06d0
    real*8, parameter :: xi = 1d0/3d0

    ! production parameters
    real*8, parameter :: k_min = 0.6d0
    real*8, parameter :: phi_k = 0.4d0
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: nu = 0.88d0

    ! size of the asset grid
    real*8, parameter :: Q_l    = 0d0
    real*8, parameter :: Q_u    = 16d0
    real*8, parameter :: Q_grow = 0.05d0

    ! size of the liquid asset grid
    real*8, parameter :: a_l    = Q_l
    real*8, parameter :: a_u    = Q_u
    real*8, parameter :: a_grow = Q_grow

    ! size of the annuity grid
    real*8, parameter :: x_l    = Q_l
    real*8, parameter :: x_u    = Q_u
    real*8, parameter :: x_grow = Q_grow

    ! size of the pension claim grid
    real*8, parameter :: p_l    = 0d0
    real*8, parameter :: p_u    = 2d0

    ! size of the capital grid
    real*8, parameter :: k_l = 0d0
    real*8, parameter :: k_u = 0.75d0*Q_u/(1d0-xi)
    real*8, parameter :: k_grow = Q_grow

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.5d0

    ! measure time
    integer :: time

    ! discretized shocks
    real*8 :: dist_eta(NW), pi_eta(NW, NW), eta(NW), dist_theta(NE), pi_theta(NE, NE), theta(NE)

    ! wages, transfer payments (old-age), survival probabilities and discount factor for housing utilty
    real*8 :: w, eff(JJ), pen(JJ, 0:NP), p_hat(JJ), psi(JJ+1)

    ! government variables
    real*8 :: lambda, phi, mu
    real*8 :: taup

    ! cohort aggregate variables
    real*8 :: c_coh(JJ, 0:1), y_coh(JJ, 0:1), l_coh(JJ, 0:1), o_coh(JJ)
    real*8 :: a_coh(JJ, 0:1), k_coh(JJ)

    ! different grids to discretize the state space
    real*8 :: Q(0:NQ), a(0:NA), x(0:NX), p(0:NP), k(0:NK)

    ! variables to store the policy functions
    real*8 :: Q_plus(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE), a_plus(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE)
    real*8 :: p_plus(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE), k_plus(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE)
    real*8 :: c(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE), l(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE)

    ! variables for temporary policy and value functions
    real*8 :: Q_plus_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO)
    real*8 :: a_plus_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO), x_plus_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO)
    real*8 :: p_plus_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO), k_plus_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO)
    real*8 :: c_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO), l_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO)
    real*8 :: V_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO)

    ! variables to store the portfolio choice decisions
    real*8 :: omega_x(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE), omega_k(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE)
    real*8 :: omega_x_t(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO)

    ! variables to store the value functions
    real*8 :: V(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE), EV(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE), S(JJ, 0:NQ, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO)

    ! weights for the different gridpoints on the discretized state space
    real*8 :: m(JJ, 0:NA, 0:NX, 0:NP, 0:NK, NW, NE)

    ! numerical variables
    integer :: ij_com, iq_com, ia_com, ix_com, ip_com, ik_com, iw_com, ie_com, ia_p_com, iq_p_com, ip_p_com, io_p_com
    integer :: iamax(JJ), ixmax(JJ), ikmax(JJ)
    real*8 :: cons_com, lab_com, x_plus_com, p_plus_com

    !$omp threadprivate(ij_com, iq_com, ia_com, ix_com, ip_com, ik_com, iw_com, ie_com, ia_p_com, ip_p_com, iq_p_com, io_p_com)
    !$omp threadprivate(cons_com, lab_com, x_plus_com, p_plus_com)


  contains

    ! solve the household's decision of how much wealth to invest into capital
    subroutine solve_worker(ij, iq_p, ix, ip_p, ik, iw, ie)

      implicit none

      integer, intent(in) :: ij, iq_p, ix, ip_p, ik, iw, ie
      integer :: ial, iar, ixl, ixr
      real*8 :: x_in, fret, a_p, x_p, EV_temp, S_temp, varphi_a

      ! set up communication variables
      ij_com = ij; iq_p_com = iq_p; ix_com = ix; ip_p_com = ip_p; ik_com = ik; iw_com = iw; ie_com = ie

       ! get best guess for the root of foc_real
       x_in = max(omega_x(ij, iq_p, ix, ip_p, ik, iw, ie), 1d-4)

       ! solve the household problem using fminsearch
       call fminsearch(x_in, fret, 0d0, 1d0, inv_w)

       ! portfolio share for capital
       omega_x_t(ij, iq_p, ix, ip_p, ik, iw, ie, 0) = x_in
       S(ij, iq_p, ix, ip_p, ik, iw, ie, 0) = -fret

    end subroutine

    ! solve the household's decision of how much wealth to invest into firm capital
    subroutine solve_entrepreneur(ij, iq_p, ix, ip_p, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, iq_p, ix, ip_p, ik, iw, ie
        real*8 :: x_in(2), fret

        ! set up communication variables
        ij_com = ij; iq_p_com = iq_p; ip_p_com = ip_p; ik_com = ik; iw_com = iw; ie_com = ie

       ! get best guess for the root of foc_real
        x_in(1) = max(omega_x(ij+1, iq_p, ix, ip_p, ik, iw, ie), 1d-4)
        x_in(2) = max(omega_k(ij+1, iq_p, ix, ip_p, ik, iw, ie), 1d-4)

       ! solve the household problem using fminsearch
       call fminsearch(x_in, fret, (/0d0, 0d0/), (/1d0, 1d0/), inv_e)

       ! portfolio share for capital
       omega_x_t(ij, iq_p, ix, ip_p, ik, iw, ie, 1) = x_in(1)
       omega_k(ij, iq_p, ix, ip_p, ik, iw, ie) = x_in(2)
       S(ij, iq_p, ix, ip_p, ik, iw, ie, 1) = -fret

    end subroutine

    ! solve the household's decision of how much wealth to invest into capital
    subroutine solve_retiree(ij, iq_p, ix, ip_p, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, iq_p, ix, ip_p, ik, iw, ie
        integer :: ial_p, iar_p, ixl_p, ixr_p
        real*8 :: a_p, x_p, EV_temp, S_temp, varphi_a, varphi_x

        a_p  = Q(iq_p)
        x_p = (1d0+r)/psi(ij)*(1d0-p_hat(ij))*x(ix)

       ! calculate linear interpolation for future assets
       call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial_p, iar_p, varphi_a)
       call linint_Grow(x_p, x_l, x_u, x_grow, NX, ixl_p, ixr_p, varphi_x)

       ! restrict values to grid just in case
       ial_p = min(ial_p, NA)
       iar_p = min(iar_p, NA)
       varphi_a = max(min(varphi_a, 1d0),0d0)

       ixl_p = min(ixl_p, NX)
       ixr_p = min(ixr_p, NX)
       varphi_x = max(min(varphi_x, 1d0),0d0)

       S_temp = (1d0-psi(ij+1))*mu_b*max(Q(iq_p), 1d-16)**egam/egam

       if (varphi_a <= varphi_k) then
         EV_temp = varphi_a             *(egam*EV(ij+1, ial_p, ixl_p, ip_p, 0, iw, ie))**(1d0/egam) + &
                   (varphi_k - varphi_a)*(egam*EV(ij+1, iar_p, ixl_p, ip_p, 0, iw, ie))**(1d0/egam) + &
                   (1d0-varphi_k)       *(egam*EV(ij+1, iar_p, ixr_p, ip_p, 0, iw, ie))**(1d0/egam)
        else
          EV_temp = varphi_k             *(egam*EV(ij+1, ial_p, ixl_p, ip_p, 0, iw, ie))**(1d0/egam) + &
                    (varphi_a - varphi_k)*(egam*EV(ij+1, ial_p, ixr_p, ip_p, 0, iw, ie))**(1d0/egam) + &
                    (1d0-varphi_a)       *(egam*EV(ij+1, iar_p, ixr_p, ip_p, 0, iw, ie))**(1d0/egam)
        endif

       S_temp = S_temp + psi(ij+1)*EV_temp**egam/egam

       S(ij, iq_p, ix, ip_p, ik, iw, ie, 0) = S_temp

    end subroutine


  ! solve the household's consumption-savings decision
  subroutine solve_consumption(ij, ia, ix, ip, ik, iw, ie, io_p)

      implicit none

      integer, intent(in) :: ij, ia, ix, ip, ik, iw, ie, io_p
      real*8 :: x_in(2), fret, varphi_q, varphi_p, k_p
      integer :: iql_p, iqr_p, ipl_p, ipr_p

      ! set up communication variables
      ij_com = ij; ia_com = ia; ix_com = ix; ip_com = ip; ik_com = ik; iw_com = iw; ie_com = ie; io_p_com = io_p

      ! get best initial guess from future period
      x_in(1) = max(Q_plus_t(ij, ia, ix, ip, ik, iw, ie, io_p), 1d-4)
      x_in(2) = max(l_t(ij, ia, ix, ip, ik, iw, ie, io_p), 0.33d0)

      ! solve the household problem using rootfinding
      call fminsearch(x_in, fret, (/Q_l, 0d0/), (/Q_u, 0.99d0/), cons_o)

      ! determine future investment
      k_p = 0d0

      if (io_p == 1) then

        call linint_Grow(x_in(1), Q_l, Q_u, Q_grow, NQ, iql_p, iqr_p, varphi_q)
        call linint_Equi(p_plus_com, p_l, p_u, NP, ipl_p, ipr_p, varphi_p)

        ! restrict values to grid just in case
        iql_p = min(iql_p, NQ)
        iqr_p = min(iqr_p, NQ)
        varphi_q = max(min(varphi_q, 1d0),0d0)

        ipl_p = min(ipl_p, NP)
        ipr_p = min(ipr_p, NP)
        varphi_p = max(min(varphi_p, 1d0),0d0)

        ! get next period's capital size
        if (varphi_q <= varphi_p) then
          k_p = ((1d0-xi)*k_min + (varphi_q             *omega_k(ij, iql_p, ix, ipl_p, ik, iw, ie) +  &
                                   (varphi_p-varphi_q)  *omega_k(ij, iqr_p, ix, ipl_p, ik, iw, ie) +  &
                                   (1d0-varphi_p)       *omega_k(ij, iqr_p, ix, ipr_p, ik, iw, ie))*(x_in(1)-(1d0-xi)*k_min))/(1d0-xi)
        else
          k_p = ((1d0-xi)*k_min + (varphi_p             *omega_k(ij, iql_p, ix, ipl_p, ik, iw, ie) +  &
                                   (varphi_q-varphi_p)  *omega_k(ij, iql_p, ix, ipr_p, ik, iw, ie) +  &
                                   (1d0-varphi_q)       *omega_k(ij, iqr_p, ix, ipr_p, ik, iw, ie))*(x_in(1)-(1d0-xi)*k_min))/(1d0-xi)
        endif

      endif

      ! copy decisions
      Q_plus_t(ij, ia, ip, ik, iw, ie, io_p) = x_in(1)
      a_plus_t(ij, ia, ip, ik, iw, ie, io_p) = x_in(1) - (1d0-xi)*k_p
      k_plus_t(ij, ia, ip, ik, iw, ie, io_p) = k_p
      c_t(ij, ia, ip, ik, iw, ie, io_p) = cons_com
      l_t(ij, ia, ip, ik, iw, ie, io_p) = lab_com
      V_t(ij, ia, ip, ik, iw, ie, io_p) = -fret

  end subroutine

  ! the first order condition with respect to next period real estate
  function inv_w(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: inv_w, a_p, x_p, k_p, EV_temp, S_temp, omega_k, varphi_a, varphi_x, a_temp
      integer :: ial_p, iar_p, ixl_p, ixr_p

      ! store real estate share
      omega_x  = x_in

      ! determine future liquid wealth and future downpayment
      x_p = (1d0+r)/psi(ij_com)*x(ix_com) + omega_x*Q(iq_p_com)
      k_p = 0d0
      a_temp = Q(iq_p_com) - omega_x*Q(iq_p_com)
      a_p = max(a_temp, 0d0)

      ! derive interpolation weights
      call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial_p, iar_p, varphi_a)
      call linint_Grow(x_p, x_l, x_u, x_grow, NX, ixl_p, ixr_p, varphi_x)

      ! restrict values to grid just in case
      ial_p = min(ial_p, NA)
      iar_p = min(iar_p, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)

      ixl_p = min(ixl_p, NX)
      ixr_p = min(ixr_p, NX)
      varphi_x = max(min(varphi_x, 1d0),0d0)

      S_temp = (1d0-psi(ij_com+1))*mu_b*max((1d0-omega_x)*Q(iq_p_com), 1d-16)**egam/egam

      ! get optimal investment strategy
      if (varphi_a <= varphi_x) then
        EV_temp = varphi_a            *(egam*EV(ij_com+1, ial_p, 0, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                  (varhpi_x-varphi_a) *(egam*EV(ij_com+1, iar_p, 0, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                  varphi_x            *(egam*EV(ij_com+1, iar_p, 0, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam)
      else
        EV_temp = varphi_x            *(egam*EV(ij_com+1, ial_p, 0, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                  (varhpi_a-varphi_x) *(egam*EV(ij_com+1, ial_p, 0, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam) + &
                  varphi_y            *(egam*EV(ij_com+1, iar_p, 0, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam)
      endif

      S_temp = S_temp + psi(ij_com+1)*EV_temp**egam/egam

      inv_w = - (S_temp + 1d-16**egam/egam*abs(a_p-a_temp))

  end function

  ! the first order condition with respect to next period real estate
  function inv_e(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in(:)

      ! variable declarations
      real*8 :: inv_e, a_p, x_p, k_p, EV_temp, S_temp, omega_k, varphi_a, varphi_x, varphi_k, a_temp
      integer :: ial_p, iar_p, ixl_p, ixr_p, ikl_p, ikr_p

      ! store real estate share
      omega_x  = x_in(1)
      omega_k  = x_in(2)

      ! determine future liquid wealth and future downpayment
      x_p = (1d0+r)/psi(ij_com)*x(ix_com) + omega_x*Q(iq_p_com)
      k_p =((1d0-xi)*k_min + omega_k*(Q(iq_p_com)-(1d0-xi)*k_min))/(1d0-xi)
      a_temp = Q(iq_p_com) - omega_x*Q(iq_p_com) - (1d0-xi)*k_p - tr(k(ik_com), k_p)
      a_p = max(a_temp, 0d0)

      ! derive interpolation weights
      call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial_p, iar_p, varphi_a)
      call linint_Grow(x_p, x_l, x_u, x_grow, NX, ixl_p, ixr_p, varphi_x)
      call linint_Grow(k_p, k_l, k_u, k_grow, NK-1, ikl_p, ikr_p, varphi_k)

      ! restrict values to grid just in case
      ial_p = min(ial_p, NA)
      iar_p = min(iar_p, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)

      ixl_p = min(ixl_p, NX)
      ixr_p = min(ixr_p, NX)
      varphi_x = max(min(varphi_x, 1d0),0d0)

      ikl_p = min(ikl_p+1, NK)
      ikr_p = min(ikr_p+1, NK)
      varphi_k = max(min(varphi_k, 1d0), 0d0)

      S_temp = (1d0-psi(ij_com+1))*mu_b*max((1d0-omega_x)*Q(iq_p_com), 1d-16)**egam/egam

      ! get optimal investment strategy
      EV_temp = varphi_a*varphi_x*varphi_k                  *(egam*EV(ij_com+1, ial_p, ixl_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                varphi_a*varphi_x*(1d0-varhpi_k)            *(egam*EV(ij_com+1, ial_p, ixl_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam) + &
                varphi_a*(1d0-varphi_x)*varhpi_k            *(egam*EV(ij_com+1, ial_p, ixr_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                varphi_a*(1d0-varphi_x)*(1d0-varhpi_k)      *(egam*EV(ij_com+1, ial_p, ixr_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam) + &
                (1d0-varphi_a)*varphi_x*varphi_k            *(egam*EV(ij_com+1, iar_p, ixl_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                (1d0-varphi_a)*varphi_x*(1d0-varphi_k)      *(egam*EV(ij_com+1, iar_p, ixl_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam) + &
                (1d0-varphi_a)*(1d0-varphi_x)*varphi_k      *(egam*EV(ij_com+1, iar_p, ixr_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                (1d0-varphi_a)*(1d0-varphi_x)*(1d0-varphi_k)*(egam*EV(ij_com+1, iar_p, ixr_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam)

      S_temp = S_temp + psi(ij_com+1)*EV_temp**egam/egam

      inv_e = - (S_temp + 1d-16**egam/egam*abs(a_p-a_temp))

  end function

    ! the first order condition regarding consumption
    function cons_o(x_in)

        implicit none

        ! input variable
        real*8, intent(in) :: x_in(:)

        ! variable declarations
        real*8 :: cons_o, Q_plus, ind_o, income, tomorrow, varphi_q, varphi_p
        integer :: iql_p, iqr_p, ipl_p, ipr_p

        ! calculate tomorrow's assets
        Q_plus  = x_in(1)
        if (ij_com < JR) then
          lab_com = max(x_in(2), 0d0)
        else
          lab_com = 0d0
        endif

        ! current occupation
        ind_o = abs(dble(ik_com > 0))

        ! calculate current income
        income = (1d0-ind_o)*w*eff(ij_com)*eta(iw_com)*lab_com + &
                 ind_o*theta(ie_com)*(k(ik_com)**alpha*(eff(ij_com)*lab_com)**(1d0-alpha))**nu + (1d0-delta_k)*k(ik_com)

        cons_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + income + pen(ij_com, ip_com) - (1d0-(1d0-phi)*ind_o)*taup*min(income, p_u)  - Q_plus

        if (ij_com >= JR) then
          p_plus_com = p(ip_com)
        else
          p_plus_com = (p(ip_com)*dble(ij_com-1) + (1d0-(1d0-phi)*ind_o)*mu*(lambda + (1d0-lambda)*min(income, p_u)))/dble(ij_com)
        endif

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(Q_plus, Q_l, Q_u, Q_grow, NQ, iql_p, iqr_p, varphi_q)
        call linint_Equi(p_plus_com, p_l, p_u, NP, ipl_p, ipr_p, varphi_p)

        ! restrict values to grid just in case
        iql_p = min(iql_p, NQ)
        iqr_p = min(iqr_p, NQ)
        varphi_q = max(min(varphi_q, 1d0),0d0)

        ! restrict values to grid just in case
        ipl_p = min(ipl_p, NP)
        ipr_p = min(ipr_p, NP)
        varphi_p = max(min(varphi_p, 1d0),0d0)

        ! get next period value function
        tomorrow = 0d0
        if (ij_com < JJ .or. mu_b /= 0d0) then

          if(varphi_q <= varphi_p) then
            tomorrow = max(varphi_q           *(egam*S(ij_com, iql_p, ix_com, ipl_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                           (varphi_p-varphi_q)*(egam*S(ij_com, iqr_p, ix_com, ipl_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                           (1d0-varphi_p)     *(egam*S(ij_com, iqr_p, ix_com, ipr_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam), 1d-10)**egam/egam
          else
            tomorrow = max(varphi_p           *(egam*S(ij_com, iql_p, ix_com, ipl_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                           (varphi_q-varphi_p)*(egam*S(ij_com, iql_p, ix_com, ipr_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                           (1d0-varphi_q)     *(egam*S(ij_com, iqr_p, ix_com, ipr_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam), 1d-10)**egam/egam
           endif

        endif

        if(cons_com <= 0d0)then
           cons_o = -1d-16**egam/egam*(1d0+abs(cons_com))
        elseif(lab_com < 0d0) then
          cons_o = -1d-16**egam/egam*(1d0+abs(lab_com))
        elseif(lab_com >= 1d0) then
          cons_o = -1d-16**egam/egam*lab_com
        else
           cons_o = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + beta*tomorrow)
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
    ! Checks for the maximum gridpoint used
    !##############################################################################
    subroutine check_grid(iamax, ixmax, ikmax)

      implicit none

      !##### INPUT/OUTPUT VARIABLES #############################################
      integer :: iamax(JJ), ixmax(JJ), ikmax(JJ)

      !##### OTHER VARIABLES ####################################################
      integer :: ij, ia, ix, ik

      iamax = 0
      ixmax = 0
      ikmax = 0

      do ij = 1, JJ

        ! check for the maximum asset grid point used at a certain age
        do ia = NA, 0, -1
          if (sum(m(ij, ia, :, :, :, :, :)) > 0d0) then
            iamax(ij) = ia
            exit
          endif
        enddo ! ia

        ! check for the maximum annuitie grid point used at a certain age
        do ix = NX, 0, -1
          if (sum(m(ij, :, ix, :, :, :, :)) > 0d0) then
            ixmax(ij) = ix
            exit
          endif
        enddo ! ix

        ! check for the maximum annuitie grid point used at a certain age
        do ik = NK, 0, -1
          if (sum(m(ij, :, :, :, ik, :, :)) > 0d0) then
            ikmax(ij) = ik
            exit
          endif
        enddo ! ik

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
