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
    integer, parameter :: NX = 32

    ! number of points on the liquid asset grid
    integer, parameter :: NA = 32

    ! number of points on the capital grid
    integer, parameter :: NK = 32

    ! number of points on the pension claim grid
    integer, parameter :: NP = 5

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
    real*8, parameter :: k_min = 0.2d0
    real*8, parameter :: phi_k = 0.2d0
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: nu = 0.88d0

    ! size of the asset grid
    real*8, parameter :: X_l    = 0d0
    real*8, parameter :: X_u    = 50d0
    real*8, parameter :: X_grow = 0.075d0

    ! size of the liquid asset grid
    real*8, parameter :: a_l    = X_l
    real*8, parameter :: a_u    = X_u
    real*8, parameter :: a_grow = X_grow

    ! size of the pension claim grid
    real*8, parameter :: p_l    = 0d0
    real*8, parameter :: p_u    = 2d0

    ! size of the capital grid
    real*8, parameter :: k_l = 0d0
    real*8, parameter :: k_u = X_u/(1d0-xi)
    real*8, parameter :: k_grow = X_grow

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.5d0

    ! measure time
    integer :: time

    ! discretized shocks
    real*8 :: dist_eta(NW), pi_eta(NW, NW), eta(NW), dist_theta(NE), pi_theta(NE, NE), theta(NE)

    ! wages, transfer payments (old-age), survival probabilities and discount factor for housing utilty
    real*8 :: w, eff(JJ), pen(JJ+1, 0:NP), psi(JJ+1)

    ! government variables
    real*8 :: lambda, phi, mu
    real*8 :: taup

    ! cohort aggregate variables
    real*8 :: c_coh(JJ, 0:1), y_coh(JJ, 0:1), l_coh(JJ, 0:1), o_coh(JJ)
    real*8 :: a_coh(JJ, 0:1), k_coh(JJ)

    ! different grids to discretize the state space
    real*8 :: X(0:NX), a(0:NA), p(0:NP), k(0:NK)

    ! variables to store the policy functions
    real*8 :: X_plus(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE), a_plus(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE)
    real*8 :: p_plus(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE), k_plus(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE)
    real*8 :: c(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE), l(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE)

    ! variables for temporary policy and value functions
    real*8 :: X_plus_t(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE, 0:NO), a_plus_t(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE, 0:NO)
    real*8 :: k_plus_t(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE, 0:NO), p_plus_t(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE, 0:NO)
    real*8 :: c_t(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE, 0:NO), l_t(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE, 0:NO)
    real*8 :: V_t(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE, 0:NO)

    ! variables to store the portfolio choice decisions
    real*8 :: omega_k(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE)

    ! variables to store the value functions
    real*8 :: V(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE), EV(JJ+1, 0:NA, 0:NP, 0:NK, NW, NE), S(JJ+1, 0:NX, 0:NP, 0:NK, NW, NE, 0:NO)

    ! weights for the different gridpoints on the discretized state space
    real*8 :: m(JJ, 0:NA, 0:NP, 0:NK, NW, NE)

    ! numerical variables
    integer :: ij_com, ix_com, ia_com, ip_com, ik_com, iw_com, ie_com, ia_p_com, ip_p_com, ix_p_com, io_p_com
    real*8 :: cons_com, lab_com, p_plus_com

    !$omp threadprivate(ij_com, ix_com, ia_com, ip_com, ik_com, iw_com, ie_com, ia_p_com, ip_p_com, ix_p_com, io_p_com)
    !$omp threadprivate(cons_com, lab_com, p_plus_com)


  contains

    ! solve the household's decision of how much wealth to invest into capital
    subroutine solve_worker(ij, ix_p, ip_p, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, ix_p, ip_p, ik, iw, ie
        integer :: ial, iar
        real*8 :: a_plus, EV_temp, S_temp, varphi_a

        a_plus  = X(ix_p)

       ! calculate linear interpolation for future assets
       call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

       ! restrict values to grid just in case
       ial = min(ial, NA)
       iar = min(iar, NA)
       varphi_a = max(min(varphi_a, 1d0),0d0)

       S_temp = (1d0-psi(ij+1))*mu_b*max(a_plus, 1d-16)**egam/egam

       EV_temp = varphi_a      *(egam*EV(ij+1, ial, ip_p, 0, iw, ie))**(1d0/egam) + &
                 (1d0-varphi_a)*(egam*EV(ij+1, iar, ip_p, 0, iw, ie))**(1d0/egam)

       S_temp = S_temp + psi(ij+1)*EV_temp**egam/egam

       S(ij, ix_p, ip_p, ik, iw, ie, 0) = S_temp

    end subroutine

    ! solve the household's decision of how much wealth to invest into firm capital
    subroutine solve_entrepreneur(ij, ix_p, ip_p, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, ix_p, ip_p, ik, iw, ie
        real*8 :: x_in, fret

        ! set up communication variables
        ij_com = ij; ix_p_com = ix_p; ip_p_com = ip_p; ik_com = ik; iw_com = iw; ie_com = ie

        if (X(ix_p) > (1d0-xi)*k_min + tr(k(ik), k_min)) then

           ! get best guess for the root of foc_real
           if(ix_p > 0)then
              x_in = omega_k(ij+1, ix_p, ip_p, ik, iw, ie)
           else
              x_in = 1d-4
           endif

           ! solve the household problem using fminsearch
           call fminsearch(x_in, fret, 0d0, 1d0, inv_o)

           ! portfolio share for capital
           omega_k(ij, ix_p, ip_p, ik, iw, ie) = x_in
           S(ij, ix_p, ip_p, ik, iw, ie, 1) = -fret

        else

          omega_k(ij, ix_p, ip_p, ik, iw, ie) = 1d0
          S(ij, ix_p, ip_p, ik, iw, ie, 1) = 1d-16**egam/egam

        endif

    end subroutine


  ! solve the household's consumption-savings decision
  subroutine solve_consumption(ij, ia, ip, ik, iw, ie, io_p)

      implicit none

      integer, intent(in) :: ij, ia, ip, ik, iw, ie, io_p
      real*8 :: x_in(2), fret, varphi_x, varphi_p, k_p
      integer :: ixl_p, ixr_p, ipl_p, ipr_p

      ! set up communication variables
      ij_com = ij; ia_com = ia; ip_com = ip; ik_com = ik; iw_com = iw; ie_com = ie; io_p_com = io_p

      ! get best initial guess from future period
      x_in(1) = max(X_plus_t(ij+1, ia, ip, ik, iw, ie, io_p), 1d-4)
      x_in(2) = max(l_t(ij+1, ia, ip, ik, iw, ie, io_p), 0.33d0)

      ! solve the household problem using rootfinding
      call fminsearch(x_in, fret, (/X_l, 0d0/), (/X_u, 0.99d0/), value_func)

      ! determine future investment
      k_p = 0d0

      if (io_p == 1) then

        call linint_Grow(x_in(1), x_l, x_u, x_grow, NX, ixl_p, ixr_p, varphi_x)
        call linint_Equi(p_plus_com, p_l, p_u, NP, ipl_p, ipr_p, varphi_p)

        ! restrict values to grid just in case
        ixl_p = min(ixl_p, NX)
        ixr_p = min(ixr_p, NX)
        varphi_x = max(min(varphi_x, 1d0),0d0)

        ipl_p = min(ipl_p, NP)
        ipr_p = min(ipr_p, NP)
        varphi_p = max(min(varphi_p, 1d0),0d0)

        ! get next period's capital size
        k_p = ((1d0-xi)*k_min + (varphi_x*varphi_p            *omega_k(ij, ixl_p, ipl_p, ik, iw, ie) +  &
                                 varphi_x*(1d0-varphi_p)      *omega_k(ij, ixl_p, ipr_p, ik, iw, ie) +  &
                                 (1d0-varphi_x)*varphi_p      *omega_k(ij, ixr_p, ipl_p, ik, iw, ie) +  &
                                 (1d0-varphi_x)*(1d0-varphi_p)*omega_k(ij, ixr_p, ipr_p, ik, iw, ie))*(x_in(1)-(1d0-xi)*k_min))/(1d0-xi)

      endif

      ! copy decisions
      X_plus_t(ij, ia, ip, ik, iw, ie, io_p) = x_in(1)
      a_plus_t(ij, ia, ip, ik, iw, ie, io_p) = x_in(1) - (1d0-xi)*k_p
      k_plus_t(ij, ia, ip, ik, iw, ie, io_p) = k_p
      c_t(ij, ia, ip, ik, iw, ie, io_p) = cons_com
      l_t(ij, ia, ip, ik, iw, ie, io_p) = x_in(2)
      V_t(ij, ia, ip, ik, iw, ie, io_p) = -fret

  end subroutine

  ! the first order condition with respect to next period real estate
  function inv_o(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: inv_o, a_p, k_p, EV_temp, S_temp, omega_k, varphi_k, varphi_a, a_temp
      integer :: ikl_p, ikr_p, ial_p, iar_p

      ! store real estate share
      omega_k  = x_in

      ! determine future liquid wealth and future downpayment
      k_p =((1d0-xi)*k_min + omega_k*(X(ix_p_com)-(1d0-xi)*k_min))/(1d0-xi)
      a_temp = X(ix_p_com) - (1d0-xi)*k_p - tr(k(ik_com), k_p)
      a_p = max(a_temp, 0d0)

      ! derive interpolation weights
      call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial_p, iar_p, varphi_a)
      call linint_Grow(k_p, k_l, k_u, k_grow, NK-1, ikl_p, ikr_p, varphi_k)

      ! restrict values to grid just in case
      ial_p = min(ial_p, NA)
      iar_p = min(iar_p, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)

      ! restrict values to grid just in case
      ikl_p = min(ikl_p+1, NK)
      ikr_p = min(ikr_p+1, NK)
      varphi_k = max(min(varphi_k, 1d0), 0d0)

      S_temp = (1d0-psi(ij_com+1))*mu_b*max(X(ix_p_com), 1d-16)**egam/egam

      ! get optimal investment strategy
      if(varphi_a <= varphi_k)then
          EV_temp = varphi_a           *(egam*EV(ij_com+1, ial_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                    (varphi_k-varphi_a)*(egam*EV(ij_com+1, iar_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                    (1d0-varphi_k)     *(egam*EV(ij_com+1, iar_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam)
      else
          EV_temp = varphi_k           *(egam*EV(ij_com+1, ial_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                    (varphi_a-varphi_k)*(egam*EV(ij_com+1, ial_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam) + &
                    (1d0-varphi_a)     *(egam*EV(ij_com+1, iar_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam)
      endif

      S_temp = S_temp + psi(ij_com+1)*EV_temp**egam/egam

      inv_o = - (S_temp + 1d-16**egam/egam*abs(a_p-a_temp))

  end function

    ! the first order condition regarding consumption
    function value_func(x_in)

        implicit none

        ! input variable
        real*8, intent(in) :: x_in(:)

        ! variable declarations
        real*8 :: value_func, X_plus, ind_o, income, tomorrow, varphi_x, varphi_p
        integer :: ixl_p, ixr_p, ipl_p, ipr_p

        ! calculate tomorrow's assets
        X_plus  = x_in(1)
        lab_com = max(x_in(2), 0d0)

        ! current occupation
        ind_o = abs(dble(ik_com > 0))

        income = (1d0-ind_o)*w*eff(ij_com)*eta(iw_com)*lab_com + &
                 ind_o*theta(ie_com)*(k(ik_com)**alpha*(eff(ij_com)*lab_com)**(1d0-alpha))**nu + (1d0-delta_k)*k(ik_com)

        cons_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + income + pen(ij_com, ip_com) - X_plus

        if (ij_com >= JR) then
          p_plus_com = p(ip_com)
        else
          p_plus_com = (p(ip_com)*dble(ij_com-1) + (1d0-(1d0-phi)*ind_o)*mu*(lambda + (1d0-lambda)*min(income, p_u)))/dble(ij_com)
        endif

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(X_plus, X_l, X_u, X_grow, NX, ixl_p, ixr_p, varphi_x)
        call linint_Equi(p_plus_com, p_l, p_u, NP, ipl_p, ipr_p, varphi_p)

        ! restrict values to grid just in case
        ixl_p = min(ixl_p, NX)
        ixr_p = min(ixr_p, NX)
        varphi_x = max(min(varphi_x, 1d0),0d0)

        ! restrict values to grid just in case
        ipl_p = min(ipl_p, NP)
        ipr_p = min(ipr_p, NP)
        varphi_p = max(min(varphi_p, 1d0),0d0)

        if (ij_com == JJ .and. mu_b == 0d0) then
          tomorrow = 0d0
        else
          ! get next period value function
          tomorrow = max(varphi_x*varphi_p              *(egam*S(ij_com, ixl_p, ipl_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                         varphi_x*(1d0-varphi_p)        *(egam*S(ij_com, ixl_p, ipr_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                         (1d0-varphi_x)*varphi_p        *(egam*S(ij_com, ixr_p, ipl_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                         (1d0-varphi_x)*(1d0-varphi_p)  *(egam*S(ij_com, ixr_p, ipr_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam), 1d-10)**egam/egam
        endif

        if(cons_com <= 0d0)then
           value_func = -1d-16**egam/egam*(1d0+abs(cons_com))
        elseif(lab_com < 0d0) then
          value_func = -1d-16**egam/egam*(1d0+abs(lab_com))
        elseif(lab_com >= 1d0) then
          value_func = -1d-16**egam/egam*lab_com
        else
           value_func = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + beta*tomorrow)
        endif

    end function

    ! function that computes adjustment cost
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

end module
