!##############################################################################
! MODULE globals
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################

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

    ! number of occupation states
    integer, parameter :: NO = 1

    ! household preference parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: beta = 0.99d0**5d0
    real*8, parameter :: mu_b = 0.50d0
    ! convert variables into per period values

    ! risk free rate and risk premium
    real*8, parameter :: r  = 0.5d0

    ! capital parameters
    real*8, parameter :: delta_k = 0.06d0
    real*8, parameter :: xi = 0d0 !1d0/3d0

    ! production parameters
    real*8, parameter :: k_min = 0.4d0
    real*8, parameter :: phi_k = 0.2d0
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: nu = 0.88d0

    ! risk processes
    real*8, parameter :: sigma_eta  = 0.0738d0
    real*8, parameter :: sigma_theta   = 0.0106d0
    real*8, parameter :: rho         = 0d0

    ! size of the asset grid
    real*8, parameter :: X_l    = 0d0
    real*8, parameter :: X_u    = 50d0
    real*8, parameter :: X_grow = 0.075d0

    ! size of the liquid asset grid
    real*8, parameter :: a_l    = X_l
    real*8, parameter :: a_u    = X_u
    real*8, parameter :: a_grow = X_grow

    ! size of the capital grid
    real*8, parameter :: k_l = 0d0
    real*8, parameter :: k_u = X_u/(1d0-xi)
    real*8, parameter :: k_grow = X_grow

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.5d0

    ! discretized shocks
    real*8 :: dist_eta(NW), pi_eta(NW, NW), eta(NW), dist_theta(NE), pi_theta(NE, NE), theta(NE)

    ! wages, transfer payments (old-age), survival probabilities and discount factor for housing utilty
    real*8 :: w, eff(JJ), pen(JJ+1), psi(JJ+1)

    ! cohort aggregate variables
    real*8 :: c_coh(JJ, 0:1), y_coh(JJ, 0:1), o_coh(JJ)
    real*8 :: a_coh(JJ, 0:1), k_coh(JJ)

    ! different grids to discretize the state space
    real*8 :: X(0:NX), a(0:NA), k(0:NK)

    ! variables to store the policy functions
    real*8 :: X_plus(JJ+1, 0:NA, 0:NK, NW, NE), a_plus(JJ+1, 0:NA, 0:NK, NW, NE), k_plus(JJ+1, 0:NA, 0:NK, NW, NE)
    real*8 :: c(JJ+1, 0:NA, 0:NK, NW, NE)

    ! variables for temporary policy and value functions
    real*8 :: X_plus_t(JJ+1, 0:NA, 0:NK, NW, NE, 0:NO), a_plus_t(JJ+1, 0:NA, 0:NK, NW, NE, 0:NO), k_plus_t(JJ+1, 0:NA, 0:NK, NW, NE, 0:NO)
    real*8 :: c_t(JJ+1, 0:NA, 0:NK, NW, NE, 0:NO)
    real*8 :: V_t(JJ+1, 0:NA, 0:NK, NW, NE, 0:NO)

    ! variables to store the portfolio choice decisions
    real*8 :: omega_k(JJ+1, 0:NA, 0:NK, NW, NE)

    ! variables to store the value functions
    real*8 :: V(JJ+1, 0:NA, 0:NK, NW, NE), EV(JJ+1, 0:NA, 0:NK, NW, NE), S(JJ+1, 0:NA, 0:NK, NW, NE, 0:NO)

    ! weights for the different gridpoints on the discretized state space
    real*8 :: m(JJ, 0:NA, 0:NK, NW, NE)

    ! numerical variables
    integer :: ij_com, ix_com, ia_com, ik_com, ia_p_com, ix_p_com, iw_com, ie_com
    real*8 :: cons_com

  contains

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
