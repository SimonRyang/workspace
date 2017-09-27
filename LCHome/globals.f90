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
    integer, parameter :: JR = 45

    ! number of years the household lives
    integer, parameter :: JJ = 80

    ! number of white noise (zeta) shocks
    integer, parameter :: NW = 5   !5

    ! number of transitory (epsilon) shocks
    integer, parameter :: NS = 5   !5

    ! number of rate of return (vtheta) shocks
    integer, parameter :: NR = 5   !5

    ! number of eps-vtheta shocks
    integer, parameter :: NSR = NS*NR

    ! number of eta shocks
    integer, parameter :: NE = 1000

    ! number of points on the cash-on-hand-grid
    integer, parameter :: NX = 50

    ! number of points on the asset grid
    integer, parameter :: NA = 50

    ! number of points on the liquid asset grid
    integer, parameter :: NL = 50

    ! number of points on the downpayment grid
    integer, parameter :: ND = 50

    ! number of points on the real estate grid
    integer, parameter :: NH = 50

    ! number of ownership states
    integer, parameter :: NO = 1

    ! number of LTC states
    integer, parameter :: NC = 1

    ! household preference parameters
    real*8, parameter :: gamma = 0.2d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: beta = 0.96d0
    real*8, parameter :: theta = 0.6d0
    real*8, parameter :: chi = 1.0d0
    real*8, parameter :: nu = 0d0

    ! risk free rate and risk premium
    real*8, parameter :: r_f  = 0.02d0
    real*8, parameter :: mu_r = 0.0d0
    real*8, parameter :: rp = 0.0d0

    ! housing parameters
    real*8, parameter :: delta_h = 0.015d0
!    real*8, parameter :: h_min = 4d0
    real*8, parameter :: h_min = 0.0001d0
    real*8, parameter :: xi = 0.7d0
    real*8, parameter :: phi_1 = 0.02d0
    real*8, parameter :: phi_2 = 0.06d0
!    real*8, parameter :: phi_1 = 0.0d0
!    real*8, parameter :: phi_2 = 0.0d0

    ! risk processes
    real*8, parameter :: sigma_zeta  = 0.0738d0
    real*8, parameter :: sigma_eps   = 0.0106d0
    real*8, parameter :: sigma_vtheta= 0.157d0**2d0
    real*8, parameter :: rho         = 0d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 500d0
    real*8, parameter :: a_grow = 0.075d0

    ! size of the liquid asset grid
    real*8, parameter :: l_l    = a_l
    real*8, parameter :: l_u    = a_u
    real*8, parameter :: l_grow = a_grow

    ! size of the downpayment grid
    real*8, parameter :: d_l    = (1d0-xi)*h_min
    real*8, parameter :: d_u    = a_u*0.5d0
    real*8, parameter :: d_grow = a_grow

    ! size of the real estate grid
    real*8, parameter :: h_l = h_min
    real*8, parameter :: h_u = d_u/(1d0-xi)
    real*8, parameter :: h_grow = a_grow

    ! growth of the cash-on-hand grid
    real*8, parameter :: X_grow = 0.175d0

    ! pension fraction of last income
    real*8, parameter :: kappa_1 = 0.5d0

    ! LTC parameters
    real*8, parameter :: p_l    = 0.01d0
    real*8, parameter :: p_u    = 0.60d0
    real*8, parameter :: p_grow = 0.05d0
    real*8, parameter :: kappa_2 = 0.25d0

    ! should cohort averages and variance be calculated analytically
    logical, parameter :: analytical = .true.

    ! lower and upper bound for X-grid
    real*8 :: X_l, X_u

    ! lower and upper bound for the eta-grid
    real*8 :: eta_l(JJ), eta_u(JJ)

    ! discretized shocks
    real*8 :: dist_zeta(NW), zeta(NW), dist_ltc(JJ+2, 0:NC, 0:NC)
    real*8 :: dist_epsvtheta(NSR), eps(NSR), vtheta(NSR)

    ! wages, transfer payments (old-age), survival probabilities and discount factor for housing utilty
    real*8 :: w, ph, eff(JJ), pen(JJ+1), ltc(JJ+1), psi(JJ+1, 0:NC), rpop(0:JJ+1)

    ! cohort aggregate variables
    real*8 :: c_coh(JJ, 0:2), ch_coh(JJ, 0:2), y_coh(JJ), o_coh(JJ), omega_coh(JJ, 0:2)
    real*8 :: a_coh(JJ, 0:2), al_coh(JJ, 0:2), h_coh(JJ) 

    ! different grids to discretize the state space
    real*8 :: a(0:NA), d(0:ND), h(0:NH), l(0:NL), X(0:NX)

    ! variables to store the policy functions
    real*8 :: a_plus(JJ+1, 0:NX, 0:NC, 0:NH), h_plus(JJ+1, 0:NX, 0:NC, 0:NH)
    real*8 :: c(JJ+1, 0:NX, 0:NC, 0:NH), ch(JJ+1, 0:NX, 0:NC, 0:NH)

    ! variables for temporary policy and value functions
    real*8 :: a_plus_t(JJ+1, 0:NX, 0:NC, 0:NH, 0:NO), h_plus_t(JJ+1, 0:NX, 0:NC, 0:NH, 0:NO)
    real*8 :: ch_t(JJ+1, 0:NX, 0:NC, 0:NH, 0:NO), c_t(JJ, 0:NX, 0:NC, 0:NH, 0:NO)
    real*8 :: V_t(JJ, 0:NX, 0:NC, 0:NH, 0:NO)

    ! variables to store the portfolio choice decisions
    real*8 :: omega_h(JJ+1, 0:NC, 0:NA, 0:NH), omega_plus(JJ+1, 0:NC, 0:NL, 0:ND)

    ! variables to store the value functions
    real*8 :: V(JJ+1, 0:NX, 0:NC, 0:NH), Q(JJ, 0:NC, 0:NL, 0:ND), S(0:JJ, 0:NC, 0:NA, 0:NH, 0:NO) 

    ! variables for the distribution of the eta shock
    real*8 :: eta(JJ, 0:NE), phi_e(JJ, 0:NE)

    ! weights for the different gridpoints on the discretized state space
    real*8 :: phi_ld(JJ, 0:NC, 0:NL, 0:ND), phi_X(JJ, 0:NX, 0:NC, 0:NH)

    ! numerical variables
    integer :: ij_com, ix_com, ia_com, ic_com, ih_com, ia_p_com, id_p_com, il_p_com
    real*8 :: cons_com, ch_com

contains


    ! function that computes adjustment cost
    function tr(h, h_p)

       implicit none

       ! input variable
       real*8, intent(in) :: h, h_p
       real*8 :: tr

       
       tr = 0d0
       if ((1-delta_h)*h > h_p) then 

          tr = phi_1*((1-delta_h)*h - h_p)

       elseif (h_p > h) then
        
          tr = phi_2*(h_p - h)

       endif
       
          
    end function

    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_a(iamax)

        implicit none

        integer :: iamax(JJ), ij, ix, ic, ih, ial, iar
        real*8 :: varphi_a 

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                do ic = 0, NC
                    do ih = 0, NH
                       if(phi_X(ij, ix, ic, ih) > 1d-8) then 
                          call linint_Grow(a_plus(ij, ix, ic, ih), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                          if (ij<JJ) iamax(ij+1) = max(min(iar, NA), iamax(ij+1))
                       endif   
                    enddo
                enddo
            enddo
        enddo

    end subroutine

    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_d(idmax)

        implicit none

        integer :: idmax(JJ), ij, ic, il, id

        idmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do id = 0, ND
                do ic = 0, NC
                    do il = 0, NL
                        if(phi_ld(ij, ic, il, id) > 1d-8)idmax(ij) = id
                    enddo
                enddo
            enddo
        enddo

    end subroutine

    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_l(ilmax)

        implicit none

        integer :: ilmax(JJ), ij, ic, il, id

        ilmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do il = 0, NL
                do ic = 0, NC
                    do id = 0, ND
                        if(phi_ld(ij, ic, il, id) > 1d-8)ilmax(ij) = il
                    enddo
                enddo
            enddo
        enddo

    end subroutine

    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_X(ixmax)

        implicit none

        integer :: ixmax(JJ), ij, ix, ic, ih

        ixmax = 0
        do ij = 1, JJ

            ! check for the maximum real estate grid point used at a certain age
            do ix = 0, NX
                do ic = 0, NC
                    do ih = 0, NH
                        if(phi_X(ij, ix, ic, ih) > 1d-8)ixmax(ij) = ix
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! generates the eta distribution
    subroutine generate_eta()

        implicit none

        integer :: ij, iel, ier, ie, isr
        real*8 :: varphi, eta_temp

        ! set bounds and grid for working ages
        eta_l(1)  = 0d0
        eta_u(1)  = 0d0
        eta(1, :) = 0d0

        do ij = 2, JR-1
            eta_l(ij)  = (ij-1)*(minval(log(eps)))
            eta_u(ij)  = (ij-1)*(maxval(log(eps)))
            call grid_Cons_Equi(eta(ij, :), eta_l(ij), eta_u(ij))
        enddo

        ! no innovations throughout retirement
        do ij = JR, JJ
            eta_l(ij)  = eta_l(JR-1)
            eta_u(ij)  = eta_u(JR-1)
            eta(ij, :) = eta(JR-1, :)
        enddo

        phi_e = 0d0

        ! initial distribution at age 1
        phi_e(1, :) = 1d0/dble(NE+1)

        ! iterate over working different years
        do ij = 2, JR-1

            ! last period's etas
            do ie = 0, NE

                ! new innovations
                do isr = 1, NSR

                    ! distribute on new grid
                    eta_temp = eta(ij-1, ie) + log(eps(isr))
                    call linint_Equi(eta_temp, eta_l(ij), eta_u(ij), NE, iel, ier, varphi)
                    phi_e(ij, iel) = phi_e(ij, iel) + dist_epsvtheta(isr)*varphi*phi_e(ij-1, ie)
                    phi_e(ij, ier) = phi_e(ij, ier) + dist_epsvtheta(isr)*(1d0-varphi)*phi_e(ij-1, ie)
                enddo
            enddo
        enddo

        ! no innovations throughout retirement
        do ij = JR, JJ
            phi_e(ij, :) = phi_e(JR-1, :)
        enddo

        ! take exponentials
        eta = exp(eta)

    end subroutine

    ! subroutine to calculate age specific quantiles of the distribution
    subroutine calculate_quantiles()

        use toolbox

        implicit none

        integer :: ij, ie, ia, ic, icount, io, it, NCOUNT
        real*8, allocatable :: a_sort(:), a_dist(:), a_cdist(:)
        integer, allocatable :: a_own(:), a_own_s(:), iorder(:)
        real*8 :: thresholds(5), quantiles(size(thresholds, 1), JJ), quantiles_own(size(thresholds, 1), JJ)
        real*8 :: own(size(thresholds, 1), JJ), slope, ages(JJ), thresholds_n(size(thresholds, 1)+1, JJ)
        real*8 :: sum, sum_2

        ! allocate arrays
        allocate(a_sort((NE+1)*(NA+1)*(NC+1)*(NO+1)), a_dist((NE+1)*(NA+1)*(NC+1)*(NO+1)), a_cdist((NE+1)*(NA+1)*(NC+1)*(NE+1)))
        allocate(iorder((NE+1)*(NA+1)*(NC+1)*(NO+1)), a_own((NE+1)*(NA+1)*(NC+1)*(NO+1)), a_own_s((NE+1)*(NA+1)*(NC+1)*(NO+1)))

        ! define quantile thresholds
        thresholds_n(:, :) = 0d0
        thresholds = (/0.05d0, 0.25d0, 0.50d0, 0.75d0, 0.95d0/)
        quantiles = 0d0
        quantiles_own = 0d0
        sum = 0d0
        sum_2 = 0d0

        ! iterate over ages
        do ij = 2, JJ

            a_sort = 0d0
            a_dist = 0d0

            ! copy savings into one-dimensional array
            icount = 1
            do ie = 0, NE
                do ic = 0, NC
                    do io = 0, NO
                        do ia = 0, NA
!                            if(phi_a(ij-1, ic, io, ia)*phi_e(ij, ie) > 1d-12)then
!                                a_sort(icount) = a(ia)*eta(ij, ie)
!                                a_dist(icount) = phi_a(ij-1, ic, io, ia)*phi_e(ij, ie)
!                                a_own(icount) = io
!                                icount = icount + 1
!                            endif
                        enddo
                    enddo
                enddo
            enddo
            NCOUNT = icount -1

            ! sort array and distribution
            call sort(a_sort(1:NCOUNT), iorder(1:NCOUNT))

            ! calculate cumulative distribution (attention ordering)
            a_cdist(1) = a_dist(iorder(1))
            do icount = 2, NCOUNT
                a_cdist(icount) = a_cdist(icount-1) + a_dist(iorder(icount))
            enddo

            ! get ownership in quantiles
            do it = 1, size(thresholds, 1)
                if(thresholds(it) <= a_cdist(1))then
                    quantiles_own(it, ij) = dble(a_own(iorder(1)))*a_dist(iorder(icount))
                else
                    do icount = 2, NCOUNT

                        sum = sum + dble(a_own(iorder(icount)))*a_dist(iorder(icount))

                        sum_2 = sum_2 + a_dist(iorder(icount))

                        if(thresholds(it) < a_cdist(icount))then
                            slope = dble(a_own(iorder(icount))-a_own(iorder(icount-1)))/(a_cdist(icount)-a_cdist(icount-1))

                            ! number of homeowners in the quantile
                            own(it, ij) = sum

                            ! compute upper threshold of quantile
                            thresholds_n(it, ij) = sum_2

                            ! compute homeownership rate in the quantile
                            quantiles_own(it, ij) = own(it, ij)/thresholds_n(it,ij)

                            exit
                        elseif(icount == NCOUNT)then
                            own(it, ij) = sum

                            ! compute upper threshold of quantile
                            thresholds_n(it, ij) = sum_2

                            ! compute homeownership rate in the quantile
                            quantiles_own(it, ij) = own(it, ij)/thresholds_n(it,ij)

                        endif
                    enddo
                endif

                ! reset sum and sum_2
                sum = 0d0
                sum_2 = 0d0

            enddo

            ! get quantiles
            do it = 1, size(thresholds, 1)
                if(thresholds(it) <= a_cdist(1))then
                    quantiles(it, ij) = a_sort(1)
                else
                    do icount = 2, NCOUNT
                        if(thresholds(it) < a_cdist(icount))then
                            slope = (a_sort(icount) - a_sort(icount-1))/(a_cdist(icount) - a_cdist(icount-1))
                            quantiles(it, ij) = a_sort(icount-1) + slope*(thresholds(it) - a_cdist(icount-1))
                            exit
                        elseif(icount == NCOUNT)then
                            quantiles(it, ij) = a_sort(NCOUNT)
                        endif
                    enddo
                endif
            enddo

        enddo

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)

        call plot(ages, quantiles_own(1, :), legend='5%')
        call plot(ages, quantiles_own(2, :), legend='25%')
        call plot(ages, quantiles_own(3, :), legend='50%')
        call plot(ages, quantiles_own(4, :), legend='75%')
        call plot(ages, quantiles_own(5, :), legend='95%')
        call execplot(title='Quantiles of Asset Distributions', &
            xlabel='Age j', ylabel='Homeownership')

        call plot(ages, quantiles(1, :), legend='50%')
        call plot(ages, quantiles(2, :), legend='25%')
        call plot(ages, quantiles(3, :), legend='50%')
        call plot(ages, quantiles(4, :), legend='75%')
        call plot(ages, quantiles(5, :), legend='95%')

        call execplot(title='Quantiles of Asset Distributions', &
            xlabel='Age j', ylabel='Assets')

    end subroutine

end module
