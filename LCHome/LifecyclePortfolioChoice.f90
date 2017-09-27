!##############################################################################
! PROGRAM Portfolio Choice with housing decision
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################
!
include "toolbox.f90"
include "globals.f90"
include "owner_foc.f90"
include "owner_solve.f90"
include "renter_foc.f90"
include "renter_solve.f90"

program PortfolioChoiceRealEstate

    ! load modules
    use globals
    use owner_solve
    use owner_foc
    use renter_foc
    use renter_solve

    implicit none

    ! initialize remaining variables
    call initialize()

    ! start the clock
    call tic()

    ! solve the household problem
    call solve_household()

    ! calculate the distribution of households over state space
    call get_distribution()

    ! aggregate individual decisions
    call aggregation()

    ! stop the clock
    call toc()

    call output()

    close(21)


contains


    ! initializes all remaining variables
    subroutine initialize

        implicit none
        real*8 :: temp(NSR, 2)


        ! wage rate for effective labor and rental price
        w = 1d0
        ph = r_f + delta_h

        ! set survival probabilities
        psi(:, 0) = (/1.00000d0, 0.99923d0, 0.99914d0, 0.99914d0, 0.99912d0, &
                      0.99906d0, 0.99908d0, 0.99906d0, 0.99907d0, 0.99901d0, &
                      0.99899d0, 0.99896d0, 0.99893d0, 0.99890d0, 0.99887d0, &
                      0.99886d0, 0.99878d0, 0.99871d0, 0.99862d0, 0.99853d0, &
                      0.99841d0, 0.99835d0, 0.99819d0, 0.99801d0, 0.99785d0, &
                      0.99757d0, 0.99735d0, 0.99701d0, 0.99676d0, 0.99650d0, &
                      0.99614d0, 0.99581d0, 0.99555d0, 0.99503d0, 0.99471d0, &
                      0.99435d0, 0.99393d0, 0.99343d0, 0.99294d0, 0.99237d0, &
                      0.99190d0, 0.99137d0, 0.99085d0, 0.99000d0, 0.98871d0, &
                      0.98871d0, 0.98721d0, 0.98612d0, 0.98462d0, 0.98376d0, &
                      0.98226d0, 0.98062d0, 0.97908d0, 0.97682d0, 0.97514d0, &
                      0.97250d0, 0.96925d0, 0.96710d0, 0.96330d0, 0.95965d0, &
                      0.95619d0, 0.95115d0, 0.94677d0, 0.93987d0, 0.93445d0, &
                      0.92717d0, 0.91872d0, 0.91006d0, 0.90036d0, 0.88744d0, &
                      0.87539d0, 0.85936d0, 0.84996d0, 0.82889d0, 0.81469d0, &
                      0.79705d0, 0.78081d0, 0.76174d0, 0.74195d0, 0.72155d0, &
                      0.00000d0/)

        ! lower survival probabilities if ltc shock
        psi(:, 1) = 0.877*psi(:, 0)

        ! initialize age earnings process
        eff(1:JR-1) = &
                     (/1.0000d0, 1.0719d0, 1.1438d0, 1.2158d0, 1.2842d0, 1.3527d0, &
                       1.4212d0, 1.4897d0, 1.5582d0, 1.6267d0, 1.6952d0, 1.7217d0, &
                       1.7438d0, 1.7748d0, 1.8014d0, 1.8279d0, 1.8545d0, 1.8810d0, &
                       1.9075d0, 1.9341d0, 1.9606d0, 1.9623d0, 1.9640d0, 1.9658d0, &
                       1.9675d0, 1.9692d0, 1.9709d0, 1.9726d0, 1.9743d0, 1.9760d0, &
                       1.9777d0, 1.9700d0, 1.9623d0, 1.9546d0, 1.9469d0, 1.9392d0, &
                       1.9315d0, 1.9238d0, 1.9161d0, 1.9084d0, 1.9007d0, 1.8354d0, &
                       1.7701d0, 1.7048d0/)

        ! earnings process is during retirement equal to zero
        eff(JR:JJ) = 0d0

        ! old-age transfers
        pen = 0d0
        ltc = 0d0
        pen(JR:JJ) = kappa_1*w*eff(JR-1)
        ltc(JR:JJ) = kappa_2*w*eff(JR-1)

        ! discretize zeta shocks
        call normal_discrete(zeta, dist_zeta, 0d0, sigma_zeta)
        zeta = exp(zeta)

        ! discretize eps-vtheta shocks
        call normal_discrete((/NS, NR/), temp, dist_epsvtheta, (/0d0, 0d0/), (/sigma_eps, sigma_vtheta/), rho)
        eps(:) = exp(temp(:, 1))
        vtheta(:)  = temp(:, 2)

        ! initialize ltc shock
        dist_ltc = 0d0

        ! during retirement ltc shocks arise
        call grid_Cons_Grow(dist_ltc(JR:JJ+2, 0, 1), p_l, p_u, p_grow)

        dist_ltc(JR:JJ+2, 0, 1) = 0d0
        dist_ltc(JR:JJ+2, 0, 0) = 1d0-dist_ltc(JR:JJ+2, 0, 1)

        ! ltc is a permanent shock
        dist_ltc(JR:JJ+2, 1, 1) = 1d0
        dist_ltc(JR:JJ+2, 1, 0) = 0d0


        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize liquid asset grid
        call grid_Cons_Grow(l, l_l, l_u, l_grow)

        ! initialize downpayment grid
        call grid_Cons_Grow(d(1:ND), d_l, d_u, d_grow)
        d(0) = 0d0

        ! endogenous upper bound of housing grid
        call grid_Cons_Grow(h(1:ND), h_l, h_u, h_grow)
        h(0) = 0d0

        ! endogenous lower and upper bound of cash-on-hand grid
        X_l = 0.0001d0
        X_u = (1d0 + r_f + mu_r + maxval(vtheta(:)))*a_u + w*maxval(eff(1:JR-1))*maxval(eps(:))*zeta(NW)
        call grid_Cons_Grow(X, X_l, X_u, X_grow)

        ! initialize value functions
        S = 1d-10**egam/egam; V = 1d-10**egam/egam; Q = 1d-10**egam/egam

        ! initialize policy functions
        a_plus = 0d0; h_plus = 0; c = 0d0; ch = 0d0
        omega_h = 0d0; omega_plus = 0d0

        ! initialize temporary policy and value functions
        a_plus_t = 0d0; h_plus_t = 0d0; ch_t = 0d0; c_t = 0d0; V_t = 0d0

        ! open files
        open(21, file='output.out')


    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none

        integer :: ij, ix, ic, ih, il_p, id_p, ia_p

        ! solve household problem recursively

        do ij = JJ, 1, -1

            if(ij == JJ)then

               omega_plus(JJ, :, :, :) = 0d0
               omega_h(JJ, :, :, :) = 0d0

               do il_p = 0, NL
                   do id_p = 1, ND
                       Q(JJ, :, il_p, id_p) = nu**(1d0/gamma)*(l(il_p)+d(id_p))**egam/egam
                   enddo
                   Q(JJ, :, il_p, 0) = Q(JJ, :, il_p, 1)
               enddo

               do ia_p = 1, NA
                   S(JJ, :, ia_p, :, :) = nu**(1d0/gamma)*a(ia_p)**egam/egam
               enddo
               S(JJ, :, 0, :, :) = 1d-10**egam/egam



               do ix = 0, NX
                   do ic = 0, NC
                       do ih = 0, NH

                           if (nu == 0d0) then

                               h_plus(JJ, ix, ic, ih) = 0d0
                               a_plus(JJ, ix, ic, ih) = 0d0

                               if (ih == 0) then
                                   ch_com = max((1d0-theta)*X(ix)/ph, 1d-10)
                                   cons_com = max(X(ix)-ph*ch_com, 1d-10)
                                   V(JJ, ix, ic, ih) = (cons_com**theta*(chi*ch_com)**(1d0-theta))**egam/egam
                                   c(JJ, ix, ic, ih) = cons_com
                                   ch(JJ, ix, ic, ih) = ch_com
                               else
                                   c(JJ, ix, ic, ih) = X(ix)
                                   ch(JJ, ix, ic, ih) = h(ih)

                                   if (ic == 0) then
                                       V(JJ, ix, ic, ih) = (max(X(ix), 1d-10)**theta*h(ih)**(1d0-theta))**egam/egam
                                   else
                                       V(JJ, ix, ic, ih) = (max(X(ix), 1d-10)**theta*(chi*h(ih))**(1d0-theta))**egam/egam
                                   endif
                               endif

                           else

                               ! with bequest motive we assume future renter
                               call solve_consumption_r(JJ, ix, ic, ih)

                               V(JJ, ix, ic, ih) = V_t(JJ, ix, ic, ih, 0)
                               h_plus(ij, ix, ic, ih) = h_plus_t(ij, ix, ic, ih, 0)
                               a_plus(ij, ix, ic, ih) = a_plus_t(ij, ix, ic, ih, 0)
                               c(ij, ix, ic, ih) = c_t(ij, ix, ic, ih, 0)
                               ch(ij, ix, ic, ih) = ch_t(ij, ix, ic, ih, 0)
                           endif
                       enddo
                   enddo
               enddo
            else


               ! solve portfolio choice problem
               do ic = 0, NC
                   do il_p = NL, 0, -1
                       do id_p = 0, ND

                           ! next period homeowner
                           if (id_p > 0) then

                               call solve_portfolio_o(ij, ic, il_p, id_p)

                           else

                               ! next period homeowner
                               call solve_portfolio_r(ij, ic, il_p, 0)
                           endif

                       enddo
                   enddo

               enddo

               ! get optimal share of wealth stored in real estate
               do ic = 0, NC
                   do ia_p = 0, NA
                       do ih = 0, NH

                           ! next period homeowner
                           call solve_realestate(ij, ic, ia_p, ih)

                           ! next period renter
                           call solve_renter(ij, ic, ia_p, ih)
                       enddo
                   enddo
               enddo


               ! solve the consumption savings problem
               do ix = 0, NX
                   do ic = 0, NC
                       do ih = 0, NH

                           ! next period homeowner
                           call solve_consumption_o(ij, ix, ic, ih)

                           ! next period renter
                           call solve_consumption_r(ij, ix, ic, ih)

                       enddo
                   enddo
               enddo

               ! decision whether to be owner or renter next period
               do ix = 0, NX
                   do ic = 0, NC
                       do ih = 0, NH

                           ! decision on whether to be homeowner or renter next period
                           if( V_t(ij, ix, ic, ih, 1) >= V_t(ij, ix, ic, ih, 0) ) then
                                 h_plus(ij, ix, ic, ih) = h_plus_t(ij, ix, ic, ih, 1)
                                 a_plus(ij, ix, ic, ih) = a_plus_t(ij, ix, ic, ih, 1)
                                 c(ij, ix, ic, ih) = c_t(ij, ix, ic, ih, 1)
                                 ch(ij, ix, ic, ih) = ch_t(ij, ix, ic, ih, 1)
                                 V(ij, ix, ic, ih) = V_t(ij, ix, ic, ih, 1)
                           else
                                 h_plus(ij, ix, ic, ih) = h_plus_t(ij, ix, ic, ih, 0)
                                 a_plus(ij, ix, ic, ih) = a_plus_t(ij, ix, ic, ih, 0)
                                 c(ij, ix, ic, ih) = c_t(ij, ix, ic, ih, 0)
                                 ch(ij, ix, ic, ih) = ch_t(ij, ix, ic, ih, 0)
                                 V(ij, ix, ic, ih) = V_t(ij, ix, ic, ih, 0)
                           endif
!                           if (h_plus(ij, ix, ic, ih) > 1d-5) write(*,'(4i5,f12.8)')ij, ix, ic, ih, h_plus(ij, ix, ic, ih)

                       enddo
                   enddo
               enddo

            endif
            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

        enddo

    end subroutine


    ! determines the invariant distribution over cash-on-hand space
    subroutine get_distribution()

        implicit none

        integer :: ij

        ! set distributions to zero
        phi_X = 0d0; phi_ld = 0d0

        do ij = 1, JJ

            ! get distribution on cash-on-hand grid
            call get_distribution_X(ij)

            ! get distribution after portolio choice decision
            call get_distribution_ld(ij)

            if (sum(phi_X(ij, :, :, :)) < 1d0-1d-10)write(*,'(a, i5, f8.2)')'X: ', ij, sum(phi_X(ij, :, :, :))
            if (sum(phi_ld(ij, :, :, :))< 1d0-1d-10)write(*,'(a, i5, f8.2)')'ld:', ij, sum(phi_ld(ij, :, :, :))


            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

        enddo

    end subroutine


    ! to calculate distribution on cash-on-hand grid
    subroutine get_distribution_X(ij)

        implicit none

        integer, intent(in) :: ij
        real*8 :: earnings, h_p, R_port, X_p, varphi_h, varphi_x, dist
        integer :: il, ic, ic_m, id, ihl, ihr, ixl, ixr, isr, iw

        if(ij == 1)then

            ! get initial distribution at age 1 of cash-on-hand
            do iw = 1, NW

                ! get initial cash-on-hand
                earnings  = w*eff(1)*zeta(iw)
                X_p = earnings

                ! derive interpolation weights
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                ! get distributional weight
                dist = dist_zeta(iw)

                ! initialize the distribution
                phi_X(1, ixl, 0, 0) = phi_X(1, ixl, 0, 0) + dist*varphi_X
                phi_X(1, ixr, 0, 0) = phi_X(1, ixr, 0, 0) + dist*(1d0-varphi_X)

            enddo

        elseif(ij <= JR-1)then

            ! iterate over yesterdays liquid asset and downpayment distribution
            do ic_m = 0, NC
                do isr = 1, NSR
                    do iw = 1, NW

                        ! derive labor earnings
                        earnings  = w*eff(ij)*zeta(iw)

                        ! get distributional weight
                        dist = dist_zeta(iw)*dist_epsvtheta(isr)

                        do il = 0, NL
                            do id = 0, ND

                                ! get today's cash-on-hand and interpolate
                                R_port = 1d0 + r_f + omega_plus(ij-1, ic_m, il, id)*(mu_r + vtheta(isr))
                                if(l(il) < xi*d(id)/(1d0-xi))R_port = 1d0 + r_f + rp

                                h_p = d(id)/(1d0-xi)/eps(isr)
                                X_p = R_port*(l(il)-xi*d(id)/(1d0-xi))/eps(isr) + earnings + (1d0 - delta_h)*h_p

                                ! derive interpolation weights
                                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                                if(id > 0)then
                                    call linint_Grow(h_p, h_l, h_u, h_grow, NH-1, ihl, ihr, varphi_h)
                                    ihl = min(ihl+1, NH)
                                    ihr = min(ihr+1, NH)
                                    varphi_h = max(min(varphi_h, 1d0), 0d0)
                                else
                                    ihl = 0 ; ihr = 0 ; varphi_h = 1d0
                                endif

                                ! restrict values to grid just in case
                                ixl = min(ixl, NX)
                                ixr = min(ixr, NX)
                                varphi_x = max(min(varphi_x, 1d0), 0d0)

                                ! distribute on today's state space
                                phi_X(ij, ixl, 0, ihl) = phi_X(ij, ixl, 0, ihl) + &
                                                         dist*varphi_X*varphi_h*phi_ld(ij-1, ic_m, il, id)
                                phi_X(ij, ixr, 0, ihl) = phi_X(ij, ixr, 0, ihl) + &
                                                   dist*(1d0-varphi_X)*varphi_h*phi_ld(ij-1, ic_m, il, id)
                                phi_X(ij, ixl, 0, ihr) = phi_X(ij, ixl, 0, ihr) + &
                                                   dist*varphi_X*(1d0-varphi_h)*phi_ld(ij-1, ic_m, il, id)
                                phi_X(ij, ixr, 0, ihr) = phi_X(ij, ixr, 0, ihr) + &
                                             dist*(1d0-varphi_X)*(1d0-varphi_h)*phi_ld(ij-1, ic_m, il, id)

                            enddo
                        enddo
                    enddo
                enddo
            enddo

        else

            ! iterate over yesterdays liquid asset and downpayment distribution
            do ic_m = 0, NC
                do isr = 1, NSR
                    do ic = 0, NC


                        ! get distributional weight
                        dist = dist_epsvtheta(isr)*dist_ltc(ij, ic_m, ic)

                        do il = 0, NL
                            do id = 0, ND

                                ! get today's cash-on-hand and interpolate
                                R_port = 1d0 + r_f + omega_plus(ij-1, ic_m, il, id)*(mu_r + vtheta(isr))
                                if(l(il) < xi*d(id)/(1d0-xi))R_port = 1d0 + r_f + rp

                                h_p = d(id)/(1d0-xi)
                                X_p = R_port*(l(il)-xi*d(id)/(1d0-xi)) + pen(ij) + (1d0 - delta_h)*h_p - dble(ic)*ltc(ij)
                                !if (X_p <= 0d0) stop

                                ! derive interpolation weights
                                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                                if(id > 0)then
                                    call linint_Grow(h_p, h_l, h_u, h_grow, NH-1, ihl, ihr, varphi_h)
                                    ihl = min(ihl+1, NH)
                                    ihr = min(ihr+1, NH)
                                    varphi_h = max(min(varphi_h, 1d0), 0d0)
                                else
                                    ihl = 0 ; ihr = 0; varphi_h = 1d0
                                endif

                                ! restrict values to grid just in case
                                ixl = min(ixl, NX)
                                ixr = min(ixr, NX)
                                varphi_x = max(min(varphi_x, 1d0), 0d0)

                                ! distribute on today's state space
                                phi_X(ij, ixl, ic, ihl) = phi_X(ij, ixl, ic, ihl) + &
                                                         dist*varphi_X*varphi_h*phi_ld(ij-1, ic_m, il, id)
                                phi_X(ij, ixr, ic, ihl) = phi_X(ij, ixr, ic, ihl) + &
                                                   dist*(1d0-varphi_X)*varphi_h*phi_ld(ij-1, ic_m, il, id)
                                phi_X(ij, ixl, ic, ihr) = phi_X(ij, ixl, ic, ihr) + &
                                                   dist*varphi_X*(1d0-varphi_h)*phi_ld(ij-1, ic_m, il, id)
                                phi_X(ij, ixr, ic, ihr) = phi_X(ij, ixr, ic, ihr) + &
                                             dist*(1d0-varphi_X)*(1d0-varphi_h)*phi_ld(ij-1, ic_m, il, id)

                            enddo
                        enddo
                    enddo
                enddo
            enddo

        endif

    end subroutine


    ! to calculate distribution after decision over wealth exposure in real estate
    subroutine get_distribution_ld(ij)

        implicit none

        integer, intent(in) :: ij
        real*8 :: al_p, ad_p, varphi_l, varphi_d
        integer :: ill, ilr, idl, idr, ix, ic, ih

        ! iterate over asset distribution
        do ix = 0, NX
            do ic = 0, NC
                do ih = 0, NH

                    ! determine future downpayment and future liquid wealth
                    ad_p = (1d0-xi)*h_plus(ij, ix, ic, ih)
                    al_p = a_plus(ij, ix, ic, ih)- ad_p - tr(h(ih), h_plus(ij, ix, ic, ih))

                    ! derive interpolation weights
                    call linint_Grow(al_p, l_l, l_u, l_grow, NL, ill, ilr, varphi_l)

                    if(ad_p > 1d-10)then
                        call linint_Grow(ad_p, d_l, d_u, d_grow, ND-1, idl, idr, varphi_d)
                        idl = min(idl+1, ND)
                        idr = min(idr+1, ND)
                        varphi_d = max(min(varphi_d, 1d0), 0d0)
                    else
                        idl = 0 ; idr = 0 ; varphi_d = 1d0
                    endif

                    ! restrict values to grid just in case
                    ill = min(ill, NL)
                    ilr = min(ilr, NL)
                    varphi_l = max(min(varphi_l, 1d0), 0d0)

                    ! get distribution over liquid asset and downpayment
                    phi_ld(ij, ic, ill, idl) = phi_ld(ij, ic, ill, idl) +   &
                                               varphi_l*varphi_d*phi_X(ij, ix, ic, ih)
                    phi_ld(ij, ic, ilr, idl) = phi_ld(ij, ic, ilr, idl) +   &
                                               (1d0-varphi_l)*varphi_d*phi_X(ij, ix, ic, ih)
                    phi_ld(ij, ic, ill, idr) = phi_ld(ij, ic, ill, idr) +   &
                                               varphi_l*(1d0-varphi_d)*phi_X(ij, ix, ic, ih)
                    phi_ld(ij, ic, ilr, idr) = phi_ld(ij, ic, ilr, idr) +   &
                                               (1d0-varphi_l)*(1d0-varphi_d)*phi_X(ij, ix, ic, ih)
                enddo
            enddo
        enddo

    end subroutine



    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none

        integer :: ij, ic, ic_m, ie, ih, iw, ix, id, il
        real*8 :: sigma_eta(JJ), mu_exp(JJ), sigma_exp(JJ)

        ! calculate cohort averages
        c_coh = 0d0; ch_coh = 0d0; y_coh = 0d0; o_coh = 0d0; a_coh = 0d0
        omega_coh = 0d0; h_coh = 0d0; al_coh = 0d0

        ! generate eta distribution if not analytical calculation
        if(.not. analytical)call generate_eta()

        ! compute cohort average wage income, consumption, assets,
        ! real estate and portfolio shares during working life

        ! analytical approach or not
        if(analytical)then

            do ij = 1, JJ


                ! wage income and pensions
                if(ij < JR)then
                    do iw = 1, NW
                        y_coh(ij) = y_coh(ij) + w*eff(ij)*zeta(iw)*dist_zeta(iw)
                    enddo
                else
                    y_coh(ij) = y_coh(ij) + pen(ij)
                endif

                ! consumption and rental consumption
                do ix = 0, NX
                    do ic = 0, NC
                        do ih = 0, NH
                            if (ih == 0) then
                               c_coh(ij, 0) = c_coh(ij, 0) + c(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)/ &
                                                   max(sum(phi_X(ij, :, :, 0)), 1d-10)
                               ch_coh(ij, 0) = ch_coh(ij, 0) + ch(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)/ &
                                                   max(sum(phi_X(ij, :, :, 0)), 1d-10)
                            else
                               c_coh(ij, 1) = c_coh(ij, 1) + c(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)/ &
                                                   max(sum(phi_X(ij, :, :, 1:NH)), 1d-10)
                               ch_coh(ij, 1) = ch_coh(ij, 1) + ch(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)/ &
                                                   max(sum(phi_X(ij, :, :, 1:NH)), 1d-10)
                            endif
                            c_coh(ij, 2) = c_coh(ij, 2) + c(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)
                            ch_coh(ij, 2) = ch_coh(ij, 2) + ch(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)

                            if (ij< JJ) then
                                if (h_plus(ij, ix, ic, ih) >= h_min) then
                                    o_coh(ij+1) = o_coh(ij+1) + phi_X(ij, ix, ic, ih)
                                    a_coh(ij+1, 1) = a_coh(ij+1, 1) + a_plus(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)
                                else
                                    a_coh(ij+1, 0) = a_coh(ij+1, 0) + a_plus(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)
                                endif
                                a_coh(ij+1, 2) = a_coh(ij+1, 2) + a_plus(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)
                                h_coh(ij+1) = h_coh(ij+1) + h_plus(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)
                            endif

                        enddo
                    enddo
                enddo
                if (ij < JJ) then
                    a_coh(ij+1, 1) = a_coh(ij+1, 1)/max(o_coh(ij+1), 1d-10)
                    a_coh(ij+1, 0) = a_coh(ij+1, 0)/max(1d0-o_coh(ij+1), 1d-10)
                endif

                ! liquid assets and housing
                if(ij>1)then

                    ! aggregated portfolio shares
                    do ic_m = 0, NC
                        do il = 0, NL
                            do id = 0, ND
                                if (id == 0) then
                                    omega_coh(ij, 0) = omega_coh(ij, 0) + omega_plus(ij-1, ic_m, il, 0)*phi_ld(ij-1, ic_m, il, 0)
                                    al_coh(ij, 0) = al_coh(ij, 0) + l(il)*phi_ld(ij-1, ic_m, il, 0)
                                else
                                    omega_coh(ij, 1) = omega_coh(ij, 1) + omega_plus(ij-1, ic_m, il, id)*phi_ld(ij-1, ic_m, il, id)
                                    al_coh(ij, 1) = al_coh(ij, 1) + l(il)*phi_ld(ij-1, ic_m, il, id)
                                endif
                                omega_coh(ij, 2) = omega_coh(ij, 2) + omega_plus(ij-1, ic_m, il, id)*phi_ld(ij-1, ic_m, il, id)
                                al_coh(ij, 2) = al_coh(ij, 2) + l(il)*phi_ld(ij-1, ic_m, il, id)
                            enddo
                        enddo
                    enddo
                    omega_coh(ij, 0) = omega_coh(ij, 0)/max(sum(phi_ld(ij-1, :, :, 0)), 1d-10)
                    omega_coh(ij, 1) = omega_coh(ij, 1)/max(sum(phi_ld(ij-1, :, :, 1:ND)), 1d-10)
                    al_coh(ij, 0) = al_coh(ij, 0)/max(sum(phi_ld(ij-1, :, :, 0)), 1d-10)
                    al_coh(ij, 1) = al_coh(ij, 1)/max(sum(phi_ld(ij-1, :, :, 1:ND)), 1d-10)
                endif


            enddo

            ! get age dependent variance of eta
            sigma_eta = sigma_eps*(/(dble(min(ij, JR-1)-1), ij=1,JJ)/)

            ! calculate age specific expectations and variance of exp(eta)
            mu_exp = exp(0.5d0*sigma_eta)
            sigma_exp = exp(sigma_eta)*(exp(sigma_eta)-1d0)

            ! add level effect to averages
            y_coh = mu_exp*y_coh
            c_coh(:,0) = mu_exp*c_coh(:,0)
            c_coh(:,1) = mu_exp*c_coh(:,1)
            c_coh(:,2) = mu_exp*c_coh(:,2)
            ch_coh(:, 0) = mu_exp*ch_coh(:, 0)
            ch_coh(:, 1) = mu_exp*ch_coh(:, 1)
            ch_coh(:, 2) = mu_exp*ch_coh(:, 2)
            a_coh(:, 0) = mu_exp*a_coh(:, 0)
            a_coh(:, 1) = mu_exp*a_coh(:, 1)
            a_coh(:, 2) = mu_exp*a_coh(:, 2)
            al_coh(:, 0) = mu_exp*al_coh(:, 0)
            al_coh(:, 1) = mu_exp*al_coh(:, 1)
            al_coh(:, 2) = mu_exp*al_coh(:, 2)
            h_coh = mu_exp*h_coh

        else

            do ij = 1, JJ

                do ie = 0, NE

                    ! wage income and pensions
                    if(ij < JR)then
                        do iw = 1, NW
                            y_coh(ij) = y_coh(ij) + w*eff(ij)*zeta(iw)*dist_zeta(iw)*eta(ij, ie)*phi_e(ij, ie)
                        enddo
                    else
                        y_coh(ij) = y_coh(ij) + pen(ij)*eta(ij, ie)*phi_e(ij, ie)
                    endif

                    ! consumption and rental consumption
                    do ix = 0, NX
                        do ic = 0, NC
                            do ih = 0, NH
                                if (ih == 0) then
                                    c_coh(ij, 0)  = c_coh(ij, 0) + c(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                    *eta(ij, ie)*phi_e(ij, ie)
                                    ch_coh(ij, 0) = ch_coh(ij, 0) + ch(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                    *eta(ij, ie)*phi_e(ij, ie)
                                else
                                    c_coh(ij, 1)  = c_coh(ij, 1) + c(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                    *eta(ij, ie)*phi_e(ij, ie)
                                    ch_coh(ij, 1) = ch_coh(ij, 1) + ch(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                    *eta(ij, ie)*phi_e(ij, ie)
                                endif
                                c_coh(ij, 2)  = c_coh(ij, 2) + c(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                    *eta(ij, ie)*phi_e(ij, ie)
                                ch_coh(ij, 2) = ch_coh(ij, 2) + ch(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                    *eta(ij, ie)*phi_e(ij, ie)

                                if (ij< JJ)o_coh(ij+1) = o_coh(ij+1) + h_plus(ij, ix, ic, ih)&
                                                             *phi_X(ij, ix, ic, ih)*phi_e(ij, ie)

                                if (ij< JJ) then
                                    if (h_plus(ij, ix, ic, ih) >= h_min) then
                                        o_coh(ij+1) = o_coh(ij+1) + phi_X(ij, ix, ic, ih)*phi_e(ij, ie)
                                        a_coh(ij+1, 1) = a_coh(ij+1, 1) + a_plus(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                            *eta(ij, ie)*phi_e(ij, ie)
                                    else
                                        a_coh(ij+1, 0) = a_coh(ij+1, 0) + a_plus(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                            *eta(ij, ie)*phi_e(ij, ie)
                                    endif
                                    a_coh(ij+1, 2) = a_coh(ij+1, 2) + a_plus(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                        *eta(ij, ie)*phi_e(ij, ie)
                                    h_coh(ij+1) = h_coh(ij+1) + h_plus(ij, ix, ic, ih)*phi_X(ij, ix, ic, ih)&
                                                                                        *eta(ij, ie)*phi_e(ij, ie)
                                endif
                            enddo
                        enddo
                    enddo

                    ! aggregated assets, liquid assets and housing
                    if(ij>1)then

                        do ic_m = 0, NC
                            do il = 0, NL
                                do id = 0, ND
                                    if (id == 0) then
                                        omega_coh(ij, 0) = omega_coh(ij, 0) + omega_plus(ij-1, ic_m, il, 0)*   &
                                                                              phi_ld(ij-1, ic_m, il, 0)*eta(ij, ie)*phi_e(ij, ie)
                                        al_coh(ij, 0) = al_coh(ij, 0) + l(il)*phi_ld(ij-1, ic_m, il, 0)*eta(ij, ie)*phi_e(ij, ie)
                                    else
                                        omega_coh(ij, 1) = omega_coh(ij, 1) + omega_plus(ij-1, ic_m, il, id)*   &
                                                                              phi_ld(ij-1, ic_m, il, id)*eta(ij, ie)*phi_e(ij, ie)
                                        al_coh(ij, 1) = al_coh(ij, 1) + l(il)*phi_ld(ij-1, ic_m, il, id)*eta(ij, ie)*phi_e(ij, ie)
                                    endif
                                    omega_coh(ij, 2) = omega_coh(ij, 2) + omega_plus(ij-1, ic_m, il, id)*   &
                                                                          phi_ld(ij-1, ic_m, il, id)*eta(ij, ie)*phi_e(ij, ie)
                                    al_coh(ij, 2) = al_coh(ij, 2) + l(il)*phi_ld(ij-1, ic_m, il, id)*eta(ij, ie)*phi_e(ij, ie)
                                enddo
                            enddo
                        enddo

                    endif

                enddo

                c_coh(ij, 0) = c_coh(ij, 0)/sum(phi_X(ij, :, :, 0))
                ch_coh(ij, 0) = ch_coh(ij, 0)/sum(phi_X(ij, :, :, 0))
                c_coh(ij, 1) = c_coh(ij, 1)/sum(phi_X(ij, :, :, 1:NH))
                ch_coh(ij, 1) = ch_coh(ij, 1)/sum(phi_X(ij, :, :, 1:NH))

                omega_coh(ij, 0) = omega_coh(ij, 0)/sum(phi_ld(ij-1, :, :, 0))
                omega_coh(ij, 1) = omega_coh(ij, 1)/sum(phi_ld(ij-1, :, :, 1:ND))
                al_coh(ij, 0) = al_coh(ij, 0)/sum(phi_ld(ij-1, :, :, 0))
                al_coh(ij, 1) = al_coh(ij, 1)/sum(phi_ld(ij-1, :, :, 1:ND))

            enddo

        endif

        ! calculate quantiles
        if(.not. analytical)call calculate_quantiles

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none

        integer :: ij, ixmax(JJ), iamax(JJ), idmax(JJ), ilmax(JJ), ages(JJ)

        ! check for the maximium grid points used
        call check_grid_X(ixmax); call check_grid_a(iamax); call check_grid_d(idmax); call check_grid_l(ilmax)

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)

        write(21, '(a,a,a)')' IJ      LINC      CONS      CONH    ASSETS     HOUSE    LIQASS  RSKSHARE  OWNSHARE',&
            '    CONS/O    CONS/R    CONH/O    CONH/R   RSHRE/O   RSHRE/R  ASSETS/O  ASSETS/R', &
            '   LIQASS/O  LIQASS/R     IXMAX     IAMAX     ILMAX     IDMAX'
        do ij = 1, JJ
            write(21,'(i3,18f10.3,4i10)')ages(ij), y_coh(ij), c_coh(ij, 2), ch_coh(ij, 2), a_coh(ij, 2), h_coh(ij), &
            al_coh(ij, 2), omega_coh(ij, 2), o_coh(ij), c_coh(ij, 1), c_coh(ij, 0), ch_coh(ij,1), ch_coh(ij,0), &
            omega_coh(ij, 1), omega_coh(ij, 0), a_coh(ij, 1), a_coh(ij, 0), al_coh(ij, 1), al_coh(ij,0), &
                    ixmax(ij), iamax(ij), ilmax(ij), idmax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! polt homeownership ratio
        call plot(dble(ages), o_coh(:), legend='Homeownership')
        call plot(dble(ages), omega_coh(:, 2), legend='Omega_p ')
        call execplot(xlabel='Age j', ylabel='Homeownership', ylim=(/0d0, 1d0/))

        ! plot consumption for homeowner
        call plot(dble(ages), c_coh(:, 2), legend='Consumption')
        call plot(dble(ages), ch_coh(:, 2), legend='ConsumptionH')
        call plot(dble(ages), a_coh(:, 2), legend='Assets')
        call plot(dble(ages), y_coh(:), legend='Labor Income)')
        call execplot(xlabel='Age j', ylabel='Consumption/assets')


        ! polt share of risky assets
        call plot(dble(ages), omega_coh(:, 0), legend='Omega_p - Renter')
        call plot(dble(ages), omega_coh(:, 1), legend='Omega_p - Owner')
        call execplot(xlabel='Age j', ylabel='Omega_p', ylim=(/0d0, 1d0/))

        ! plot consumption for homeowner
        call plot(dble(ages), c_coh(:, 1), legend='Consumption - Owner')
        call plot(dble(ages), ch_coh(:, 1), legend='Housing     - Owner')
        call plot(dble(ages), a_coh(:, 1), legend='Assets      - Owner')
        call plot(dble(ages), al_coh(:, 1)-xi*h_coh(ij), legend='debt')
        call execplot(xlabel='Age j', ylabel='Consumption/assets')

        ! polt consumption for renter
        call plot(dble(ages), c_coh(:, 0), legend='Consumption - Renter')
        call plot(dble(ages), ph*ch_coh(:, 0), legend='Housing    - Renter')
        call plot(dble(ages), a_coh(:, 0), legend='Assets      - Renter')
!        call plot(dble(ages), y_coh(:), legend='Labor Income (Mean)')
        call execplot(xlabel='Age j', ylabel='Consumption/Income')

        ! plot assets for homeowner
!        call plot(dble(ages), h_coh(:, 1), legend='Real Estate (Mean) - Owner')
!        call execplot(xlabel='Age j', ylabel='Assets', ylim=(/0d0, 50d0/))

        ! plot assets for renter
!        call plot(dble(ages), a_coh(:, 0), legend='Total Assets (Mean) - Renter')
!        call execplot(xlabel='Age j', ylabel='Assets', ylim=(/0d0, 50d0/))

    end subroutine

end program
