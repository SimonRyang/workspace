!##############################################################################
! PROGRAM Portfolio Choice with housing decision
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################
!
include "toolbox.f90"
include "globals.f90"
include "entrepreneur_foc.f90"
include "entrepreneur_solve.f90"
include "worker_foc.f90"
include "worker_solve.f90"

program LifecycleEntrepreneur

    ! load modules
    use globals
    use entrepreneur_solve
    use entrepreneur_foc
    use worker_solve
    use worker_foc

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

        ! wage rate for effective labor and rental price
        w = 1d0

        ! set survival probabilities
        psi(:) = (/1.00000d0, 0.99923d0, 0.99914d0, 0.99914d0, 0.99912d0, &
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

        ! initialize age earnings process
        eff(1:JR-1) = (/1.0000d0, 1.0719d0, 1.1438d0, 1.2158d0, 1.2842d0, 1.3527d0, &
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
        pen(JR:JJ) = kappa*w*eff(JR-1)

        ! discretize eta shocks
        call discretize_AR(rho, 0d0, sigma_eta, eta, pi_eta, dist_eta)
        eta = exp(eta)

        ! discretize theta shocks
        call discretize_AR(rho, 0d0, sigma_theta, theta, pi_theta, dist_theta)
        theta = exp(theta)
        theta = 0d0

        ! initialize asset grid
        call grid_Cons_Grow(X, X_l, X_u, X_grow)

        ! initialize liquid asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! endogenous upper bound of housing grid
        call grid_Cons_Grow(k, k_l, k_u, k_grow)

        ! initialize value functions
        S = 1d-10**egam/egam; V = 1d-10**egam/egam

        ! initialize policy functions
        X_plus = 0d0; a_plus = 0d0; k_plus = 0; c = 0d0
        omega_k = 0d0

        ! initialize temporary policy and value functions
        a_plus_t = 0d0; k_plus_t = 0d0; c_t = 0d0; V_t = 0d0

        ! open files
        open(21, file='output.out')


    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none

        integer :: ij, ix, ix_p, ia, ik, iw, ie

        ! solve household problem recursively

        do ij = JJ, 1, -1

            if(ij == JJ)then

              omega_k(JJ, :, :, :, :) = 0d0

              do ia = 0, NA
                do ik = 0, NA

                if (mu_b == 0d0) then

                  S(JJ, :, :, :, :, :) = 0d0

                  X_plus(JJ, ia, ik, :, :) = 0d0
                  a_plus(JJ, ia, ik, :, :) = 0d0
                  k_plus(JJ, ia, ik, :, :) = 0d0
                  cons_com = max((1d0+r)*a(ia) + pen(JJ), 1d-10)
                  c(JJ, ia, ik, :, :) = cons_com
                  V(JJ, ia, ik, :, :) = cons_com**egam/egam

                else


                  do ix_p = 1, NX
                      S(JJ, ix_p, :, :, :, :) = mu_b*X(ix_p)**egam/egam
                  enddo
                  S(JJ, 0, :, :, :, :) = 1d-10**egam/egam


                  ! with bequest motive we assume future worker
                  call solve_consumption_w(JJ, ia, ik, 1, 1)

                  X_plus(JJ, ia, ik, :, :) = X_plus_t(JJ, ia, ik, 1, 1, 0)
                  a_plus(JJ, ia, ik, :, :) = a_plus_t(JJ, ia, ik, 1, 1, 0)
                  k_plus(JJ, ia, ik, :, :) = k_plus_t(JJ, ia, ik, 1, 1, 0)
                  c(JJ, ia, ik, :, :) = c_t(JJ, ia, ik, 1, 1, 0)
                  V(JJ, ia, ik, :, :) = V_t(JJ, ia, ik, 1, 1, 0)

                 endif

               enddo
             enddo

             else

               ! get optimal share of wealth invested into capital

                   do ix_p = 0, NX
                       do ik = 0, NK
                         do iw = 1, NW
                           do ie = 1, NE

                           ! next period entrepreneur
                           !call solve_entrepreneur(ij, ix_p, ik, iw, ie)

                           ! next period worker
                           call solve_worker(ij, ix_p, ik, iw, ie)
                       enddo
                     enddo
                   enddo
               enddo


               ! solve the consumption savings problem
               do ia = 0, NA
                   do ik = 0, NK
                       do iw = 1, NW
                         do ie = 1, NE

                           ! next period entrpreneur
                           !call solve_consumption_e(ij, ia, ik, iw, ie)

                           ! next period worker
                           call solve_consumption_w(ij, ia, ik, iw, ie)


                         enddo

                       enddo
                   enddo
               enddo

               ! decision whether to be owner or renter next period
               do ia = 0, NA
                   do ik = 0, NK
                       do iw = 1, NW
                         do ie = 1, NE

                           ! decision on whether to be homeowner or renter next period
                           !if( V_t(ij, ia, ik, iw, ie, 1) > V_t(ij, ia, ik, iw, ie, 0) ) then
                                 k_plus(ij, ia, ik, iw, ie) = k_plus_t(ij, ia, ik, iw, ie, 1)
                                 a_plus(ij, ia, ik, iw, ie) = a_plus_t(ij, ia, ik, iw, ie, 1)
                                 c(ij, ia, ik, iw, ie) = c_t(ij, ia, ik, iw, ie, 1)
                                 V(ij, ia, ik, iw, ie) = V_t(ij, ia, ik, iw, ie, 1)
                          !  else
                          !    k_plus(ij, ia, ik, iw, ie) = k_plus_t(ij, ia, ik, iw, ie, 0)
                          !    a_plus(ij, ia, ik, iw, ie) = a_plus_t(ij, ia, ik, iw, ie, 0)
                          !    c(ij, ia, ik, iw, ie) = c_t(ij, ia, ik, iw, ie, 0)
                          !    V(ij, ia, ik, iw, ie) = V_t(ij, ia, ik, iw, ie, 0)
                          !  endif

                          enddo

                       enddo
                   enddo
               enddo

             endif

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

            call plot(a, V(ij, :, 0, 1, 1))
            call plot(a, V(ij, :, 3, 1, 1))
            call plot(a, V(ij, :, 5, 1, 1))
            call execplot

            call plot(X(1:), S(ij, 1:, 0, 1, 1, 0))
            call plot(X(1:), S(ij, 1:, 3, 1, 1, 0))
            call plot(X(1:), S(ij, 1:, 5, 1, 1, 0))
            call execplot

            call plot(X(1:), S(ij, 1:, 0, 1, 1, 1))
            call plot(X(1:), S(ij, 1:, 3, 1, 1, 1))
            call plot(X(1:), S(ij, 1:, 5, 1, 1, 1))
            call execplot

            call plot(X, X_plus(ij, :, 0, 1, 1))
            call plot(X, X_plus(ij, :, 3, 1, 1))
            call plot(X, X_plus(ij, :, 5, 1, 1))
            call execplot

            call plot(a, a_plus(ij, :, 0, 1, 1))
            call plot(a, a_plus(ij, :, 3, 1, 1))
            call plot(a, a_plus(ij, :, 5, 1, 1))
            call execplot

            call plot(a, c(ij, :, 0, 1, 1))
            call plot(a, c(ij, :, 3, 1, 1))
            call plot(a, c(ij, :, 5, 1, 1))
            call execplot

        enddo

    end subroutine


    ! determines the invariant distribution over cash-on-hand space
    subroutine get_distribution()

        implicit none

        integer :: ij, ia, ik, iw, iw_p, ie, ie_p
        integer :: ial, iar, ikl, ikr
        real*8 :: varphi_a, varphi_k

        m(:, :, :, :, :) = 0d0

        do iw = 1, NW
          do ie = 1, NE
              m(1, 0, 0, iw, ie) = dist_eta(iw)*dist_theta(ie)
          enddo
        enddo

        do ij = 2, JJ

          do ia = 0, NA
            do ik = 0, NK
              do iw = 1, NW
                do ie = 1, NE

                  ! skip if there is no household
                  if (m(ij-1, ia, ik, iw, ie) <= 0d0) cycle

                  ! derive interpolation weights
                  call linint_Grow(a_plus(ij-1, ia, ik, iw, ie), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                  call linint_Grow(k_plus(ij-1, ia, ik, iw, ie), k_l, k_u, k_grow, NK, ikl, ikr, varphi_k)

                  ! restrict values to grid just in case
                  ial = min(ial, NA)
                  iar = min(iar, NA)
                  varphi_a = max(min(varphi_a, 1d0),0d0)

                  ! restrict values to grid just in case
                  ikl = min(ikl, NK)
                  ikr = min(ikr, NK)
                  varphi_k = max(min(varphi_k, 1d0), 0d0)

                  do iw_p = 1, NW
                    do ie_p = 1, NE

                      if(varphi_a <= varphi_k)then
                          m(ij, ial, ikl, iw_p, ie_p) = m(ij, ial, ikl, iw_p, ie_p) + &
                                varphi_a*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ik, iw, ie)
                          m(ij, iar, ikl, iw_p, ie_p) = m(ij, iar, ikl, iw_p, ie_p) + &
                                (varphi_k-varphi_a)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ik, iw, ie)
                          m(ij, iar, ikr, iw_p, ie_p) = m(ij, iar, ikr, iw_p, ie_p) + &
                                (1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ik, iw, ie)
                      else
                        m(ij, ial, ikl, iw_p, ie_p) = m(ij, ial, ikl, iw_p, ie_p) + &
                              varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ik, iw, ie)
                        m(ij, ial, ikr, iw_p, ie_p) = m(ij, ial, ikr, iw_p, ie_p) + &
                              (varphi_a-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ik, iw, ie)
                        m(ij, iar, ikr, iw_p, ie_p) = m(ij, iar, ikr, iw_p, ie_p) + &
                              (1d0-varphi_a)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ik, iw, ie)
                      endif
                    enddo
                  enddo

                enddo
              enddo
            enddo
          enddo

        enddo

        write(*,*)sum(m(1, :, 0, :, :)), sum(m(2, :, 0, :, :)), sum(m(10, :, 0, :, :)), sum(m(40, :, 0, :, :))
        write(*,*)sum(m(1, :, 1:, :, :)), sum(m(2, :, 1:, :, :)), sum(m(10, :, 1:, :, :)), sum(m(40, :, 1:, :, :))
        write(*,*)sum(m(1, 1:, :, :, :)), sum(m(2, 1:, :, :, :)), sum(m(10, 1:, :, :, :)), sum(m(40, 1:, :, :, :))
    end subroutine



    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none

        integer :: ij, ia, ik, iw, ie

        ! calculate cohort averages
        c_coh = 0d0; y_coh = 0d0; o_coh = 0d0; a_coh = 0d0; k_coh = 0d0

        do ij = 1, JJ

            do ia = 0, NA
              do ik = 0, NK
                do iw = 1, NW
                  do ie = 1, NE
                    if(ik == 0) then
                      c_coh(ij, 0) = c_coh(ij,0) + c(ij, ia, ik, iw, ie)*m(ij, ia, ik, iw, ie)
                      a_coh(ij+1, 0) = a_coh(ij+1,0) + a_plus(ij, ia, ik, iw, ie)*m(ij, ia, ik, iw, ie)
                      y_coh(ij, 0) = y_coh(ij, 0) + w*eff(ij)*eta(iw)*m(ij, ia, ik, iw, ie)
                    else
                      c_coh(ij, 1) = c_coh(ij,1) + c(ij, ia, ik, iw, ie)*m(ij, ia, ik, iw, ie)
                      a_coh(ij+1, 1) = a_coh(ij+1,1) + a_plus(ij, ia, ik, iw, ie)*m(ij, ia, ik, iw, ie)
                      k_coh(ij) = k_coh(ij) + k(ik)*m(ij, ia, ik, iw, ie)
                      y_coh(ij, 1) = y_coh(ij, 1) + theta(ie)*k(ik)**nu*m(ij, ia, ik, iw, ie)
                      o_coh(ij) = o_coh(ij) + m(ij, ia, ik, iw, ie)
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo


            c_coh(ij, 0) = c_coh(ij, 0)/sum(m(ij, :, 0, :, :))
            c_coh(ij, 1) = c_coh(ij, 1)/sum(m(ij, :, 1:NK, :, :))
            a_coh(ij, 0) = a_coh(ij, 0)/sum(m(ij, :, 0, :, :))
            a_coh(ij, 1) = a_coh(ij, 1)/sum(m(ij, :, 1:NK, :, :))
            y_coh(ij, 0) = y_coh(ij, 0)/sum(m(ij, :, 0, :, :))
            y_coh(ij, 1) = y_coh(ij, 1)/sum(m(ij, :, 1:NK, :, :))
            o_coh(ij) = o_coh(ij)/sum(m(ij, :, :, :, :))
            k_coh(ij) = k_coh(ij)/sum(m(ij, :, :, :, :))


    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none

        integer :: ij, ages(JJ)
        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)


        ! polt homeownership ratio
        call plot(dble(ages), o_coh(:), legend='Homeownership')
        call execplot(xlabel='Age j', ylabel='Homeownership', ylim=(/0d0, 1d0/))

        ! plot consumption for homeowner
        call plot(dble(ages), c_coh(:, 1), legend='Consumption  - Owner')
        call plot(dble(ages), a_coh(:, 1), legend='Assets       - Owner')
        call plot(dble(ages), y_coh(:, 1), legend='Labor Income - Owner)')
        call execplot(xlabel='Age j', ylabel='Consumption/assets')

        ! polt consumption for renter
        call plot(dble(ages), c_coh(:, 0), legend='Consumption  - Renter')
        call plot(dble(ages), a_coh(:, 0), legend='Assets       - Renter')
        call plot(dble(ages), y_coh(:, 0), legend='Labor Income - Renter)')
        call execplot(xlabel='Age j', ylabel='Consumption/Income')

    end subroutine

end program
