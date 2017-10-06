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

        integer :: ij

        ! wage rate for effective labor and rental price
        w = 1d0

        ! set survival probabilities
        open(301, file='sp.dat')
        do ij = 1, JJ+1
          read(301,'(f13.8)')psi(ij)
        enddo
        close(301)

        ! initialize age earnings process
        eff(1:JR-1) = (/1.4327164d0, 1.8210024d0, 1.9747812d0, 2.0647004d0, 2.1559744d0, &
                        2.2020510d0, 2.2484878d0, 2.2359332d0, 2.1737906d0/)/1.4327164d0

        ! earnings process is during retirement equal to zero
        eff(JR:JJ) = 0d0

        ! old-age transfers
        pen = 0d0
        pen(JR:JJ) = kappa*w*eff(JR-1)

        ! discretize eta shocks
        call discretize_AR(rho, 0d0, sigma_eta, eta, pi_eta, dist_eta)
        call discretize_AR(0.95687d0**5d0, 0.0d0, sigma5(0.95687d0, 0.02812d0), eta, pi_eta, dist_eta)
        eta = exp(eta)

        ! discretize theta shocks
        call discretize_AR(rho, 0d0, sigma_theta, theta, pi_theta, dist_theta)
        call discretize_AR(0.920d0**5d0, 0.0d0, sigma5(0.920d0, 0.0375d0), theta, pi_theta, dist_theta)

        theta = exp(theta)

        ! initialize asset grid
        call grid_Cons_Grow(X, X_l, X_u, X_grow)

        ! initialize liquid asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! endogenous upper bound of housing grid
        call grid_Cons_Grow(k, k_l, k_u, k_grow)

        ! initialize value functions
        V = 1d-10**egam/egam; EV = 1d-10**egam/egam; S = 1d-10**egam/egam

        ! initialize policy functions
        X_plus = 0d0; a_plus = 0d0; k_plus = 0d0; c = 0d0; l = 0d0
        omega_k = 0d0

        ! initialize temporary policy and value functions
        X_plus_t = 0d0; a_plus_t = 0d0; k_plus_t = 0d0; c_t = 0d0; V_t = 0d0

        ! open files
        open(21, file='output.out')


    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none

        integer :: ij, ix, ix_p, ia, ik, iw, ie
        real*8:: y_plot(0:NX)

        ! solve household problem recursively

        do ij = JJ, 1, -1

            if(ij == JJ)then

              omega_k(JJ, :, :, :, :) = 0d0

              do ix_p = 0, NX
                  S(JJ, ix_p, :, :, :, :) = mu_b*max(X(ix_p), 1d-10)**egam/egam
              enddo

              do ia = 0, NA
                do ik = 0, NA

                  ! with bequest motive we assume future worker
                  call solve_consumption_w(JJ, ia, ik, 1, 1)

                  X_plus(JJ, ia, ik, :, :) = X_plus_t(JJ, ia, ik, 1, 1, 0)
                  a_plus(JJ, ia, ik, :, :) = a_plus_t(JJ, ia, ik, 1, 1, 0)
                  k_plus(JJ, ia, ik, :, :) = k_plus_t(JJ, ia, ik, 1, 1, 0)
                  c(JJ, ia, ik, :, :) = c_t(JJ, ia, ik, 1, 1, 0)
                  l(JJ, ia, ik, :, :) = l_t(JJ, ia, ik, 1, 1, 0)
                  V(JJ, ia, ik, :, :) = V_t(JJ, ia, ik, 1, 1, 0)

               enddo
             enddo

             else

               ! get optimal share of wealth invested into capital

                   do ix_p = 0, NX
                       do ik = 0, NK
                         do iw = 1, NW
                           do ie = 1, NE

                           ! next period entrepreneur
                           call solve_entrepreneur(ij, ix_p, ik, iw, ie)

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
                           call solve_consumption_e(ij, ia, ik, iw, ie)

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
                            if( V_t(ij, ia, ik, iw, ie, 1) > V_t(ij, ia, ik, iw, ie, 0) ) then
                                  X_plus(ij, ia, ik, iw, ie) = X_plus_t(ij, ia, ik, iw, ie, 1)
                                  k_plus(ij, ia, ik, iw, ie) = k_plus_t(ij, ia, ik, iw, ie, 1)
                                  a_plus(ij, ia, ik, iw, ie) = a_plus_t(ij, ia, ik, iw, ie, 1)
                                  c(ij, ia, ik, iw, ie) = c_t(ij, ia, ik, iw, ie, 1)
                                  l(ij, ia, ik, iw, ie) = l_t(ij, ia, ik, iw, ie, 1)
                                  V(ij, ia, ik, iw, ie) = V_t(ij, ia, ik, iw, ie, 1)
                           else
                             X_plus(ij, ia, ik, iw, ie) = X_plus_t(ij, ia, ik, iw, ie, 0)
                             k_plus(ij, ia, ik, iw, ie) = k_plus_t(ij, ia, ik, iw, ie, 0)
                             a_plus(ij, ia, ik, iw, ie) = a_plus_t(ij, ia, ik, iw, ie, 0)
                             c(ij, ia, ik, iw, ie) = c_t(ij, ia, ik, iw, ie, 0)
                             l(ij, ia, ik, iw, ie) = l_t(ij, ia, ik, iw, ie, 0)
                             V(ij, ia, ik, iw, ie) = V_t(ij, ia, ik, iw, ie, 0)
                           endif

                          enddo

                       enddo
                   enddo
               enddo

             endif

            call interpolate(ij)

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

            ! if (ij>40) cycle
            !
            ! ia_com = 4
            ! ik_com = 0
            ! iw_com = 3
            ! ie_com = 3
            ! do ix = 0, NX
            !   y_plot(ix) = cons_w(X(ix))
            ! enddo
            ! call plot(X(2:4), y_plot(2:4))
            ! call execplot
            !
            ! call solve_consumption_w(ij,ia_com,ik_com,iw_com,ie_com)
            ! write(*,*)X_plus_t(ij, ia_com, ik_com, iw_com, ie_com, 0)
            ! write(*,*)a_plus_t(ij, ia_com, ik_com, iw_com, ie_com, 0)
            ! write(*,*)c_t(ij, ia_com, ik_com, iw_com, ie_com, 0)
            !
            ! call plot(k(:), S(ij, 5, :, 1, 5, 0))
            ! call plot(k(:), S(ij, 5, :, 1, 5, 1))
            ! call execplot
            !
            ! call plot(a(1:), egam*V_t(ij, 1:, 0, 2, 1, 0)**(1d0/egam))
            ! call plot(a(1:), egam*V_t(ij, 1:, 0, 2, 1, 1)**(1d0/egam))
            ! call execplot
            !
            ! call plot(k(1:), egam*V_t(ij, 0, 1:, 1, 5, 0)**(1d0/egam))
            ! call plot(k(1:), egam*V_t(ij, 0, 1:, 1, 5, 1)**(1d0/egam))
            ! call execplot
            !
            ! call plot(a, EV(ij, :, 0, 3, 1))
            ! call plot(a, EV(ij, :, 0, 1, 1))
            ! call plot(a, EV(ij, :, 3, 1, 1))
            ! call plot(a, EV(ij, :, 5, 1, 1))
            ! call execplot
            !
            ! call plot(X(1:), S(ij, 1:, 0, 1, 1, 0))
            ! call plot(X(1:), S(ij, 1:, 3, 1, 1, 0))
            ! call plot(X(1:), S(ij, 1:, 5, 1, 1, 0))
            ! call execplot
            !
            ! call plot(X(1:), S(ij, 1:, 0, 1, 1, 1))
            ! call plot(X(1:), S(ij, 1:, 3, 1, 1, 1))
            ! call plot(X(1:), S(ij, 1:, 5, 1, 1, 1))
            ! call execplot
            !
            ! call plot(X, X_plus(ij, :, 0, 1, 1))
            ! call plot(X, X_plus(ij, :, 3, 1, 1))
            ! call plot(X, X_plus(ij, :, 5, 1, 1))
            ! call execplot
            !
            ! call plot(a, a_plus(ij, :, 0, 1, 1))
            ! call plot(a, a_plus(ij, :, 3, 1, 1))
            ! call plot(a, a_plus(ij, :, 5, 1, 1))
            ! call execplot
            !
            ! call plot(X, omega_k(ij, :, 0, 1, 1))
            ! call plot(X, omega_k(ij, :, 3, 1, 1))
            ! call plot(X, omega_k(ij, :, 5, 1, 1))
            ! call execplot
            !
            ! call plot(a, k_plus(ij, :, 0, 1, 1))
            ! call plot(a, k_plus(ij, :, 3, 1, 1))
            ! call plot(a, k_plus(ij, :, 5, 1, 1))
            ! call execplot
            !
            ! call plot(a, c(ij, :, 0, 1, 1))
            ! call plot(a, c(ij, :, 3, 1, 1))
            ! call plot(a, c(ij, :, 5, 1, 1))
            ! call execplot

        enddo

    end subroutine

    ! calculates the expected valuefunction of cohort ij
    subroutine interpolate(ij)

      implicit none

      integer, intent(in) :: ij

      integer :: ia, ik, iw, iw_p, ie, ie_p

      do ia = 0, NA
        do ik = 0, NK
          do iw = 1, NW
            do ie = 1, NE

              EV(ij, ia, ik, iw, ie) = 0d0
              do ie_p = 1, NE
                do iw_p = 1, NW
                  EV(ij, ia, ik, iw, ie) = EV(ij, ia, ik, iw, ie) &
                    + pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*V(ij, ia, ik, iw, ie)
                enddo ! iw_p
              enddo ! ie_p

            enddo
          enddo
        enddo
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

                    ! skip if there is no household
                    if (m(ij, ia, ik, iw, ie) <= 0d0) cycle

                    if(ik == 0) then
                      c_coh(ij, 0) = c_coh(ij,0) + c(ij, ia, ik, iw, ie)*m(ij, ia, ik, iw, ie)
                      a_coh(ij, 0) = a_coh(ij,0) + a(ia)*m(ij, ia, ik, iw, ie)
                      y_coh(ij, 0) = y_coh(ij, 0) + w*eff(ij)*eta(iw)*m(ij, ia, ik, iw, ie)
                    else
                      c_coh(ij, 1) = c_coh(ij,1) + c(ij, ia, ik, iw, ie)*m(ij, ia, ik, iw, ie)
                      a_coh(ij, 1) = a_coh(ij,1) + a(ia)*m(ij, ia, ik, iw, ie)
                      k_coh(ij) = k_coh(ij) + k(ik)*m(ij, ia, ik, iw, ie)
                      y_coh(ij, 1) = y_coh(ij, 1) + theta(ie)*(k(ik)**alpha*eff(ij)**(1d0-alpha))**nu*m(ij, ia, ik, iw, ie)
                      o_coh(ij) = o_coh(ij) + m(ij, ia, ik, iw, ie)
                    endif
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
            k_coh(ij) = k_coh(ij)/sum(m(ij, :, 1:NK, :, :))

        enddo

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none

        integer :: ij, ages(JJ)
        ! set up age variable
        ages = 20 + 5*(/(ij-1, ij=1,JJ)/)


        ! polt homeownership ratio
        call plot(dble(ages), o_coh(:), legend='Entrepreneurship')
        call execplot(xlabel='Age j', ylabel='Entrepreneurship', ylim=(/0d0, 1d0/))

        ! plot consumption for homeowner
        call plot(dble(ages), c_coh(:, 1), legend='Consumption  - Entrepreneur')
        call plot(dble(ages), a_coh(:, 1), legend='Assets       - Entrepreneur')
        call plot(dble(ages), y_coh(:, 1), legend='Income       - Entrepreneur')
        call plot(dble(ages), k_coh(:),    legend='Investment   - Entrepreneur')
        call execplot(xlabel='Age j', ylabel='Consumption/Assets')

        ! polt consumption for renter
        call plot(dble(ages), c_coh(:, 0), legend='Consumption  - Worker')
        call plot(dble(ages), a_coh(:, 0), legend='Assets       - Worker')
        call plot(dble(ages), y_coh(:, 0), legend='Labor Income - Worker')
        call execplot(xlabel='Age j', ylabel='Consumption/Assets')

    end subroutine

end program
