include "toolbox.f90"
include "globals.f90"

program main

    ! load modules
    use globals
    use omp_lib

    implicit none

    integer, parameter :: numthreads = 8

    ! set government variables
    mu     = 1d0
    lambda = 0d0
    phi    = 0d0
    taup   = 0.184d0

    ! initialize remaining variables
    call initialize()

    ! start the clock
    call tick(time)

    ! solve the household problem
    call solve_household()

    ! calculate the distribution of households over state space
    call get_distribution()

    ! aggregate individual decisions
    call aggregation()

    ! stop the clock
    call tock(time)

    call output()

    close(21)


contains


    ! initializes all remaining variables
    subroutine initialize

        implicit none

        integer :: ij, ip

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
                        2.2020510d0, 2.2484878d0, 2.2359332d0, 2.1737906d0/) !/1.4327164d0

        ! earnings process is during retirement equal to zero
        eff(JR:JJ) = 0d0

        ! discretize eta shocks
        call discretize_AR(0.95687d0**5d0, 0.0d0, sigma5(0.95687d0, 0.02812d0), eta, pi_eta, dist_eta)
        eta = exp(eta)

        ! discretize theta shocks
        call discretize_AR(0.920d0**5d0, 0.0d0, sigma5(0.920d0, 0.0375d0), theta, pi_theta, dist_theta)
        theta = exp(theta)

        ! initialize asset grid
        call grid_Cons_Grow(Q, Q_l, Q_u, Q_grow)

        ! initialize liquid asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize liquid annuity grid
        call grid_Cons_Grow(x, x_l, x_u, x_grow)

        ! initialize pension claim grid
        call grid_Cons_Equi(p, p_l, p_u)

        ! endogenous upper bound of housing grid
        call grid_Cons_Grow(k(1:NK), k_l, k_u, k_grow)
        k(0) = 0d0

        ! old-age transfers
        pen = 0d0
        do ip = 0, NP
          pen(JR:JJ, ip) = p(ip)*kappa*w*eff(JR-1)
        enddo

        ! initialize value functions
        V = 1d-16**egam/egam; EV = 1d-16**egam/egam; S = 1d-16**egam/egam

        ! initialize policy functions
        Q_plus = 0d0; a_plus = 0d0; x_plus = 0d0; p_plus = 0d0; k_plus = 0d0; c = 0d0; l = 0d0

        ! initialize temporary policy and value functions
        Q_plus_t = 0d0; a_plus_t = 0d0; x_plus_t = 0d0; p_plus_t = 0d0; k_plus_t = 0d0; c_t = 0d0
        omega_x_t = 0d0; omega_k_t = 0d0
        V_t = 0d0

        ! open files
        open(21, file='output.out')

        p_hat = 1d0
        do ij = JJ-1, JR, -1
          p_hat(ij) = p_hat(ij+1)/(1d0+r)*psi(ij) + 1d0
        enddo
        p_hat = 1d0/p_hat
        p_hat(1:JR-1) = 0d0

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none

        integer :: ij, iq, iq_p, ia, ix, ip, ip_p, ik, iw, ie
        real*8:: y_plot(0:1000)

        ! solve household problem recursively

        do ij = JJ, 1, -1

            if(ij == JJ)then

              omega_k_t(JJ, :, :, :, :, :, :, :) = 0d0
              omega_x_t(JJ, :, :, :, :, :, :, :) = 0d0

              do iq_p = 0, NQ
                  S(JJ, iq_p, :, :, :, :, :, 0) = mu_b*max(Q(iq_p), 1d-16)**egam/egam
              enddo

              do ia = 0, NA
                do ix = 0, NX
                  do ip = 0, NP
                    do ik = 0, NA

                      ! with bequest motive we assume future worker
                      call solve_consumption(JJ, ia, ix, ip, ik, 1, 1, 0)

                      Q_plus(JJ, ia, ix, ip, ik, :, :) = Q_plus_t(JJ, ia, ix, ip, ik, 1, 1, 0)
                      a_plus(JJ, ia, ix, ip, ik, :, :) = a_plus_t(JJ, ia, ix, ip, ik, 1, 1, 0)
                      x_plus(JJ, ia, ix, ip, ik, :, :) = x_plus_t(JJ, ia, ix, ip, ik, 1, 1, 0)
                      p_plus(JJ, ia, ix, ip, ik, :, :) = p_plus_t(JJ, ia, ix, ip, ik, 1, 1, 0)
                      k_plus(JJ, ia, ix, ip, ik, :, :) = k_plus_t(JJ, ia, ix, ip, ik, 1, 1, 0)
                      c(JJ, ia, ix, ip, ik, :, :) = c_t(JJ, ia, ix, ip, ik, 1, 1, 0)
                      l(JJ, ia, ix, ip, ik, :, :) = l_t(JJ, ia, ix, ip, ik, 1, 1, 0)
                      V(JJ, ia, ix, ip, ik, :, :) = V_t(JJ, ia, ix, ip, ik, 1, 1, 0)

                    enddo
                  enddo
               enddo
             enddo

             else

               ! get optimal share of wealth invested into capital

                  !$omp parallel do schedule(dynamic) num_threads(numthreads) shared(ij) default(none)
                   do iq_p = 0, NQ
                     do ix = 0, NX
                       do ip_p = 0, NP
                         do ik = 0, NK
                           do iw = 1, NW
                             do ie = 1, NE

                               if(ij >= JR) then

                                 ! next period retiree
                                 call solve_retiree(ij, iq_p, ix, ip_p, ik, iw, ie)

                               else

                                 ! next period worker
                                 call solve_worker(ij, iq_p, ix, ip_p, ik, iw, ie)

                                 ! next period entrepreneur
                                 call solve_entrepreneur(ij, iq_p, ix, ip_p, ik, iw, ie)

                               endif

                            enddo
                         enddo
                       enddo
                     enddo
                   enddo
               enddo
               !$omp end parallel do

                !$omp parallel do schedule(dynamic) num_threads(numthreads) shared(ij)
               ! solve the consumption savings problem
               do ia = 0, NA
                 do ix = 0, NX
                   do ip = 0, NP
                     do ik = 0, NK
                      do iw = 1, NW
                         do ie = 1, NE

                           ! next period worker
                           call solve_consumption(ij, ia, ix, ip, ik, iw, ie, 0)

                           ! next period entrpreneur
                           if(ij < JR-1) call solve_consumption(ij, ia, ix, ip, ik, iw, ie, 1)

                           ! decision on whether to be homeowner or renter next period
                            if(ij < JR-1 .and. V_t(ij, ia, ix, ip, ik, iw, ie, 1) > V_t(ij, ia, ix, ip, ik, iw, ie, 0)) then
                                  Q_plus(ij, ia, ix, ip, ik, iw, ie) = Q_plus_t(ij, ia, ix, ip, ik, iw, ie, 1)
                                  a_plus(ij, ia, ix, ip, ik, iw, ie) = a_plus_t(ij, ia, ix, ip, ik, iw, ie, 1)
                                  x_plus(ij, ia, ix, ip, ik, iw, ie) = x_plus_t(ij, ia, ix, ip, ik, iw, ie, 1)
                                  p_plus(ij, ia, ix, ip, ik, iw, ie) = p_plus_t(ij, ia, ix, ip, ik, iw, ie, 1)
                                  k_plus(ij, ia, ix, ip, ik, iw, ie) = k_plus_t(ij, ia, ix, ip, ik, iw, ie, 1)
                                  c(ij, ia, ix, ip, ik, iw, ie) = c_t(ij, ia, ix, ip, ik, iw, ie, 1)
                                  l(ij, ia, ix, ip, ik, iw, ie) = l_t(ij, ia, ix, ip, ik, iw, ie, 1)
                                  V(ij, ia, ix, ip, ik, iw, ie) = V_t(ij, ia, ix, ip, ik, iw, ie, 1)
                             else
                               Q_plus(ij, ia, ix, ip, ik, iw, ie) = Q_plus_t(ij, ia, ix, ip, ik, iw, ie, 0)
                               a_plus(ij, ia, ix, ip, ik, iw, ie) = a_plus_t(ij, ia, ix, ip, ik, iw, ie, 0)
                               x_plus(ij, ia, ix, ip, ik, iw, ie) = x_plus_t(ij, ia, ix, ip, ik, iw, ie, 0)
                               p_plus(ij, ia, ix, ip, ik, iw, ie) = p_plus_t(ij, ia, ix, ip, ik, iw, ie, 0)
                               k_plus(ij, ia, ix, ip, ik, iw, ie) = k_plus_t(ij, ia, ix, ip, ik, iw, ie, 0)
                               c(ij, ia, ix, ip, ik, iw, ie) = c_t(ij, ia, ix, ip, ik, iw, ie, 0)
                               l(ij, ia, ix, ip, ik, iw, ie) = l_t(ij, ia, ix, ip, ik, iw, ie, 0)
                               V(ij, ia, ix, ip, ik, iw, ie) = V_t(ij, ia, ix, ip, ik, iw, ie, 0)
                             endif

                         enddo
                       enddo
                     enddo
                   enddo
                 enddo
               enddo
             !$omp end parallel do

             endif

            call interpolate(ij)

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

        enddo

    end subroutine

    ! calculates the expected valuefunction of cohort ij
    subroutine interpolate(ij)

      implicit none

      integer, intent(in) :: ij

      integer :: ia, ix, ip, ik, iw, iw_p, ie, ie_p

      !$omp parallel do schedule(dynamic,1) private(ie_p, iw_p) num_threads(numthreads)
      do ia = 0, NA
        do ix = 0, NX
          do ip = 0, NP
            do ik = 0, NK
              do iw = 1, NW
                do ie = 1, NE

                  EV(ij, ia, ix, ip, ik, iw, ie) = 0d0
                  do ie_p = 1, NE
                    do iw_p = 1, NW
                      EV(ij, ia, ix, ip, ik, iw, ie) = EV(ij, ia, ix, ip, ik, iw, ie) &
                        + pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*V(ij, ia, ix, ip, ik, iw, ie)
                    enddo ! iw_p
                  enddo ! ie_p

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end parallel do

    end subroutine


    ! determines the invariant distribution over cash-on-hand space
    subroutine get_distribution()

        implicit none

        integer :: ij, ia, ix, ip, ik, iw, iw_p, ie, ie_p
        integer :: ial, iar, ixl, ixr, ipl, ipr, ikl, ikr
        real*8 :: varphi_a, varphi_x, varphi_p, varphi_k

        m(:, :, :, :, :, :, :) = 0d0

        do iw = 1, NW
          do ie = 1, NE
              m(1, 0, 0, 0, 0, iw, ie) = dist_eta(iw)*dist_theta(ie)
          enddo
        enddo

        do ij = 2, JJ

          do ia = 0, NA
            do ix = 0, NX
              do ip = 0, NP
                do ik = 0, NK
                  do iw = 1, NW
                    do ie = 1, NE

                      ! skip if there is no household
                      if (m(ij-1, ia, ix, ip, ik, iw, ie) <= 0d0) cycle

                      ! derive interpolation weights
                      call linint_Grow(a_plus(ij-1, ia, ix, ip, ik, iw, ie), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                      call linint_Grow(x_plus(ij-1, ia, ix, ip, ik, iw, ie), x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)
                      call linint_Equi(p_plus(ij-1, ia, ix, ip, ik, iw, ie), p_l, p_u, NP, ipl, ipr, varphi_p)
                      call linint_Grow(k_plus(ij-1, ia, ix, ip, ik, iw, ie), k_l, k_u, k_grow, NK-1, ikl, ikr, varphi_k)

                      ! restrict values to grid just in case
                      ial = min(ial, NA)
                      iar = min(iar, NA)
                      varphi_a = max(min(varphi_a, 1d0),0d0)

                      ! restrict values to grid just in case
                      ixl = min(ixl, NX)
                      ixr = min(ixr, NX)
                      varphi_x = max(min(varphi_x, 1d0),0d0)

                      ! restrict values to grid just in case
                      ipl = min(ipl, NP)
                      ipr = min(ipr, NP)
                      varphi_p = max(min(varphi_p, 1d0), 0d0)

                      ! restrict values to grid just in case
                      if (k_plus(ij-1, ia, ix, ip, ik, iw, ie) > 0d0) then
                        ikl = min(ikl+1, NK)
                        ikr = min(ikr+1, NK)
                        varphi_k = max(min(varphi_k, 1d0), 0d0)
                      else
                        ikl = 0; ikr = 0; varphi_k = 1d0
                      endif

                      do iw_p = 1, NW
                        do ie_p = 1, NE

                            m(ij, ial, ixl, ipl, ikl, iw_p, ie_p) = m(ij, ial, ixl, ipl, ikl, iw_p, ie_p) + &
                                  varphi_a*varphi_x*varphi_p*varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, ial, ixl, ipl, ikr, iw_p, ie_p) = m(ij, ial, ixl, ipl, ikr, iw_p, ie_p) + &
                                  varphi_a*varphi_x*varphi_p*(1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, ial, ixl, ipr, ikl, iw_p, ie_p) = m(ij, ial, ixl, ipr, ikl, iw_p, ie_p) + &
                                  varphi_a*varphi_x*(1d0-varphi_p)*varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, ial, ixl, ipr, ikr, iw_p, ie_p) = m(ij, ial, ixl, ipr, ikr, iw_p, ie_p) + &
                                  varphi_a*varphi_x*(1d0-varphi_p)*(1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, ial, ixr, ipl, ikl, iw_p, ie_p) = m(ij, ial, ixr, ipl, ikl, iw_p, ie_p) + &
                                  varphi_a*(1d0-varphi_x)*varphi_p*varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, ial, ixr, ipl, ikr, iw_p, ie_p) = m(ij, ial, ixr, ipl, ikr, iw_p, ie_p) + &
                                  varphi_a*(1d0-varphi_x)*varphi_p*(1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, ial, ixr, ipr, ikl, iw_p, ie_p) = m(ij, ial, ixr, ipr, ikl, iw_p, ie_p) + &
                                  varphi_a*(1d0-varphi_x)*(1d0-varphi_p)*varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, ial, ixr, ipr, ikr, iw_p, ie_p) = m(ij, ial, ixr, ipr, ikr, iw_p, ie_p) + &
                                  varphi_a*(1d0-varphi_x)*(1d0-varphi_p)*(1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, iar, ixl, ipl, ikl, iw_p, ie_p) = m(ij, iar, ixl, ipl, ikl, iw_p, ie_p) + &
                                  (1d0-varphi_a)*varphi_x*varphi_p*varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, iar, ixl, ipl, ikr, iw_p, ie_p) = m(ij, iar, ixl, ipl, ikr, iw_p, ie_p) + &
                                  (1d0-varphi_a)*varphi_x*varphi_p*(1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, iar, ixl, ipr, ikl, iw_p, ie_p) = m(ij, iar, ixl, ipr, ikl, iw_p, ie_p) + &
                                  (1d0-varphi_a)*varphi_x*(1d0-varphi_p)*varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, iar, ixl, ipr, ikr, iw_p, ie_p) = m(ij, iar, ixl, ipr, ikr, iw_p, ie_p) + &
                                  (1d0-varphi_a)*varphi_x*(1d0-varphi_p)*(1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, iar, ixr, ipl, ikl, iw_p, ie_p) = m(ij, iar, ixr, ipl, ikl, iw_p, ie_p) + &
                                  (1d0-varphi_a)*(1d0-varphi_x)*varphi_p*varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, iar, ixr, ipl, ikr, iw_p, ie_p) = m(ij, iar, ixr, ipl, ikr, iw_p, ie_p) + &
                                  (1d0-varphi_a)*(1d0-varphi_x)*varphi_p*(1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, iar, ixr, ipr, ikl, iw_p, ie_p) = m(ij, iar, ixr, ipr, ikl, iw_p, ie_p) + &
                                  (1d0-varphi_a)*(1d0-varphi_x)*(1d0-varphi_p)*varphi_k*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)
                            m(ij, iar, ixr, ipr, ikr, iw_p, ie_p) = m(ij, iar, ixr, ipr, ikr, iw_p, ie_p) + &
                                  (1d0-varphi_a)*(1d0-varphi_x)*(1d0-varphi_p)*(1d0-varphi_k)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ij-1, ia, ix, ip, ik, iw, ie)

                        enddo
                      enddo

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

        integer :: ij, ia, ix, ip, ik, iw, ie

        ! calculate cohort averages
        c_coh = 0d0; y_coh = 0d0; l_coh = 0d0; o_coh = 0d0; a_coh = 0d0; x_coh = 0d0; k_coh = 0d0

        do ij = 1, JJ

            do ia = 0, NA
              do ix = 0, NX
                do ip = 0, NP
                  do ik = 0, NK
                    do iw = 1, NW
                      do ie = 1, NE

                        ! skip if there is no household
                        if (m(ij, ia, ix, ip, ik, iw, ie) <= 0d0) cycle

                        if (k_plus(ij, ia, ix, ip, ik, iw, ie) > 0d0 .and. k_plus(ij, ia, ix, ip, ik, iw, ie) < k_min) write(*,*)k_plus(ij, ia, ix, ip, ik, iw, ie), a_plus(ij, ia, ix, ip, ik, iw, ie), Q_plus(ij, ia, ix, ip, ik, iw, ie)

                        if(ik == 0) then
                          c_coh(ij, 0) = c_coh(ij, 0) + c(ij, ia, ix, ip, ik, iw, ie)*m(ij, ia, ix, ip, ik, iw, ie)
                          a_coh(ij, 0) = a_coh(ij, 0) + a(ia)*m(ij, ia, ix, ip, ik, iw, ie)
                          x_coh(ij, 0) = x_coh(ij, 0) + x(ix)*m(ij, ia, ix, ip, ik, iw, ie)
                          l_coh(ij, 0) = l_coh(ij, 0) + l(ij, ia, ix, ip, ik, iw, ie)*m(ij, ia, ix, ip, ik, iw, ie)
                          y_coh(ij, 0) = y_coh(ij, 0) + w*eff(ij)*eta(iw)*l(ij, ia, ix, ip, ik, iw, ie)*m(ij, ia, ix, ip, ik, iw, ie)
                        else
                          c_coh(ij, 1) = c_coh(ij, 1) + c(ij, ia, ix, ip, ik, iw, ie)*m(ij, ia, ix, ip, ik, iw, ie)
                          a_coh(ij, 1) = a_coh(ij, 1) + a(ia)*m(ij, ia, ix, ip, ik, iw, ie)
                          x_coh(ij, 1) = x_coh(ij, 1) + x(ix)*m(ij, ia, ix, ip, ik, iw, ie)
                          k_coh(ij) = k_coh(ij) + k(ik)*m(ij, ia, ix, ip, ik, iw, ie)
                          y_coh(ij, 1) = y_coh(ij, 1) + theta(ie)*(k(ik)**alpha*(eff(ij)*l(ij, ia, ix, ip, ik, iw, ie))**(1d0-alpha))**nu*m(ij, ia, ix, ip, ik, iw, ie)
                          l_coh(ij, 1) = l_coh(ij, 1) + l(ij, ia, ix, ip, ik, iw, ie)*m(ij, ia, ix, ip, ik, iw, ie)
                          o_coh(ij) = o_coh(ij) + m(ij, ia, ix, ip, ik, iw, ie)
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo

            c_coh(ij, 0) = c_coh(ij, 0)/max(sum(m(ij, :, :, :, 0, :, :)), 1d-16)
            c_coh(ij, 1) = c_coh(ij, 1)/max(sum(m(ij, :, :, :, 1:NK, :, :)), 1d-16)
            a_coh(ij, 0) = a_coh(ij, 0)/max(sum(m(ij, :, :, :, 0, :, :)), 1d-16)
            a_coh(ij, 1) = a_coh(ij, 1)/max(sum(m(ij, :, :, :, 1:NK, :, :)), 1d-16)
            y_coh(ij, 0) = y_coh(ij, 0)/max(sum(m(ij, :, :, :, 0, :, :)), 1d-16)
            y_coh(ij, 1) = y_coh(ij, 1)/max(sum(m(ij, :, :, :, 1:NK, :, :)), 1d-16)
            l_coh(ij, 0) = l_coh(ij, 0)/max(sum(m(ij, :, :, :, 0, :, :)), 1d-16)
            l_coh(ij, 1) = l_coh(ij, 1)/max(sum(m(ij, :, :, :, 1:NK, :, :)), 1d-16)
            o_coh(ij) = o_coh(ij)/max(sum(m(ij, :, :, :, :, :, :)), 1d-16)
            k_coh(ij) = k_coh(ij)/max(sum(m(ij, :, :, :, 1:NK, :, :)), 1d-16)

        enddo

        call check_grid(iamax, ixmax, ikmax)
        write(*,*) maxval(iamax), maxval(ixmax), maxval(ikmax)

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
        call plot(dble(ages), x_coh(:, 1), legend='Annuities    - Entrepreneur')
        call plot(dble(ages), y_coh(:, 1), legend='Income       - Entrepreneur')
        call plot(dble(ages), l_coh(:, 1), legend='Labor        - Entrepreneur')
        call plot(dble(ages), k_coh(:),    legend='Investment   - Entrepreneur')
        call execplot(xlabel='Age j', ylabel='Consumption/Assets')

        ! polt consumption for renter
        call plot(dble(ages), c_coh(:, 0), legend='Consumption  - Worker')
        call plot(dble(ages), a_coh(:, 0), legend='Assets       - Worker')
        call plot(dble(ages), x_coh(:, 0), legend='Annuities    - Worker')
        call plot(dble(ages), y_coh(:, 0), legend='Labor Income - Worker')
        call plot(dble(ages), l_coh(:, 0), legend='Labor        - Worker')
        call execplot(xlabel='Age j', ylabel='Consumption/Assets')

    end subroutine

end program
