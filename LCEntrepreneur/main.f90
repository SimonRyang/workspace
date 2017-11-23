include "toolbox.f90"
include "globals.f90"

program main

    ! load modules
    use globals
    use omp_lib

    implicit none

    integer, parameter :: numthreads = 4

    ! set government variables
    mu     = 1d0
    lambda = 0d0
    phi    = 1d0
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

        integer :: ij, ix, ip
        real*8 :: ann_temp

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

        theta = 0d0

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

        ! annuity payments
        ann = 0d0
        ann_temp = 1d0

        do ij = JJ-1, JR, -1
          ann_temp = ann_temp/(1d0+r)*psi(ij) + 1d0
        enddo
        do ix = 0, NX
          ann(ix, JR:JJ) = (1d0+r)/psi(JR)*x(ix)/ann_temp
        enddo

        ! old-age transfers
        pen = 0d0
        do ip = 0, NP
          pen(ip, JR:JJ) = p(ip)*kappa*w*eff(JR-1)
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

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none

        integer :: ij, iq, ia, ik, ix, ip, iw, ie, iq_p, ip_p

        ! solve household problem recursively

        do ij = JJ, 1, -1

            if(ij == JJ)then

              omega_k_t(:, :, :, :, :, :, :, JJ) = 0d0
              omega_x_t(:, :, :, :, :, :, :, JJ) = 0d0

              do iq_p = 0, NQ
                  S(0, iq_p, :, :, :, :, :, JJ) = mu_b*max(Q(iq_p), 1d-16)**egam/egam
              enddo

              do ip = 0, NP
                do ix = 0, NX
                  do ik = 0, NK
                    do ia = 0, NA

                      ! with bequest motive we assume future worker
                      call solve_consumption(0, ia, ik, ix, ip, 1, 1, JJ)

                      Q_plus(ia, ik, ix, ip, :, :, JJ) = Q_plus_t(0, ia, ik, ix, ip, 1, 1, JJ)
                      a_plus(ia, ik, ix, ip, :, :, JJ) = a_plus_t(0, ia, ik, ix, ip, 1, 1, JJ)
                      x_plus(ia, ik, ix, ip, :, :, JJ) = x_plus_t(0, ia, ik, ix, ip, 1, 1, JJ)
                      p_plus(ia, ik, ix, ip, :, :, JJ) = p_plus_t(0, ia, ik, ix, ip, 1, 1, JJ)
                      k_plus(ia, ik, ix, ip, :, :, JJ) = k_plus_t(0, ia, ik, ix, ip, 1, 1, JJ)
                      c(ia, ik, ix, ip, :, :, JJ) = c_t(0, ia, ik, ix, ip, 1, 1, JJ)
                      l(ia, ik, ix, ip, :, :, JJ) = l_t(0, ia, ik, ix, ip, 1, 1, JJ)
                      V(ia, ik, ix, ip, :, :, JJ) = V_t(0, ia, ik, ix, ip, 1, 1, JJ)

                    enddo
                  enddo
               enddo
             enddo

             else

               ! get optimal share of wealth invested into capital

                  !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads) shared(ij)
                   do ie = 1, NE
                     do iw = 1, NW
                       do ip_p = 0, NP
                         do ix = 0, NX
                           do ik = 0, NK
                             do iq_p = 0, NQ

                               if(ij >= JR) then

                                 ! next period retiree
                                 call solve_retiree(iq_p, ik, ix, ip_p, iw, ie, ij)

                               else

                                 ! next period worker
                                 call solve_worker(iq_p, ik, ix, ip_p, iw, ie, ij)

                                 ! next period entrepreneur
                                 call solve_entrepreneur(iq_p, ik, ix, ip_p, iw, ie, ij)

                               endif

                            enddo
                         enddo
                       enddo
                     enddo
                   enddo
               enddo
               !$omp end parallel do

               !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads) shared(ij)
               ! solve the consumption savings problem
               do ie = 1, NE
                 do iw = 1, NW
                   do ip = 0, NP
                     do ix = 0, NX
                       do ik = 0, NK
                         do ia = 0, NA

                           ! next period worker
                           call solve_consumption(0, ia, ik, ix, ip, iw, ie, ij)

                           ! next period entrpreneur
                           if(ij < JR-1) call solve_consumption(1, ia, ik, ix, ip, iw, ie, ij)

                           ! decision on whether to be homeowner or renter next period
                            if(ij < JR-1 .and. V_t(1, ia, ik, ix, ip, iw, ie, ij) > V_t(0, ia, ik, ix, ip, iw, ie, ij)) then
                                  Q_plus(ia, ik, ix, ip, iw, ie, ij) = Q_plus_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  a_plus(ia, ik, ix, ip, iw, ie, ij) = a_plus_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  x_plus(ia, ik, ix, ip, iw, ie, ij) = x_plus_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  p_plus(ia, ik, ix, ip, iw, ie, ij) = p_plus_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  k_plus(ia, ik, ix, ip, iw, ie, ij) = k_plus_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  c(ia, ik, ix, ip, iw, ie, ij) = c_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  l(ia, ik, ix, ip, iw, ie, ij) = l_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  V(ia, ik, ix, ip, iw, ie, ij) = V_t(1, ia, ik, ix, ip, iw, ie, ij)
                             else
                               Q_plus(ia, ik, ix, ip, iw, ie, ij) = Q_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               a_plus(ia, ik, ix, ip, iw, ie, ij) = a_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               x_plus(ia, ik, ix, ip, iw, ie, ij) = x_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               p_plus(ia, ik, ix, ip, iw, ie, ij) = p_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               k_plus(ia, ik, ix, ip, iw, ie, ij) = k_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               c(ia, ik, ix, ip, iw, ie, ij) = c_t(0, ia, ik, ix, ip, iw, ie, ij)
                               l(ia, ik, ix, ip, iw, ie, ij) = l_t(0, ia, ik, ix, ip, iw, ie, ij)
                               V(ia, ik, ix, ip, iw, ie, ij) = V_t(0, ia, ik, ix, ip, iw, ie, ij)
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

      integer :: ia, ik, ix, ip, iw, ie, iw_p, ie_p

      !$omp parallel do collapse(3) schedule(dynamic,1) private(ie_p, iw_p) num_threads(numthreads) shared(ij)
      do ie = 1, NE
        do iw = 1, NW
          do ip = 0, NP
            do ix = 0, NX
              do ik = 0, NK
                do ia = 0, NA

                  EV(ia, ik, ix, ip, iw, ie, ij) = 0d0
                  do ie_p = 1, NE
                    do iw_p = 1, NW
                      EV(ia, ik, ix, ip, iw, ie, ij) = EV(ia, ik, ix, ip, iw, ie, ij) &
                        + pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*V(ia, ik, ix, ip, iw, ie, ij)
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

        integer :: ia, ik, ix, ip, iw, ie, ij, iw_p, ie_p
        integer :: iql, iqr, ial, iar, ikl, ikr, ixl, ixr, ipl, ipr
        real*8 :: varphi_q, varphi_a, varphi_k, varphi_x, varphi_p

        m(:, :, :, :, :, :, :) = 0d0
        m_Q(:, :, :, :, :, :, :) = 0d0

        do iw = 1, NW
          do ie = 1, NE
              m(0, 0, 0, 0, iw, ie, 1) = dist_eta(iw)*dist_theta(ie)
          enddo
        enddo

        do ij = 2, JJ

          do ie = 1, NE
            do iw = 1, NW
              do ip = 0, NP
                do ix = 0, NX
                  do ik = 0, NK
                    do ia = 0, NA

                      ! skip if there is no household
                      if (m(ia, ik, ix, ip, iw, ie, ij-1) <= 0d0) cycle

                      ! derive interpolation weights
                      call linint_Grow(Q_plus(ia, ik, ix, ip, iw, ie, ij-1), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                      call linint_Grow(a_plus(ia, ik, ix, ip, iw, ie, ij-1), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                      call linint_Grow(k_plus(ia, ik, ix, ip, iw, ie, ij-1), k_l, k_u, k_grow, NK-1, ikl, ikr, varphi_k)
                      call linint_Grow(x_plus(ia, ik, ix, ip, iw, ie, ij-1), x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)
                      call linint_Equi(p_plus(ia, ik, ix, ip, iw, ie, ij-1), p_l, p_u, NP, ipl, ipr, varphi_p)

                      ! restrict values to grid just in case
                      iql = min(iql, NA)
                      iqr = min(iqr, NA)
                      varphi_q = max(min(varphi_q, 1d0),0d0)

                      ! restrict values to grid just in case
                      ial = min(ial, NA)
                      iar = min(iar, NA)
                      varphi_a = max(min(varphi_a, 1d0),0d0)

                      ! restrict values to grid just in case
                      if (k_plus(ia, ik, ix, ip, iw, ie, ij-1) >= k_min) then
                        ikl = min(ikl+1, NK)
                        ikr = min(ikr+1, NK)
                        varphi_k = max(min(varphi_k, 1d0), 0d0)
                      else
                        ikl = 0; ikr = 0; varphi_k = 1d0
                      endif

                      ! restrict values to grid just in case
                      ixl = min(ixl, NX)
                      ixr = min(ixr, NX)
                      varphi_x = max(min(varphi_x, 1d0),0d0)

                      ! restrict values to grid just in case
                      ipl = min(ipl, NP)
                      ipr = min(ipr, NP)
                      varphi_p = max(min(varphi_p, 1d0), 0d0)

                      do iw_p = 1, NW
                        do ie_p = 1, NE

                            m(ial, ikl, ixl, ipl, iw_p, ie_p, ij) = m(ial, ikl, ixl, ipl, iw_p, ie_p, ij) + &
                                  varphi_a*varphi_k*varphi_x*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(ial, ikl, ixl, ipr, iw_p, ie_p, ij) = m(ial, ikl, ixl, ipr, iw_p, ie_p, ij) + &
                                  varphi_a*varphi_k*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(ial, ikl, ixr, ipl, iw_p, ie_p, ij) = m(ial, ikl, ixr, ipl, iw_p, ie_p, ij) + &
                                  varphi_a*varphi_k*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(ial, ikl, ixr, ipr, iw_p, ie_p, ij) = m(ial, ikl, ixr, ipr, iw_p, ie_p, ij) + &
                                  varphi_a*varphi_k*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(ial, ikr, ixl, ipl, iw_p, ie_p, ij) = m(ial, ikr, ixl, ipl, iw_p, ie_p, ij) + &
                                  varphi_a*(1d0-varphi_k)*varphi_x*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(ial, ikr, ixl, ipr, iw_p, ie_p, ij) = m(ial, ikr, ixl, ipr, iw_p, ie_p, ij) + &
                                  varphi_a*(1d0-varphi_k)*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(ial, ikr, ixr, ipl, iw_p, ie_p, ij) = m(ial, ikr, ixr, ipl, iw_p, ie_p, ij) + &
                                  varphi_a*(1d0-varphi_k)*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(ial, ikr, ixr, ipr, iw_p, ie_p, ij) = m(ial, ikr, ixr, ipr, iw_p, ie_p, ij) + &
                                  varphi_a*(1d0-varphi_k)*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(iar, ikl, ixl, ipl, iw_p, ie_p, ij) = m(iar, ikl, ixl, ipl, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*varphi_k*varphi_x*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(iar, ikl, ixl, ipr, iw_p, ie_p, ij) = m(iar, ikl, ixl, ipr, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*varphi_k*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(iar, ikl, ixr, ipl, iw_p, ie_p, ij) = m(iar, ikl, ixr, ipl, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*varphi_k*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(iar, ikl, ixr, ipr, iw_p, ie_p, ij) = m(iar, ikl, ixr, ipr, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*varphi_k*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(iar, ikr, ixl, ipl, iw_p, ie_p, ij) = m(iar, ikr, ixl, ipl, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*(1d0-varphi_k)*varphi_x*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(iar, ikr, ixl, ipr, iw_p, ie_p, ij) = m(iar, ikr, ixl, ipr, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*(1d0-varphi_k)*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(iar, ikr, ixr, ipl, iw_p, ie_p, ij) = m(iar, ikr, ixr, ipl, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)
                            m(iar, ikr, ixr, ipr, iw_p, ie_p, ij) = m(iar, ikr, ixr, ipr, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)

                        enddo
                      enddo

                      m_Q(iql, ik, ix, ipl, iw, ie, ij-1) = m_Q(iql, ik, ix, ipl, iw, ie, ij-1) + &
                                  varphi_q*varphi_p*m(ia, ik, ix, ip, iw, ie, ij-1)
                      m_Q(iql, ik, ix, ipr, iw, ie, ij-1) = m_Q(iql, ik, ix, ipr, iw, ie, ij-1) + &
                                  varphi_q*varphi_p*m(ia, ik, ix, ip, iw, ie, ij-1)
                      m_Q(iqr, ik, ix, ipl, iw, ie, ij-1) = m_Q(iqr, ik, ix, ipl, iw, ie, ij-1) + &
                                  (1d0-varphi_q)*varphi_p*m(ia, ik, ix, ip, iw, ie, ij-1)
                      m_Q(iqr, ik, ix, ipr, iw, ie, ij-1) = m_Q(iqr, ik, ix, ipr, iw, ie, ij-1) + &
                                  (1d0-varphi_q)*(1d0-varphi_p)*m(ia, ik, ix, ip, iw, ie, ij-1)

                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo

          write(*,*)sum(m(:, :, :, :, :, :, ij-1))
          write(*,*)sum(m_Q(:, :, :, :, :, :, ij-1))

        enddo

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none

        integer :: ia, ik, ix, ip, iw, ie, ij

        ! calculate cohort averages
        c_coh = 0d0; y_coh = 0d0; l_coh = 0d0; o_coh = 0d0; a_coh = 0d0; x_coh = 0d0; k_coh = 0d0

        do ij = 1, JJ

          do ie = 1, NE
            do iw = 1, NW
              do ip = 0, NP
                do ix = 0, NX
                  do ik = 0, NK
                    do ia = 0, NA

                        ! skip if there is no household
                        if (m(ia, ik, ix, ip, iw, ie, ij) <= 0d0) cycle

                        if (k_plus(ia, ik, ix, ip, iw, ie, ij) > 0d0 .and. k_plus(ia, ik, ix, ip, iw, ie, ij) < k_min) write(*,*)k_plus(ia, ik, ix, ip, iw, ie, ij), a_plus(ia, ik, ix, ip, iw, ie, ij), Q_plus(ia, ik, ix, ip, iw, ie, ij)

                        if(ik == 0) then
                          c_coh(0, ij) = c_coh(0, ij) + c(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                          a_coh(0, ij) = a_coh(0, ij) + a(ia)*m(ia, ik, ix, ip, iw, ie, ij)
                          x_coh(0, ij) = x_coh(0, ij) + x(ix)*m(ia, ik, ix, ip, iw, ie, ij)
                          l_coh(0, ij) = l_coh(0, ij) + l(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                          y_coh(0, ij) = y_coh(0, ij) + w*eff(ij)*eta(iw)*l(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        else
                          c_coh(1, ij) = c_coh(1, ij) + c(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                          a_coh(1, ij) = a_coh(1, ij) + a(ia)*m(ia, ik, ix, ip, iw, ie, ij)
                          x_coh(1, ij) = x_coh(1, ij) + x(ix)*m(ia, ik, ix, ip, iw, ie, ij)
                          k_coh(ij) = k_coh(ij) + k(ik)*m(ia, ik, ix, ip, iw, ie, ij)
                          y_coh(1, ij) = y_coh(1, ij) + theta(ie)*(k(ik)**alpha*(eff(ij)*l(ia, ik, ix, ip, iw, ie, ij))**(1d0-alpha))**nu*m(ia, ik, ix, ip, iw, ie, ij)
                          l_coh(1, ij) = l_coh(1, ij) + l(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                          o_coh(ij) = o_coh(ij) + m(ia, ik, ix, ip, iw, ie, ij)
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo

            c_coh(0, ij) = c_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-16)
            c_coh(1, ij) = c_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-16)
            a_coh(0, ij) = a_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-16)
            a_coh(1, ij) = a_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-16)
            x_coh(0, ij) = x_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-16)
            x_coh(1, ij) = x_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-16)
            y_coh(0, ij) = y_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-16)
            y_coh(1, ij) = y_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-16)
            l_coh(0, ij) = l_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-16)
            l_coh(1, ij) = l_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-16)
            o_coh(ij) = o_coh(ij)/max(sum(m(:, :, :, :, :, :, ij)), 1d-16)
            k_coh(ij) = k_coh(ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-16)

        enddo

        call check_grid(iamax, ikmax, ixmax)
        write(*,*) maxval(iamax), maxval(ikmax), maxval(ixmax)

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
        call plot(dble(ages), c_coh(1, :), legend='Consumption  - Entrepreneur')
        call plot(dble(ages), a_coh(1, :), legend='Assets       - Entrepreneur')
        call plot(dble(ages), x_coh(1, :), legend='Annuities    - Entrepreneur')
        call plot(dble(ages), y_coh(1, :), legend='Income       - Entrepreneur')
        call plot(dble(ages), l_coh(1, :), legend='Labor        - Entrepreneur')
        call plot(dble(ages), k_coh(:),    legend='Investment   - Entrepreneur')
        call execplot(xlabel='Age j', ylabel='Consumption/Assets')

        ! polt consumption for renter
        call plot(dble(ages), c_coh(0, :), legend='Consumption  - Worker')
        call plot(dble(ages), a_coh(0, :), legend='Assets       - Worker')
        call plot(dble(ages), x_coh(0, :), legend='Annuities    - Worker')
        call plot(dble(ages), y_coh(0, :), legend='Labor Income - Worker')
        call plot(dble(ages), l_coh(0, :), legend='Labor        - Worker')
        call execplot(xlabel='Age j', ylabel='Consumption/Assets')

    end subroutine

end program
