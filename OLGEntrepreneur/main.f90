include "toolbox.f90"
include "globals.f90"

program main

    ! load modules
    use globals
    use omp_lib

    implicit none

    integer, parameter :: numthreads = 28

    ! set government variables
    mu     = 0d0
    lambda = 0d0
    phi    = 0d0

    ! initialize remaining variables
    call initialize()

    ! start the clock
    call tick(time)

    ! calculate initial equilibrium
    call get_SteadyState()

    ! stop the clock
    call tock(time)

    call output()

    close(21)


contains

    ! computes the initial steady state of the economy
    subroutine get_SteadyState()

        implicit none

        ! iterate until value function converges
        do iter = 1, itermax

            ! get new prices
            call get_prices()

            ! solve the household problem
            call solve_household()

            ! calculate the distribution of households over state space
            call get_distribution()

            ! aggregate individual decisions
            call aggregation()

            ! determine the government parameters
            call government()

            ! check maximum grid points used
            call check_grid(iqmax, iamax, ikmax, ixmax)

            write(*,'(i4,4i7,5f8.2,f16.8)')iter, maxval(iqmax), maxval(iamax), maxval(ikmax), maxval(ixmax),&
                                            (/5d0*KK, CC, II/)/YY*100d0, &
                                            ((1d0+r)**0.2d0-1d0)*100d0, w, DIFF/YY*100d0

            if(abs(DIFF/YY)*100d0 < sig) return

        enddo

        write(*,*)'No Convergence'

    end subroutine

    ! initializes all remaining variables
    subroutine initialize

        implicit none

        integer :: ij, ix, ip
        real*8 :: ann_temp

        ! set survival probabilities
        open(301, file='sp.dat')
        do ij = 1, JJ+1
          read(301,'(f13.8)')psi(ij)
        enddo
        close(301)

        ! set up population structure
        rpop(0) = 1d0+n_p
        do ij = 1, JJ+1
            rpop(ij) = rpop(ij-1)*psi(ij)/(1d0+n_p)
        enddo

        workpop = 0d0
        ! initialize workforce
        do ij = 1, JJ
           if(ij < JR)workpop = workpop + rpop(ij)
        enddo

        ! initialize age earnings process
        eff(1:JR-1) = (/1.4327164d0, 1.8210024d0, 1.9747812d0, 2.0647004d0, 2.1559744d0, &
                        2.2020510d0, 2.2484878d0, 2.2359332d0, 2.1737906d0/)/1.4327164d0

        ! earnings process is during retirement equal to zero
        eff(JR:JJ) = 0d0

        ! discretize eta shocks
        call discretize_AR(0.95687d0**5d0, 0.0d0, sigma5(0.95687d0, 0.02812d0), eta, pi_eta, dist_eta)
        eta = exp(eta) !/sum(dist_eta*exp(eta))

        ! discretize theta shocks
        call discretize_AR(0.920d0**5d0, 0.0d0, sigma5(0.920d0, 0.0375d0), theta, pi_theta, dist_theta)
        theta = exp(theta)

        theta = 0d0

        ! initialize asset grid
        call grid_Cons_Grow(Q, Q_l, Q_u, Q_grow)

        ! initialize liquid asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! endogenous upper bound of housing grid
        call grid_Cons_Grow(k(1:NK), k_l, k_u, k_grow)
        k(0) = 0d0

        ! initialize annuity grid
        call grid_Cons_Grow(x, x_l, x_u, x_grow)
        write(*,*)x

        ! initialize pension claim grid
        call grid_Cons_Equi(p, p_l, p_u)

        ! annuity payments
        ann = 0d0
        ann_temp = 1d0

        do ij = JJ-1, JR, -1
          ann_temp = ann_temp/(1d0+r)*psi(ij) + 1d0
        enddo
        do ix = 0, NX
          ann(ix, JR:JJ) = (1d0+r)/psi(JR)*x(ix)/ann_temp
        enddo

        write(*,*) ann(:, JR)

        ! initialize tax rates
        taup  = 0.164d0

        ! set starting values
        KC = 7.30d0
        LC = 4.70d0
        BQ = 0.70d0

        ! initial guess bequests
        b = 0d0
        do ij = 1, JR-1
           b(ij) = BQ/workpop
        enddo

        ! initial guess average income
        ybar = 0.4d0

        ! initialize value functions
        V = 1d-13**egam/egam; EV = 1d-13**egam/egam; S = 1d-13**egam/egam

        ! initialize policy functions
        Q_plus = 0d0; a_plus = 0d0; x_plus = 0d0; p_plus = 0d0; k_plus = 0d0; c = 0d0; l = 0d0

        ! initialize temporary policy and value functions
        Q_plus_t = 0d0; a_plus_t = 0d0; x_plus_t = 0d0; p_plus_t = 0d0; k_plus_t = 0d0; c_t = 0d0
        omega_x_t = 0d0; omega_k_t = 0d0
        V_t = 0d0

        ! open files
        open(21, file='output.out')

    end subroutine


    ! compute prices and distribution of bequests for next iteration step
    subroutine get_prices()

        implicit none

        integer :: ij, ip

        ! calculate new prices
        r = Omega*alpha*(KC/LC)**(alpha-1d0)-delta_k
        w = Omega*(1d0-alpha)*(KC/LC)**alpha

        ! compute bequest per capita within workforce for next iteration step
        b(:) = 0d0
        do ij = 1, JR-1
           b(ij) = BQ/workpop
        enddo

        ! old-age transfers
        pen = 0d0
        do ip = 0, NP
          pen(ip, JR:JJ) = p(ip)*kappa*ybar
        enddo

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none

        integer :: ij, iq, ia, ik, ix, ip, iw, ie, iq_p, ip_p

        ! solve household problem recursively

        do ij = JJ, 1, -1

            if(ij == JJ)then

              omega_x_t(:, :, :, :, :, :, :, JJ) = 0d0
              omega_k_t(:, :, :, :, :, :, :, JJ) = 0d0

              do iq_p = 0, NQ
                  S(:, iq_p, :, :, :, :, :, JJ) = mu_b*max(Q(iq_p), 1d-13)**egam/egam
              enddo

              !$omp parallel do collapse(2) schedule(dynamic) num_threads(numthreads) shared(ij)
              do ip = 0, NP
                do ix = 0, NX
                  do ia = 0, NA

                    ! with bequest motive we assume future worker
                    call solve_consumption(0, ia, 0, ix, ip, 1, 1, JJ)

                    Q_plus(ia, :, ix, ip, :, :, JJ) = Q_plus_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    a_plus(ia, :, ix, ip, :, :, JJ) = a_plus_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    x_plus(ia, :, ix, ip, :, :, JJ) = x_plus_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    p_plus(ia, :, ix, ip, :, :, JJ) = p_plus_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    k_plus(ia, :, ix, ip, :, :, JJ) = k_plus_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    c(ia, :, ix, ip, :, :, JJ) = c_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    l(ia, :, ix, ip, :, :, JJ) = l_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    penb(ia, :, ix, ip, :, :, JJ) = penb_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    penc(ia, :, ix, ip, :, :, JJ) = penc_t(0, ia, 0, ix, ip, 1, 1, JJ)
                    V(ia, :, ix, ip, :, :, JJ) = V_t(0, ia, 0, ix, ip, 1, 1, JJ)

                  enddo
               enddo
             enddo
             !$omp end parallel do

           else
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
                                  penb(ia, ik, ix, ip, iw, ie, ij) = penb_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  penc(ia, ik, ix, ip, iw, ie, ij) = penc_t(1, ia, ik, ix, ip, iw, ie, ij)
                                  V(ia, ik, ix, ip, iw, ie, ij) = V_t(1, ia, ik, ix, ip, iw, ie, ij)
                             else
                               Q_plus(ia, ik, ix, ip, iw, ie, ij) = Q_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               a_plus(ia, ik, ix, ip, iw, ie, ij) = a_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               x_plus(ia, ik, ix, ip, iw, ie, ij) = x_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               p_plus(ia, ik, ix, ip, iw, ie, ij) = p_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               k_plus(ia, ik, ix, ip, iw, ie, ij) = k_plus_t(0, ia, ik, ix, ip, iw, ie, ij)
                               c(ia, ik, ix, ip, iw, ie, ij) = c_t(0, ia, ik, ix, ip, iw, ie, ij)
                               l(ia, ik, ix, ip, iw, ie, ij) = l_t(0, ia, ik, ix, ip, iw, ie, ij)
                               penb(ia, ik, ix, ip, iw, ie, ij) = penb_t(0, ia, ik, ix, ip, iw, ie, ij)
                               penc(ia, ik, ix, ip, iw, ie, ij) = penc_t(0, ia, ik, ix, ip, iw, ie, ij)
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

            !write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

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
                      call linint_Grow(Q_plus(ia, ik, ix, ip, iw, ie, ij-1), Q_l, Q_u, Q_grow, NQ, iql, iqr, varphi_q)
                      call linint_Grow(a_plus(ia, ik, ix, ip, iw, ie, ij-1), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                      call linint_Grow(k_plus(ia, ik, ix, ip, iw, ie, ij-1), k_l, k_u, k_grow, NK-1, ikl, ikr, varphi_k)
                      call linint_Grow(x_plus(ia, ik, ix, ip, iw, ie, ij-1), x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)
                      call linint_Equi(p_plus(ia, ik, ix, ip, iw, ie, ij-1), p_l, p_u, NP, ipl, ipr, varphi_p)

                      ! restrict values to grid just in case
                      iql = min(iql, NQ)
                      iqr = min(iqr, NQ)
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
                                  varphi_a*varphi_k*varphi_x*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(ial, ikl, ixl, ipr, iw_p, ie_p, ij) = m(ial, ikl, ixl, ipr, iw_p, ie_p, ij) + &
                                  varphi_a*varphi_k*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(ial, ikl, ixr, ipl, iw_p, ie_p, ij) = m(ial, ikl, ixr, ipl, iw_p, ie_p, ij) + &
                                  varphi_a*varphi_k*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(ial, ikl, ixr, ipr, iw_p, ie_p, ij) = m(ial, ikl, ixr, ipr, iw_p, ie_p, ij) + &
                                  varphi_a*varphi_k*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(ial, ikr, ixl, ipl, iw_p, ie_p, ij) = m(ial, ikr, ixl, ipl, iw_p, ie_p, ij) + &
                                  varphi_a*(1d0-varphi_k)*varphi_x*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(ial, ikr, ixl, ipr, iw_p, ie_p, ij) = m(ial, ikr, ixl, ipr, iw_p, ie_p, ij) + &
                                  varphi_a*(1d0-varphi_k)*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(ial, ikr, ixr, ipl, iw_p, ie_p, ij) = m(ial, ikr, ixr, ipl, iw_p, ie_p, ij) + &
                                  varphi_a*(1d0-varphi_k)*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(ial, ikr, ixr, ipr, iw_p, ie_p, ij) = m(ial, ikr, ixr, ipr, iw_p, ie_p, ij) + &
                                  varphi_a*(1d0-varphi_k)*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(iar, ikl, ixl, ipl, iw_p, ie_p, ij) = m(iar, ikl, ixl, ipl, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*varphi_k*varphi_x*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(iar, ikl, ixl, ipr, iw_p, ie_p, ij) = m(iar, ikl, ixl, ipr, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*varphi_k*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(iar, ikl, ixr, ipl, iw_p, ie_p, ij) = m(iar, ikl, ixr, ipl, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*varphi_k*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(iar, ikl, ixr, ipr, iw_p, ie_p, ij) = m(iar, ikl, ixr, ipr, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*varphi_k*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(iar, ikr, ixl, ipl, iw_p, ie_p, ij) = m(iar, ikr, ixl, ipl, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*(1d0-varphi_k)*varphi_x*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(iar, ikr, ixl, ipr, iw_p, ie_p, ij) = m(iar, ikr, ixl, ipr, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*(1d0-varphi_k)*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(iar, ikr, ixr, ipl, iw_p, ie_p, ij) = m(iar, ikr, ixr, ipl, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)
                            m(iar, ikr, ixr, ipr, iw_p, ie_p, ij) = m(iar, ikr, ixr, ipr, iw_p, ie_p, ij) + &
                                  (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij)*m(ia, ik, ix, ip, iw, ie, ij-1)/(1d0+n_p)

                        enddo
                      enddo

                      m_Q(iql, ik, ix, ipl, iw, ie, ij-1) = m_Q(iql, ik, ix, ipl, iw, ie, ij-1) + &
                                  varphi_q*varphi_p*m(ia, ik, ix, ip, iw, ie, ij-1)
                      m_Q(iql, ik, ix, ipr, iw, ie, ij-1) = m_Q(iql, ik, ix, ipr, iw, ie, ij-1) + &
                                  varphi_q*(1d0-varphi_p)*m(ia, ik, ix, ip, iw, ie, ij-1)
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

        enddo

       !  do ip = 0, NP
       !    write(*,*)sum(m(:, :, :, ip, :, :, 1))/sum(m(:, :, :, :, :, :, 1))
       !    write(*,*)sum(m(:, :, :, ip, :, :, 2))/sum(m(:, :, :, :, :, :, 2))
       !    write(*,*)sum(m(:, :, :, ip, :, :, 3))/sum(m(:, :, :, :, :, :, 3))
       !    write(*,*)sum(m(:, :, :, ip, :, :, 4))/sum(m(:, :, :, :, :, :, 4))
       !    write(*,*)sum(m(:, :, :, ip, :, :, 5))/sum(m(:, :, :, :, :, :, 5))
       !    write(*,*)sum(m(:, :, :, ip, :, :, 6))/sum(m(:, :, :, :, :, :, 6))
       !    write(*,*)sum(m(:, :, :, ip, :, :, 7))/sum(m(:, :, :, :, :, :, 7))
       !    write(*,*)sum(m(:, :, :, ip, :, :, 8))/sum(m(:, :, :, :, :, :, 8))
       !    write(*,*)sum(m(:, :, :, ip, :, :, 9))/sum(m(:, :, :, :, :, :, 9))
       !    write(*,*)sum(m(:, :, :, ip, :, :, JR))/sum(m(:, :, :, :, :, :, JR))
       !    write(*,*)sum(m(:, :, :, ip, :, :, JJ))/sum(m(:, :, :, :, :, :, JJ))
       ! enddo

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none

        integer :: ia, ik, ix, ip, iw, ie, ij
        real*8 :: LC_old, BQ_old
        real*8 :: Q_tmp(JJ), KC_tmp(JJ), Y_tmp(JJ), KE_tmp(JJ), BQ_tmp(JJ), PEN_tmp(JJ), TAUP_tmp(JJ), C_tmp(JJ)

        ! copy labor supply
        LC_old = LC
        BQ_old = BQ

        ! calculate cohort averages
        c_coh = 0d0; y_coh = 0d0; l_coh = 0d0; o_coh = 0d0; a_coh = 0d0; x_coh = 0d0; k_coh = 0d0; penb_coh = 0d0; penc_coh = 0d0

        ! reset macroeconomic aggregates in each iteration step
        AA = 0d0; BQ = 0d0; CC = 0d0; LC = 0d0; YE = 0d0; KE = 0d0; TC = 0d0; PBEN = 0d0; PCON = 0d0

        Q_tmp = 0d0; KC_tmp = 0d0; Y_tmp = 0d0; KE_tmp = 0d0; BQ_tmp = 0d0; PEN_tmp = 0d0; C_tmp = 0d0

        do ij = 1, JJ

          do ie = 1, NE
            do iw = 1, NW
              do ip = 0, NP
                do ix = 0, NX
                  do ik = 0, NK
                    do ia = 0, NA

                        Q_tmp(ij) = Q_tmp(ij) + Q_plus(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        KC_tmp(ij) = KC_tmp(ij) + (a(ia)-xi*k(ik))*m(ia, ik, ix, ip, iw, ie, ij)
                        KE_tmp(ij) = KE_tmp(ij) + k(ik)*m(ia, ik, ix, ip, iw, ie, ij)
                        BQ_tmp(ij) = BQ_tmp(ij) + b(ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        PEN_tmp(ij) = PEN_tmp(ij) + pen(ip, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        if (ik == 0) Y_tmp(ij) = Y_tmp(ij) + w*eff(ij)*eta(iw)*l(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        if (ik == 0) TAUP_tmp(ij) = TAUP_tmp(ij) + taup*w*eff(ij)*eta(iw)*l(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        if (ik > 0) Y_tmp(ij) = Y_tmp(ij) + theta(ie)*(k(ik)**alpha*(eff(ij)*l(ia, ik, ix, ip, iw, ie, ij))**(1d0-alpha))**nu*m(ia, ik, ix, ip, iw, ie, ij)
                        C_tmp(ij) = C_tmp(ij) + c(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)

                        ! skip if there is no household
                        if (m(ia, ik, ix, ip, iw, ie, ij) <= 0d0) cycle

                        if (a_plus(ia, ik, ix, ip, iw, ie, ij) < 0d0) write(*,*)'a_plus:', k_plus(ia, ik, ix, ip, iw, ie, ij), a_plus(ia, ik, ix, ip, iw, ie, ij), Q_plus(ia, ik, ix, ip, iw, ie, ij)
                        if (a_plus(ia, ik, ix, ip, iw, ie, ij) < 0d0) write(*,*) ia, ik, ix, ip, iw, ie, ij, m(ia, ik, ix, ip, iw, ie, ij)
                        if (k_plus(ia, ik, ix, ip, iw, ie, ij) > 0d0 .and. k_plus(ia, ik, ix, ip, iw, ie, ij) < k_min) write(*,*)'k_plus:', k_plus(ia, ik, ix, ip, iw, ie, ij), a_plus(ia, ik, ix, ip, iw, ie, ij), Q_plus(ia, ik, ix, ip, iw, ie, ij)

                        AA = AA + (a_plus(ia, ik, ix, ip, iw, ie, ij)-xi*k_plus(ia, ik, ix, ip, iw, ie, ij))*psi(ij+1)*m(ia, ik, ix, ip, iw, ie, ij)/(1d0+n_p)
                        CC = CC + c(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        BQ = BQ + Q_plus(ia, ik, ix, ip, iw, ie, ij)*(1d0-psi(ij+1))*m(ia, ik, ix, ip, iw, ie, ij)
                        KE = KE + k(ik)*m(ia, ik, ix, ip, iw, ie, ij)
                        TC = TC + tr(k(ik), k_plus(ia, ik, ix, ip, iw, ie, ij))*psi(ij+1)*m(ia, ik,ix, ip, iw, ie, ij)
                        PBEN = PBEN + penb(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        PCON = PCON + penc(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)

                        penb_coh(ij) = penb_coh(ij) + penb(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        penc_coh(ij) = penc_coh(ij) + penc(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)

                        if(ik == 0) then
                          LC = LC + eff(ij)*eta(iw)*l(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                          c_coh(0, ij) = c_coh(0, ij) + c(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                          a_coh(0, ij) = a_coh(0, ij) + a(ia)*m(ia, ik, ix, ip, iw, ie, ij)
                          x_coh(0, ij) = x_coh(0, ij) + x(ix)*m(ia, ik, ix, ip, iw, ie, ij)
                          l_coh(0, ij) = l_coh(0, ij) + l(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                          y_coh(0, ij) = y_coh(0, ij) + w*eff(ij)*eta(iw)*l(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                        else
                          YE = YE + theta(ie)*(k(ik)**alpha*(eff(ij)*l(ia, ik, ix, ip, iw, ie, ij))**(1d0-alpha))**nu*m(ia, ik, ix, ip, iw, ie, ij)
                          c_coh(1, ij) = c_coh(1, ij) + c(ia, ik, ix, ip, iw, ie, ij)*m(ia, ik, ix, ip, iw, ie, ij)
                          a_coh(1, ij) = a_coh(1, ij) + (a(ia)-xi*k(ik))*m(ia, ik, ix, ip, iw, ie, ij)
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

            c_coh(0, ij) = c_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-13)
            c_coh(1, ij) = c_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-13)
            a_coh(0, ij) = a_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-13)
            a_coh(1, ij) = a_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-13)
            x_coh(0, ij) = x_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-13)
            x_coh(1, ij) = x_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-13)
            y_coh(0, ij) = y_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-13)
            y_coh(1, ij) = y_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-13)
            l_coh(0, ij) = l_coh(0, ij)/max(sum(m(:, 0, :, :, :, :, ij)), 1d-13)
            l_coh(1, ij) = l_coh(1, ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-13)
            penb_coh(ij) = penb_coh(ij)/max(sum(m(:, :, :, :, :, :, ij)), 1d-13)
            penc_coh(ij) = penc_coh(ij)/max(sum(m(:, :, :, :, :, :, ij)), 1d-13)
            o_coh(ij) = o_coh(ij)/max(sum(m(:, :, :, :, :, :, ij)), 1d-13)
            k_coh(ij) = k_coh(ij)/max(sum(m(:, 1:NK, :, :, :, :, ij)), 1d-13)

            !write(*,*)'tmp:', Q_tmp(ij), KC_tmp(ij), KE_tmp(ij), Y_tmp(ij), BQ_tmp(ij), PEN_tmp(ij), TAUP_tmp(ij), C_tmp(ij)
            !write(*,*)Q_tmp(ij) - (1d0+r)*KC_tmp(ij) - (1d0-delta_k)*KE_tmp(ij) - Y_tmp(ij) - BQ_tmp(ij) - PEN_tmp(ij) + TAUP_tmp(ij) + C_tmp(ij)

        enddo ! ij

        write(*,*)'sum:', sum(Q_tmp(:) - (1d0+r)*KC_tmp(:) - (1d0-delta_k)*KE_tmp(:) - Y_tmp(:) - BQ_tmp(:) - PEN_tmp(:) + TAUP_tmp(:) + C_tmp(:))

        ! get average income
        ybar = w*LC/workpop

        ! compute stock of capital
        KC = damp*AA +(1d0-damp)*KC
        KK =  KC + KE
        LC = damp*LC +(1d0-damp)*LC_old
        BQ = damp*BQ +(1d0-damp)*BQ_old
        II = (n_p+delta_k)*KK
        YC = Omega*KC**alpha*LC**(1d0-alpha)
        YY = YC + YE

        ! compute gap on goods market
        DIFF = YY-CC-II-TC

        write(*,*)'pen', PBEN, PCON, w*LC, ybar, PBEN/PCON, taup

        write(*,*) YC, YE, KC, KE, LC, BQ

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government()

        implicit none

        real*8 :: taup_old

        taup_old = taup

        ! get budget balancing pension contribution rate
        taup = PBEN/PCON

        !taup = damp*taup + (1d0-damp)*taup_old

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
