include "toolbox.f90"
include "globals.f90"

program main

  ! load modules
  use globals
  use omp_lib

  implicit none

  integer, parameter :: numthreads = 28

  ! set government variables
  mu     = 1d0
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

  ! write output
  call output()

  call get_SteadyState()

  close(21)


contains


  !#############################################################################
  ! FUNCTION get_SteadyState
  !
  ! computes the initial steady state of the economy
  !#############################################################################
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

      write(*,'(i4,4i5,5f8.2,f16.8)')iter, maxval(iqmax), maxval(iamax), maxval(ikmax), maxval(ixmax),&
                                      (/5d0*KK, CC, II/)/YY*100d0, &
                                      ((1d0+r)**0.2d0-1d0)*100d0, w, DIFF/YY*100d0

      if(abs(DIFF/YY)*100d0 < sig) return

    enddo

    write(*,*)'No Convergence'

  end subroutine


  !#############################################################################
  ! SUBROUTINE initialize
  !
  ! initializes all remaining variables
  !#############################################################################
  subroutine initialize()

    implicit none

    !##### OTHER VARIABLES #####################################################
    integer :: iq, ia, ik, ix, ip, iw, ie, is, ij

    ! set survival probabilities
    open(301, file='sp.dat')
    do ij = 1, JJ+1
      read(301,'(f13.8)')psi(2, ij)
    enddo
    close(301)

    ! compute survival probabilities for high/low skilled
    psi(:, 1) = psi(2, 1)
    psi(:, JJ+1) = 0d0
    do ij = 2, JJ
      psi(1, ij) = psi(2, ij) - exp(0.33d0*(dble(ij-1)-22d0))
      psi(3, ij) = psi(2, ij) + exp(0.33d0*(dble(ij-1)-22d0))
    enddo

    ! set up population structure
    rpop(:, 1) = dist_skill(:)
    do ij = 2, JJ
      rpop(:, ij) = psi(:, ij)*rpop(:, ij-1)/(1d0+n_p)
    enddo

    ! set distribution of bequests
    Gama(1:JR-1) = 1d0
    Gama(JR:JJ) = 0d0
    Gama = Gama/sum(Gama)

    ! initialize age earnings process
    eff(1, 1:JR-1) = (/1.2987778d0, 1.5794954d0, 1.6404434d0, 1.6908550d0, 1.7507724d0, &
                       1.7586790d0, 1.7611338d0, 1.8054554d0, 1.7423268d0/)
    eff(2, 1:JR-1) = (/1.4327164d0, 1.8210024d0, 1.9747812d0, 2.0647004d0, 2.1559744d0, &
                      2.2020510d0, 2.2484878d0, 2.2359332d0, 2.1737906d0/)
    eff(3, 1:JR-1) = (/1.3882564d0, 2.1841104d0, 2.9655702d0, 3.3290738d0, 3.4171474d0, &
                       3.4497238d0, 3.4046532d0, 3.3062074d0, 3.1235630d0/)

    ! earnings process during retirement is equal to zero
    eff(:, JR:JJ) = 0d0

    ! initialize productivity process
    call discretize_AR(0.95666d0**5d0, 0.0d0, sigma5(0.95666d0, 0.02321d0), eta(:, 1), pi_eta(:, :, 1), dist_eta(:, 1))
    eta(:, 1) = exp(eta(:, 1))/sum(dist_eta(:, 1)*exp(eta(:, 1)))

    call discretize_AR(0.95687d0**5d0, 0.0d0, sigma5(0.95687d0, 0.02812d0), eta(:, 2), pi_eta(:, :, 2), dist_eta(:, 2))
    eta(:, 2) = exp(eta(:, 2))/sum(dist_eta(:, 2)*exp(eta(:, 2)))

    call discretize_AR(0.95828d0**5d0, 0.0d0, sigma5(0.95828d0, 0.03538d0), eta(:, 3), pi_eta(:, :, 3), dist_eta(:, 3))
    eta(:, 3) = exp(eta(:, 3))/sum(dist_eta(:, 3)*exp(eta(:, 3)))

    ! initialize entrepreneurial ability process
    theta(:, 1) = (/0d0, 0.95d0/)
    theta(:, 2) = (/0d0, 0.95d0/)
    theta(:, 3) = (/0d0, 0.95d0/)
    dist_theta(:, 1) = (/0d0, 1d0/)
    dist_theta(:, 2) = (/0d0, 1d0/)
    dist_theta(:, 3) = (/0d0, 1d0/)
    pi_theta(1, 1, :) = 1d0
    pi_theta(1, 2, :) = 0d0
    pi_theta(2, 1, :) = 0.2d0
    pi_theta(2, 2, :) = 0.8d0

    ! initialize asset grid
    call grid_Cons_Grow(Q, Q_l, Q_u, Q_grow)

    ! initialize liquid asset grid
    call grid_Cons_Grow(a, a_l, a_u, a_grow)

    ! endogenous upper bound of housing grid
    if (NK > 0) call grid_Cons_Grow(k(1:NK), k_l, k_u, k_grow)
    k(0) = 0d0

    ! initialize annuity grid
    if (NX > 0) call grid_Cons_Grow(x, x_l, x_u, x_grow)
    x(0) = x_l

    ! initialize pension claim grid
    call grid_Cons_Equi(p, p_l, p_u)

    ! get initial guess for household decisions
    omega_x_t = 0.05d0
    omega_k_t = 0.05d0
    l_t = 0.33d0
    do ia = 0, NA
      Q_plus_t(:, ia, :, :, :, :, :, :, :) = a(ia)/2d0
    enddo ! ia

    ! initialize tax rates
    tauc = 0.190d0
    taup = 0.189d0

    ! initial guesses for macro variables
    KC = 3.400d0
    LC = 3.604d0
    bqs(:) = (/4.610d-2, 0.180d0, 0.106d0/)
    BB = 2.964d0
    ybar = 0.555d0

    ! initialize value functions
    V = 1d-13**egam/egam; EV = 1d-13**egam/egam; S = 1d-13**egam/egam

    ! initialize policy functions
    Q_plus = 0d0; a_plus = 0d0; x_plus = 0d0; p_plus = 0d0; k_plus = 0d0; c = 0d0; l = 0d0

    ! initialize temporary policy and value functions
    a_plus_t = 0d0; x_plus_t = 0d0; p_plus_t = 0d0; k_plus_t = 0d0; c_t = 0d0
    V_t = 0d0

    ! open files
    open(21, file='output.out')

  end subroutine


  !#############################################################################
  ! SUBROUTINE get_prices
  !
  ! computes prices and distribution of bequests for next iteration step
  !#############################################################################
  subroutine get_prices()

      implicit none

      real*8 :: ann_tmp(NS)
      integer :: ix, ip, is, ij

      ! calculate new prices
      r = (1d0-tauy)*(Omega*alpha*(KC/LC)**(alpha-1d0)-delta_k)
      w = Omega*(1d0-alpha)*(KC/LC)**alpha

      ! set prices in case of life-cycle model
      ! r = 0.393280506035032d0
      ! w = 0.877841532937879d0
      ! bqs = (/4.608543623547606d-2, 0.181029882698876d0, 0.106845332164835d0/)
      ! ybar = 0.555719715351030d0
      ! tauc = 0.128579256047982d0
      ! taup = 7.867802841513299d-2

      ! calculate gross price of consumption (inverse)
      pinv = 1d0/(1d0+tauc)

      ! calculate individual bequests
      beq(1, :) = Gama(:)*bqs(1)/rpop(1, :)
      beq(2, :) = Gama(:)*bqs(2)/rpop(2, :)
      beq(3, :) = Gama(:)*bqs(3)/rpop(3, :)

      ! determine the income tax system
      r1 = 0.286d0*ybar*2d0
      r2 = 0.456d0*ybar*2d0
      r3 = 1.786d0*ybar*2d0

      b1 = (t2-t1)/(r2-r1)
      b2 = (t3-t2)/(r3-r2)

      ! calculate annuity payments
      ann = 0d0
      ans = 0d0
      ann_tmp = 1d0

      do ij = JJ-1, JR, -1
        ann_tmp(:) = ann_tmp(:)/(1d0+r)*psi(:, ij+1) + 1d0
      enddo

      do is = 1, NS
        do ix = 0, NX
          ann(ix, is, JR:JJ) = (1d0+r)/psi(is, JR)*x(ix)/ann_tmp(is)
        enddo

        do ij = 1, JR
          ans(:, is, ij) = x(:)
        enddo
        do ij = JR+1, JJ
          ans(:, is, ij) = (1d0+r)/psi(is, ij-1)*ans(:, is, ij-1)-ann(:, is, ij-1)
        enddo
      enddo

      ! calculate old-age transfers
      pen = 0d0
      do ip = 0, NP
        pen(ip, JR:JJ) = p(ip)*kappa*ybar
      enddo

  end subroutine


  !#############################################################################
  ! FUNCTION solve_household
  !
  ! determines the solution to the household optimization problem
  !#############################################################################
  subroutine solve_household()

    implicit none

    !##### OTHER VARIABLES #####################################################
    integer :: ij, iq, ia, ik, ix, ip, iw, ie, is, iq_p, ip_p, io_p

    ! solve household problem recursively

    omega_x_t(:, :, :, :, :, :, :, :, JJ) = 0d0
    omega_k_t(:, :, :, :, :, :, :, :, JJ) = 0d0

    do iq_p = 0, NQ
        S(:, iq_p, :, :, :, :, :, :, JJ) = mu_b*max(Q(iq_p), 1d-13)**egam/egam
    enddo ! iq_p

    !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads)
    do is = 1, NS
      do ip = 0, NP
        do ix = 0, NX
          do ia = 0, NA

            call solve_consumption(0, ia, 0, ix, ip, 1, 1, is, JJ)

            ! copy decisions
            Q_plus(ia, :, ix, ip, :, :, is, JJ) = Q_plus_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            a_plus(ia, :, ix, ip, :, :, is, JJ) = a_plus_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            x_plus(ia, :, ix, ip, :, :, is, JJ) = x_plus_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            p_plus(ia, :, ix, ip, :, :, is, JJ) = p_plus_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            k_plus(ia, :, ix, ip, :, :, is, JJ) = k_plus_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            c(ia, :, ix, ip, :, :, is, JJ) = c_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            l(ia, :, ix, ip, :, :, is, JJ) = l_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            inctax(ia, :, ix, ip, :, :, is, JJ) = inctax_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            captax(ia, :, ix, ip, :, :, is, JJ) = captax_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            penben(ia, :, ix, ip, :, :, is, JJ) = penben_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            pencon(ia, :, ix, ip, :, :, is, JJ) = pencon_t(0, ia, 0, ix, ip, 1, 1, is, JJ)
            V(ia, :, ix, ip, :, :, is, JJ) = V_t(0, ia, 0, ix, ip, 1, 1, is, JJ)

          enddo ! ia
        enddo ! ix
      enddo ! ip
    enddo ! is
    !$omp end parallel do

    call interpolate(JJ)

    ! solve for retirement age
    do ij = JJ-1, JR, -1

      !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads) shared(ij)
      do is = 1, NS
        do ip_p = 0, NP
          do ix = 0, NX
            do iq_p = 0, NQ

              ! next period retiree
              call solve_retiree(iq_p, 0, ix, ip_p, 1, 1, is, ij)

              omega_x_t(:, iq_p, :, ix, ip_p, :, :, is, ij) = omega_x_t(0, iq_p, 0, ix, ip_p, 1, 1, is, ij)
              omega_k_t(:, iq_p, :, ix, ip_p, :, :, is, ij) = omega_k_t(0, iq_p, 0, ix, ip_p, 1, 1, is, ij)
              S(:, iq_p, :, ix, ip_p, :, :, is, ij) = S(0, iq_p, 0, ix, ip_p, 1, 1, is, ij)

            enddo ! iq_p
          enddo ! ix
        enddo ! ip_p
      enddo ! is
      !$omp end parallel do

      !$omp parallel do collapse(2) schedule(dynamic) num_threads(numthreads) shared(ij)
      do is = 1, NS
        do ip = 0, NP
          do ix = 0, NX
            do ia = 0, NA

              ! next period worker
              call solve_consumption(0, ia, 0, ix, ip, 1, 1, is, ij)

              ! copy decisions
              Q_plus(ia, :, ix, ip, :, :, is, ij) = Q_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              a_plus(ia, :, ix, ip, :, :, is, ij) = a_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              x_plus(ia, :, ix, ip, :, :, is, ij) = x_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              p_plus(ia, :, ix, ip, :, :, is, ij) = p_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              k_plus(ia, :, ix, ip, :, :, is, ij) = k_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              c(ia, :, ix, ip, :, :, is, ij) = c_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              l(ia, :, ix, ip, :, :, is, ij) = l_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              inctax(ia, :, ix, ip, :, :, is, ij) = inctax_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              captax(ia, :, ix, ip, :, :, is, ij) = captax_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              penben(ia, :, ix, ip, :, :, is, ij) = penben_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              pencon(ia, :, ix, ip, :, :, is, ij) = pencon_t(0, ia, 0, ix, ip, 1, 1, is, ij)
              V(ia, :, ix, ip, :, :, is, ij) = V_t(0, ia, 0, ix, ip, 1, 1, is, ij)

            enddo ! ia
          enddo ! ix
        enddo ! ip
      enddo ! is
      !$omp end parallel do

      call interpolate(ij)
      !write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

    enddo ! ij

    do ij = JR-1, 1, -1

      !$omp parallel do collapse(4) schedule(dynamic) num_threads(numthreads) shared(ij)
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip_p = 0, NP
              do ix = 0, NX
                do ik = 0, NK
                  do iq_p = 0, NQ

                    ! next period worker
                    call solve_worker(iq_p, ik, ix, ip_p, iw, ie, is, ij)

                    ! next period entrepreneur
                    call solve_entrepreneur(iq_p, ik, ix, ip_p, iw, ie, is, ij)

                  enddo ! iq_p
                enddo ! ik
              enddo ! ix
            enddo ! ip_p
          enddo ! iw
        enddo ! ie
      enddo ! is
      !$omp end parallel do

      !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads) shared(ij)
      ! solve the consumption savings problem
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ik = 0, NK
                  do ia = 0, NA

                    ! next period worker
                    call solve_consumption(0, ia, ik, ix, ip, iw, ie, is, ij)

                    ! next period entrpreneur
                    if(ij < JR-1) call solve_consumption(1, ia, ik, ix, ip, iw, ie, is, ij)

                    ! decision on whether to be homeowner or renter next period
                    io_p = 0
                    if(ij < JR-1 .and. V_t(1, ia, ik, ix, ip, iw, ie, is, ij) > V_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)) io_p = 1

                    ! copy decisions
                    Q_plus(ia, ik, ix, ip, iw, ie, is, ij) = Q_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    a_plus(ia, ik, ix, ip, iw, ie, is, ij) = a_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    x_plus(ia, ik, ix, ip, iw, ie, is, ij) = x_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    p_plus(ia, ik, ix, ip, iw, ie, is, ij) = p_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    k_plus(ia, ik, ix, ip, iw, ie, is, ij) = k_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    c(ia, ik, ix, ip, iw, ie, is, ij) = c_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    l(ia, ik, ix, ip, iw, ie, is, ij) = l_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    inctax(ia, ik, ix, ip, iw, ie, is, ij) = inctax_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    captax(ia, ik, ix, ip, iw, ie, is, ij) = captax_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    penben(ia, ik, ix, ip, iw, ie, is, ij) = penben_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    pencon(ia, ik, ix, ip, iw, ie, is, ij) = pencon_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)
                    V(ia, ik, ix, ip, iw, ie, is, ij) = V_t(io_p, ia, ik, ix, ip, iw, ie, is, ij)

                  enddo ! ia
                enddo ! ik
              enddo ! ix
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is
      !$omp end parallel do

      call interpolate(ij)
      !write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

    enddo ! ij

  end subroutine


  !#############################################################################
  ! SUBROUTINE interpolate
  !
  ! calculates the expected valuefunction of cohort ij
  !#############################################################################
  subroutine interpolate(ij)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: ij

    !##### OTHER VARIABLES #####################################################
    integer :: ia, ik, ix, ip, iw, ie, is, iw_p, ie_p

    !$omp parallel do collapse(4) schedule(dynamic,1) private(ie_p, iw_p) num_threads(numthreads) shared(ij)
    do is = 1, NS
      do ie = 1, NE
        do iw = 1, NW
          do ip = 0, NP
            do ix = 0, NX
              do ik = 0, NK
                do ia = 0, NA

                  EV(ia, ik, ix, ip, iw, ie, is, ij) = 0d0
                  do ie_p = 1, NE
                    do iw_p = 1, NW
                      EV(ia, ik, ix, ip, iw, ie, is, ij) = EV(ia, ik, ix, ip, iw, ie, is, ij) &
                        + pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*V(ia, ik, ix, ip, iw, ie, is, ij)
                    enddo ! iw_p
                  enddo ! ie_p

                enddo ! ia
              enddo ! ik
            enddo ! ix
          enddo ! ip
        enddo ! iw
      enddo ! iw
    enddo ! is
    !$omp end parallel do

  end subroutine


  !#############################################################################
  ! SUBROUTINE get_distribution
  !
  ! determines the invariant distribution
  !#############################################################################
  subroutine get_distribution()

    implicit none

    !##### OTHER VARIABLES #####################################################
    integer :: ia, ik, ix, ip, iw, ie, is, ij, iw_p, ie_p
    integer :: iql, iqr, ial, iar, ikl, ikr, ixl, ixr, ipl, ipr
    real*8 :: varphi_q, varphi_a, varphi_k, varphi_x, varphi_p

    m(:, :, :, :, :, :, :, :) = 0d0
    m_Q(:, :, :, :, :, :, :, :) = 0d0

    do is = 1, NS
      do iw = 1, NW
        do ie = 1, NE
            m(0, 0, 0, 0, iw, ie, is, 1) = dist_eta(iw, is)*dist_theta(ie, is)*dist_skill(is)
        enddo
      enddo
    enddo

    do ij = 2, JJ

      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ik = 0, NK
                  do ia = 0, NA

                    ! skip if there is no household
                    if (m(ia, ik, ix, ip, iw, ie, is, ij-1) <= 0d0) cycle

                    ! derive interpolation weights
                    call linint_Grow(Q_plus(ia, ik, ix, ip, iw, ie, is, ij-1), Q_l, Q_u, Q_grow, NQ, iql, iqr, varphi_q)
                    call linint_Grow(a_plus(ia, ik, ix, ip, iw, ie, is, ij-1), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                    if (NK > 0) then
                      call linint_Grow(k_plus(ia, ik, ix, ip, iw, ie, is, ij-1), k_l, k_u, k_grow, NK-1, ikl, ikr, varphi_k)
                    else
                      ikl = 0; ikr = 0; varphi_k = 1d0
                    endif
                    if (NX > 0) then
                      call linint_Grow(x_plus(ia, ik, ix, ip, iw, ie, is, ij-1), x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)
                    else
                      ixl = 0; ixr = 0; varphi_x = 1d0
                    endif
                    call linint_Equi(p_plus(ia, ik, ix, ip, iw, ie, is, ij-1), p_l, p_u, NP, ipl, ipr, varphi_p)

                    ! restrict values to grid just in case
                    iql = min(iql, NQ)
                    iqr = min(iqr, NQ)
                    varphi_q = max(min(varphi_q, 1d0),0d0)

                    ! restrict values to grid just in case
                    ial = min(ial, NA)
                    iar = min(iar, NA)
                    varphi_a = max(min(varphi_a, 1d0),0d0)

                    ! restrict values to grid just in case
                    if (k_plus(ia, ik, ix, ip, iw, ie, is, ij-1) >= k_min) then
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

                        m(ial, ikl, ixl, ipl, iw_p, ie_p, is, ij) = m(ial, ikl, ixl, ipl, iw_p, ie_p, is, ij) + &
                              varphi_a*varphi_k*varphi_x*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(ial, ikl, ixl, ipr, iw_p, ie_p, is, ij) = m(ial, ikl, ixl, ipr, iw_p, ie_p, is, ij) + &
                              varphi_a*varphi_k*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(ial, ikl, ixr, ipl, iw_p, ie_p, is, ij) = m(ial, ikl, ixr, ipl, iw_p, ie_p, is, ij) + &
                              varphi_a*varphi_k*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(ial, ikl, ixr, ipr, iw_p, ie_p, is, ij) = m(ial, ikl, ixr, ipr, iw_p, ie_p, is, ij) + &
                              varphi_a*varphi_k*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(ial, ikr, ixl, ipl, iw_p, ie_p, is, ij) = m(ial, ikr, ixl, ipl, iw_p, ie_p, is, ij) + &
                              varphi_a*(1d0-varphi_k)*varphi_x*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(ial, ikr, ixl, ipr, iw_p, ie_p, is, ij) = m(ial, ikr, ixl, ipr, iw_p, ie_p, is, ij) + &
                              varphi_a*(1d0-varphi_k)*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(ial, ikr, ixr, ipl, iw_p, ie_p, is, ij) = m(ial, ikr, ixr, ipl, iw_p, ie_p, is, ij) + &
                              varphi_a*(1d0-varphi_k)*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(ial, ikr, ixr, ipr, iw_p, ie_p, is, ij) = m(ial, ikr, ixr, ipr, iw_p, ie_p, is, ij) + &
                              varphi_a*(1d0-varphi_k)*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(iar, ikl, ixl, ipl, iw_p, ie_p, is, ij) = m(iar, ikl, ixl, ipl, iw_p, ie_p, is, ij) + &
                              (1d0-varphi_a)*varphi_k*varphi_x*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(iar, ikl, ixl, ipr, iw_p, ie_p, is, ij) = m(iar, ikl, ixl, ipr, iw_p, ie_p, is, ij) + &
                              (1d0-varphi_a)*varphi_k*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(iar, ikl, ixr, ipl, iw_p, ie_p, is, ij) = m(iar, ikl, ixr, ipl, iw_p, ie_p, is, ij) + &
                              (1d0-varphi_a)*varphi_k*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(iar, ikl, ixr, ipr, iw_p, ie_p, is, ij) = m(iar, ikl, ixr, ipr, iw_p, ie_p, is, ij) + &
                              (1d0-varphi_a)*varphi_k*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(iar, ikr, ixl, ipl, iw_p, ie_p, is, ij) = m(iar, ikr, ixl, ipl, iw_p, ie_p, is, ij) + &
                              (1d0-varphi_a)*(1d0-varphi_k)*varphi_x*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(iar, ikr, ixl, ipr, iw_p, ie_p, is, ij) = m(iar, ikr, ixl, ipr, iw_p, ie_p, is, ij) + &
                              (1d0-varphi_a)*(1d0-varphi_k)*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(iar, ikr, ixr, ipl, iw_p, ie_p, is, ij) = m(iar, ikr, ixr, ipl, iw_p, ie_p, is, ij) + &
                              (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)
                        m(iar, ikr, ixr, ipr, iw_p, ie_p, is, ij) = m(iar, ikr, ixr, ipr, iw_p, ie_p, is, ij) + &
                              (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1)/(1d0+n_p)

                      enddo
                    enddo

                    m_Q(iql, ik, ix, ipl, iw, ie, is, ij-1) = m_Q(iql, ik, ix, ipl, iw, ie, is, ij-1) + &
                                varphi_q*varphi_p*m(ia, ik, ix, ip, iw, ie, is, ij-1)
                    m_Q(iql, ik, ix, ipr, iw, ie, is, ij-1) = m_Q(iql, ik, ix, ipr, iw, ie, is, ij-1) + &
                                varphi_q*(1d0-varphi_p)*m(ia, ik, ix, ip, iw, ie, is, ij-1)
                    m_Q(iqr, ik, ix, ipl, iw, ie, is, ij-1) = m_Q(iqr, ik, ix, ipl, iw, ie, is, ij-1) + &
                                (1d0-varphi_q)*varphi_p*m(ia, ik, ix, ip, iw, ie, is, ij-1)
                    m_Q(iqr, ik, ix, ipr, iw, ie, is, ij-1) = m_Q(iqr, ik, ix, ipr, iw, ie, is, ij-1) + &
                                (1d0-varphi_q)*(1d0-varphi_p)*m(ia, ik, ix, ip, iw, ie, is, ij-1)

                  enddo ! ia
                enddo ! ik
              enddo ! ix
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is

    enddo ! ij

  end subroutine


  !#############################################################################
  ! SUBROUTINE aggregation
  !
  ! calculate aggregated quantities
  !#############################################################################
  subroutine aggregation()

    implicit none

    !##### OTHER VARIABLES #####################################################
    integer :: ia, ik, ix, ip, iw, ie, is, ij
    real*8 :: LC_old

    ! copy labor supply
    LC_old = LC

    ! reset macroeconomic aggregates in each iteration step
    AA = 0d0; AX = 0d0; BQ = 0d0; bqs(:) = 0d0; PBEN = 0d0; PCON = 0d0
    CC = 0d0; LC = 0d0; YE = 0d0; KE = 0d0; TC = 0d0
    TAc = 0d0; TAr = 0d0; TAw = 0d0; TAy = 0d0

    do ij = 1, JJ

      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ik = 0, NK
                  do ia = 0, NA

                    ! skip if there is no household
                    if (m(ia, ik, ix, ip, iw, ie, is, ij) <= 0d0) cycle

                    AA = AA + (a_plus(ia, ik, ix, ip, iw, ie, is, ij)-xi*k_plus(ia, ik, ix, ip, iw, ie, is, ij))*psi(is, ij+1)*m(ia, ik, ix, ip, iw, ie, is, ij)/(1d0+n_p)
                    AX = AX + ans(ix, is, ij)/psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij)
                    CC = CC + c(ia, ik, ix, ip, iw, ie, is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij)
                    bqs(is) = bqs(is) + (a_plus(ia, ik, ix, ip, iw, ie, is, ij)+(1d0-xi)*k_plus(ia, ik, ix, ip, iw, ie, is, ij))*(1d0-psi(is, ij+1))*m(ia, ik, ix, ip, iw, ie, is, ij)
                    KE = KE + k(ik)*m(ia, ik, ix, ip, iw, ie, is, ij)
                    TC = TC + tr(k(ik), k_plus(ia, ik, ix, ip, iw, ie, is, ij))*m(ia, ik, ix, ip, iw, ie, is, ij)
                    TAc = TAc + tauc*c(ia, ik, ix, ip, iw, ie, is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij)
                    TAw = TAw + inctax(ia, ik, ix, ip, iw, ie, is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij)
                    TAr = TAr + captax(ia, ik, ix, ip, iw, ie, is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij)
                    PBEN = PBEN + penben(ia, ik, ix, ip, iw, ie, is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij)
                    PCON = PCON + pencon(ia, ik, ix, ip, iw, ie, is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij)

                    if(ik == 0) then
                      LC = LC + eff(is, ij)*eta(iw, is)*l(ia, ik, ix, ip, iw, ie, is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij)
                    else
                       YE = YE + theta(ie, is)*(k(ik)**alpha*(eff(is, ij)*l(ia, ik, ix, ip, iw, ie, is, ij))**(1d0-alpha))**nu*m(ia, ik, ix, ip, iw, ie, is, ij)
                    endif

                  enddo ! ia
                enddo ! ik
              enddo ! ix
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is

    enddo ! ij

    ! get average income
    ybar = w*LC/sum(m(:, :, :, :, :, :, :, 1:JR-1))

    ! compute stock of capital
    KC = damp*(AA+AX-BB) +(1d0-damp)*KC
    KK = KC + KE

    ! update work supply
    LC = damp*LC +(1d0-damp)*LC_old

    ! compute total bequests
    BQ = sum(bqs)

    ! commpute investment
    II = (n_p+delta_k)*KK

    ! compute output
    YC = Omega*KC**alpha*LC**(1d0-alpha)
    YY = YC + YE

    ! compute corporate tax incom
    TAy = tauy*(YC-delta_k*KC-w*LC)

  end subroutine


  !#############################################################################
  ! SUBROUTINE government
  !
  ! calculates government parameters
  !#############################################################################
  subroutine government()

    implicit none

    !##### OTHER VARIABLES #####################################################
    real*8 :: expend

    ! computes government expenditures
    GG = gy*YY
    BB = by*YY
    expend = GG + (1d0+r)*BB - (1d0+n_p)*BB

    ! calculates consumption tax rate
    tauc = (expend-TAy-TAw-TAr)/CC

    ! get budget balancing pension contribution rate
    taup = PBEN/PCON

    ! compute gap on goods market
    DIFF = YY-CC-II-TC-GG

  end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none

        integer :: is, ij
        real*8 :: life_exp(NS), punb(NS, JJ)

        ! integer :: ij, ages(JJ)
        ! ! set up age variable
        ! ages = 20 + 5*(/(ij-1, ij=1,JJ)/)
        !
        ! ! polt homeownership ratio
        ! call plot(dble(ages), o_coh(:), legend='Entrepreneurship')
        ! call execplot(xlabel='Age j', ylabel='Entrepreneurship', ylim=(/0d0, 1d0/))
        !
        ! ! plot consumption for homeowner
        ! call plot(dble(ages), c_coh(1, :), legend='Consumption  - Entrepreneur')
        ! call plot(dble(ages), a_coh(1, :), legend='Assets       - Entrepreneur')
        ! call plot(dble(ages), x_coh(1, :), legend='Annuities    - Entrepreneur')
        ! call plot(dble(ages), y_coh(1, :), legend='Income       - Entrepreneur')
        ! call plot(dble(ages), l_coh(1, :), legend='Labor        - Entrepreneur')
        ! call plot(dble(ages), k_coh(:),    legend='Investment   - Entrepreneur')
        ! call execplot(xlabel='Age j', ylabel='Consumption/Assets')
        !
        ! ! polt consumption for renter
        ! call plot(dble(ages), c_coh(0, :), legend='Consumption  - Worker')
        ! call plot(dble(ages), a_coh(0, :), legend='Assets       - Worker')
        ! call plot(dble(ages), x_coh(0, :), legend='Annuities    - Worker')
        ! call plot(dble(ages), y_coh(0, :), legend='Labor Income - Worker')
        ! call plot(dble(ages), l_coh(0, :), legend='Labor        - Worker')
        ! call execplot(xlabel='Age j', ylabel='Consumption/Assets')

        life_exp = 0d0
        do is = 1, NS
          punb(is, 1) = psi(is, 1)
          life_exp(is) = life_exp(is) + 22d0*punb(is, 1)*(1d0-psi(is, 2))
          do ij = 2, JJ
            punb(is, ij) = punb(is, ij-1)*psi(is, ij)
            life_exp(is) = life_exp(is) + (22d0 + 5d0*dble(ij-1))*punb(is, ij)*(1d0-psi(is, ij+1))
          enddo ! ij
        enddo ! is

        write(*,'(/, a, /)')     '******* CALIBRATION *******'
        write(*,'(a, 3f10.4)')   '- life_exp:            ', life_exp
        write(*,'(a, f10.4)')    '- life_exp(avg):       ', sum(life_exp*dist_skill)
        write(*,'(a, f10.4, /)') '- dep. ratio:          ', sum(m(:, :, :, :, :, :, :, JR:JJ))/sum(m(:, :, :, :, :, :, :, 1:JR-1))*100d0
        write(*,'(a, f10.4)')    '- avg. lab. supply (h):', sum(l(:, :, :, :, :, :, :, 1:JR-1)*m(:, :, :, :, :, :, :, 1:JR-1))/sum(m(:, :, :, :, :, :, :, 1:JR-1))
        write(*,'(a, f10.4)')    '  + corp. sector:      ', sum(l(:, 0, :, :, :, :, :, 1:JR-1)*m(:, 0, :, :, :, :, :, 1:JR-1))/sum(m(:, 0, :, :, :, :, :, 1:JR-1))
        write(*,'(a, f10.4, /)') '  + non-corp. sector:  ', sum(l(:, 1:NK, :, :, :, :, :, 1:JR-1)*m(:, 1:NK, :, :, :, :, :, 1:JR-1))/max(sum(m(:, 1:NK, :, :, :, :, :, 1:JR-1)), 1d-4)
        write(*,'(a, f10.4)')    '- pen. ben. (%):       ', PBEN/YY*100d0
        write(*,'(a, f10.4)')    '- pen. con. rate (%):  ', taup*100d0
        write(*,'(a, f10.4)')    '- gov. expend. (%):    ', (GG+(1d0+r)*BB-(1d0+n_p)*BB)/YY*100d0
        write(*,'(a, f10.4)')    '- tax rev. (%):        ', (TAc+TAw+TAr+TAy)/YY*100d0
        write(*,'(a, f10.4)')    '  + cons. tax (%):     ', TAc/(TAc+TAw+TAr+TAy)*100d0
        write(*,'(a, f10.4)')    '  + inc. tax (%):      ', TAw/(TAc+TAw+TAr+TAy)*100d0
        write(*,'(a, f10.4)')    '  + cap. tax (%):      ', TAr/(TAc+TAw+TAr+TAy)*100d0
        write(*,'(a, f10.4)')    '  + corp. tax (%):     ', TAy/(TAc+TAw+TAr+TAy)*100d0
        write(*,'(a, f10.4)')    '- cons. tax rate (%):  ', tauc*100d0
        write(*,'(a, f10.4)')    '- cap.-output ratio:   ', 5*KK/YY
        write(*,'(a, f10.4)')    '  + corp. sector:      ', 5*KC/YC
        write(*,'(a, f10.4, /)') '  + non-corp. sector:  ', 5*KE/max(YE, 1d-4)
        write(*,'(a, f10.4)')    '- int. rate p.a. (%):  ', ((1d0+r)**0.2d0-1d0)*100d0
        write(*,'(a, f10.4)')    '- bequests (%):        ', BQ/YY*100d0
        write(*,*)

        ! write(*,'(a, f10.4)')'KK:', KK
        ! write(*,'(a, f10.4)')'AA:', AA
        ! write(*,'(a, f10.4)')'LC:', LC
        ! write(*,'(a, f10.4)')'YY:', YY
        ! write(*,'(a, f10.4)')'CC:', CC
        ! write(*,'(a, f10.4)')'II:', II
        ! write(*,'(a, f10.4)')'GG:', GG
        ! write(*,'(a, f10.4)')'BB:', BB
        ! write(*,'(a, f10.4)')'r: ', r
        ! write(*,'(a, f10.4)')'w: ', w

    end subroutine

end program
