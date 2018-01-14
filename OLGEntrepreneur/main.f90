include "toolbox.f90"
include "globals.f90"

program main

  ! load modules
  use globals
  use omp_lib

  implicit none

  ! set government variables
  mu     = 1d0
  lambda = 0d0
  phi    = 0d0

  ! calculate initial equilibrium
  call get_SteadyState()

  ! set reforms
  ! mu(1:TT) = 0d0
  ! lambda (1:TT) = 1d0
  ! phi(1:TT) = 1d0

  !tauk = 0d0

  ! calculate transition path
  call get_transition()

  ! close files
  close(21)

contains


  !#############################################################################
  ! FUNCTION get_SteadyState
  !
  ! computes the initial steady state of the economy
  !#############################################################################
  subroutine get_SteadyState()

    implicit none

    ! initialize remaining variables
    call initialize()

    ! start the clock
    call tick(time)

    ! iterate until value function converges
    do iter = 1, itermax

      ! get new prices
      call get_prices(0)

      ! solve the household problem
      call solve_household(1, 0)

      ! calculate the distribution of households over state space
      call get_distribution(0)

      ! aggregate individual decisions
      call aggregation(0)

      ! determine the government parameters
      call government(0)

      ! check maximum grid points used
      call check_grid(iqmax, iamax, ikmax, ixmax, 0)

      write(*,'(i4,4i5,5f8.2,f16.8)')iter, maxval(iqmax), maxval(iamax), maxval(ikmax), maxval(ixmax),&
                                      (/5d0*KK(0), CC(0), II(0)/)/YY(0)*100d0, &
                                      ((1d0+r(0))**0.2d0-1d0)*100d0, w(0), DIFF(0)/YY(0)*100d0

      ! check for convergence
      if(abs(DIFF(0)/YY(0))*100d0 < tol) exit

    enddo ! iter

    ! stop the clock
    call tock(time)

    ! write output
    call output(0)

    if (iter > itermax) write(*,*)'No Convergence'

  end subroutine


  !#############################################################################
  ! SUBROUTINE get_transition
  !
  ! computes the transition path of the economy
  !#############################################################################
  subroutine get_transition()

    implicit none

    !##### OTHER VARIABLES #####################################################
    integer :: iter, ij, it, itmax
    logical :: check

    ! initialize remaining variables
    call initialize_trn()

    ! start the clock
    call tick(time)

    ! iterate until value function converges
    do iter = 1, itermax

      ! get factor and other prices
      do it = 1, TT
        call get_prices(it)
      enddo

      ! solve the household problem
      if (TT > 1) then
        do ij = JJ, 2, -1
          call solve_household(ij, 1)
        enddo
      endif
      do it = 1, TT
        call solve_household(1, it)
      enddo

      ! calculate the distribution of households over state space
      do it = 1, TT
        call get_distribution(it)
      enddo

      ! aggregate individual decisions
      do it = 1, TT
        call aggregation(it)
      enddo

      ! determine the government parameters
      do it = 1, TT
        call government(it)
      enddo

      ! write screen output
      itmax = maxloc(abs(DIFF(1:TT)/YY(1:TT)), 1)
      write(*,'(i4,5f8.2,f14.8)')iter, (/5d0*KK(TT), CC(TT), II(TT)/)/YY(TT)*100d0, &
         ((1d0+r(TT))**0.2d0-1d0)*100d0, w(TT), DIFF(itmax)/YY(itmax)*100d0

      check = abs(DIFF(itmax)/YY(itmax))*100d0 < tol

      ! check for convergence
      if (check) exit

    enddo ! iter

    ! stop the clock
    call tock(time)

    ! write output
    do it = 1, TT
      call output(it)
    enddo

    if (iter > itermax) write(*,*)'No Convergence'

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

    ! write screen output
    write(*,'(/a/)')'INITIAL EQUILIBRIUM'
    write(*,'(a)')'ITER   IQ   IA   IK   IX     K/Y     C/Y     I/Y       r       w            DIFF'

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

    psi(:, 1:JJ) = 1d0

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
    theta(:, 1) = (/0d0, 0.79d0/)
    theta(:, 2) = (/0d0, 0.79d0/)
    theta(:, 3) = (/0d0, 0.98d0/)
    dist_theta(:, 1) = (/0d0, 1d0/)
    dist_theta(:, 2) = (/0d0, 1d0/)
    dist_theta(:, 3) = (/0d0, 1d0/)
    pi_theta(1, 1, :) = 1d0
    pi_theta(1, 2, :) = 0d0
    pi_theta(2, 1, 1) = 0.005d0
    pi_theta(2, 2, 1) = 1d0-pi_theta(2, 1, 1)
    pi_theta(2, 1, 2) = 0.012d0
    pi_theta(2, 2, 2) = 1d0-pi_theta(2, 1, 2)
    pi_theta(2, 1, 3) = 0.017d0
    pi_theta(2, 2, 3) = 1d0-pi_theta(2, 1, 3)

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

    ! initialize tax rates
    tauc(0) = 0.190d0
    taup(0) = 0.189d0

    ! initial guesses for macro variables
    KC(0) = 3.400d0
    LC(0) = 3.604d0
    BQS(:, 0) = (/4.610d-2, 0.180d0, 0.106d0/)
    BB(0) = 2.964d0
    ybar(0) = 0.555d0

    ! initialize value functions
    V = 1d-13**egam/egam; EV = 1d-13**egam/egam; S = 1d-13**egam/egam

    ! initialize policy functions
    Q_plus = 0d0; a_plus = 0d0; k_plus = 0d0; x_plus = 0d0; p_plus = 0d0
    inctax = 0d0; captax = 0d0; penben = 0d0; pencon = 0d0; c = 0d0; l = 0d0

    ! initialize temporary policy and value functions
    Q_plus_t = 0d0; a_plus_t = 0d0; k_plus_t = 0d0; x_plus_t = 0d0; p_plus_t = 0d0
    inctax_t = 0d0; captax_t = 0d0; penben_t = 0d0; pencon_t = 0d0; c_t = 0d0; l_t = 0d0
    omega_x_t = 0d0; omega_k_t = 0d0
    V_t = 0d0

    ! get initial guess for household decisions
    omega_x_t(:, :, :, :, :, :, :, :, 1:JR-1, 0) = 0.05d0
    omega_k_t(:, :, :, :, :, :, :, :, 1:JR-2, 0) = 0.05d0
    l_t(:, :, :, :, :, :, :, :, 1:JR-1, 0) = 0.33d0
    do ia = 0, NA
      Q_plus_t(:, ia, :, :, :, :, :, :, :, 0) = a(ia)/2d0
    enddo ! ia

    ! open files
    open(21, file='output.out')

  end subroutine


  !#############################################################################
  ! SUBROUTINE initialize_trn
  !
  ! initializes transitional variables
  !#############################################################################
  subroutine initialize_trn()

    implicit none

    !##### OTHER VARIABLES #####################################################
    integer :: it

    write(*,'(/a/)')'TRANSITION PATH'
    write(*,'(a)')'ITER     K/Y     C/Y     I/Y       r       w          DIFF'

    do it = 1, TT

      r(it) = r(0)
      w(it) = w(0)
      ybar(it) = ybar(0)
      pinv(it) = pinv(0)

      AA(it) = AA(0)
      AX(it) = AX(0)
      BQ(it) = BQ(0)
      PBEN(it) = PBEN(0)
      PCON(it) = PCON(0)
      KK(it) = KK(0)
      KC(it) = KC(0)
      KE(it) = KE(0)
      LC(it) = LC(0)
      BB(it) = BB(0)
      YY(it) = YY(0)
      YC(it) = YC(0)
      YE(it) = YE(0)
      CC(it) = CC(0)
      II(it) = II(0)
      TC(it) = TC(0)

      TAc(it) = TAc(0)
      TAr(it) = TAr(0)
      TAw(it) = TAw(0)
      TAk(it) = TAk(0)

      BQS(:, it) = BQS(:, 0)

      tauc(it) = tauc(0)
      taup(it) = taup(0)

      ann(:, :, :, it) = ann(:, :, :, 0)
      pen(:, :, it) = pen(:, :, 0)
      beq(:, :, it) = beq(:, :, 0)

      Q_plus(:, :, :, :, :, :, :, :, it) = Q_plus(:, :, :, :, :, :, :, :, 0)
      a_plus(:, :, :, :, :, :, :, :, it) = a_plus(:, :, :, :, :, :, :, :, 0)
      k_plus(:, :, :, :, :, :, :, :, it) = k_plus(:, :, :, :, :, :, :, :, 0)
      p_plus(:, :, :, :, :, :, :, :, it) = p_plus(:, :, :, :, :, :, :, :, 0)
      inctax(:, :, :, :, :, :, :, :, it) = inctax(:, :, :, :, :, :, :, :, 0)
      captax(:, :, :, :, :, :, :, :, it) = captax(:, :, :, :, :, :, :, :, 0)
      penben(:, :, :, :, :, :, :, :, it) = penben(:, :, :, :, :, :, :, :, 0)
      pencon(:, :, :, :, :, :, :, :, it) = pencon(:, :, :, :, :, :, :, :, 0)
      c(:, :, :, :, :, :, :, :, it) = c(:, :, :, :, :, :, :, :, 0)
      l(:, :, :, :, :, :, :, :, it) = l(:, :, :, :, :, :, :, :, 0)
      V(:, :, :, :, :, :, :, :, it) = V(:, :, :, :, :, :, :, :, 0)
      EV(:, :, :, :, :, :, :, :, it) = EV(:, :, :, :, :, :, :, :, 0)

      Q_plus_t(:, :, :, :, :, :, :, :, :, it) = Q_plus_t(:, :, :, :, :, :, :, :, :, 0)
      a_plus_t(:, :, :, :, :, :, :, :, :, it) = a_plus_t(:, :, :, :, :, :, :, :, :, 0)
      k_plus_t(:, :, :, :, :, :, :, :, :, it) = k_plus_t(:, :, :, :, :, :, :, :, :, 0)
      p_plus_t(:, :, :, :, :, :, :, :, :, it) = p_plus_t(:, :, :, :, :, :, :, :, :, 0)
      inctax_t(:, :, :, :, :, :, :, :, :, it) = inctax_t(:, :, :, :, :, :, :, :, :, 0)
      captax_t(:, :, :, :, :, :, :, :, :, it) = captax_t(:, :, :, :, :, :, :, :, :, 0)
      penben_t(:, :, :, :, :, :, :, :, :, it) = penben_t(:, :, :, :, :, :, :, :, :, 0)
      pencon_t(:, :, :, :, :, :, :, :, :, it) = pencon_t(:, :, :, :, :, :, :, :, :, 0)
      c_t(:, :, :, :, :, :, :, :, :, it) = c_t(:, :, :, :, :, :, :, :, :, 0)
      l_t(:, :, :, :, :, :, :, :, :, it) = l_t(:, :, :, :, :, :, :, :, :, 0)
      V_t(:, :, :, :, :, :, :, :, :, it) = V_t(:, :, :, :, :, :, :, :, :, 0)

      omega_x_t(:, :, :, :, :, :, :, :, :, it) = omega_x_t(:, :, :, :, :, :, :, :, :, 0)
      omega_k_t(:, :, :, :, :, :, :, :, :, it) = omega_k_t(:, :, :, :, :, :, :, :, :, 0)
      S(:, :, :, :, :, :, :, :, :, it) = S(:, :, :, :, :, :, :, :, :, 0)

      m(:, :, :, :, :, :, :, :, it) = m(:, :, :, :, :, :, :, :, 0)
      m_Q(:, :, :, :, :, :, :, :, it) = m_Q(:, :, :, :, :, :, :, :, 0)

    enddo

  end subroutine


  !#############################################################################
  ! SUBROUTINE get_prices
  !
  ! computes prices and distribution of bequests for next iteration step
  !#############################################################################
  subroutine get_prices(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES #####################################################
    real*8 :: ann_tmp(NS)
    integer :: ix, ip, is, ij, iij, itm, itp

    ! calculate new prices
    r(it) = (1d0-tauk)*(Omega*alpha*(KC(it)/LC(it))**(alpha-1d0)-delta_k)
    w(it) = Omega*(1d0-alpha)*(KC(it)/LC(it))**alpha

    ! set prices in case of life-cycle model
    ! r = 0.393280506035032d0
    ! w = 0.877841532937879d0
    ! BQS = (/4.608543623547606d-2, 0.181029882698876d0, 0.106845332164835d0/)
    ! ybar = 0.555719715351030d0
    ! tauc = 0.128579256047982d0
    ! taup = 7.867802841513299d-2

    ! calculate gross price of consumption (inverse)
    pinv(it) = 1d0/(1d0+tauc(it))

    ! calculate individual bequests
    beq(1, :, it) = Gama(:)*BQS(1, it)/rpop(1, :)
    beq(2, :, it) = Gama(:)*BQS(2, it)/rpop(2, :)
    beq(3, :, it) = Gama(:)*BQS(3, it)/rpop(3, :)

    ! determine the income tax system
    r1 = 0.278d0*ybar(0)*2d0 !  8,354.00 Euro
    r2 = 0.449d0*ybar(0)*2d0 ! 13,469.00 Euro
    r3 = 1.763d0*ybar(0)*2d0 ! 52,881.00 Euro

    ! calculate annuity payments
    ann(:, :, :, it) = 0d0

    do ij = JR, JJ
      ann_tmp = 1d0
      do iij = JJ, ij+1, -1
        itp = year(it, ij, iij)
        ann_tmp(:) = ann_tmp(:)/(1d0+r(itp))*psi(:, iij) + 1d0
      enddo ! iij
      do ix = 0, NX
        ann(ix, :, ij, it) = (1d0+r(it))/psi(:, ij)*x(ix)/ann_tmp(:)
      enddo ! ix
    enddo ! ij6

    ! calculate old-age transfers
    pen(:, :, it) = 0d0
    do ip = 0, NP
      pen(ip, JR:JJ, it) = p(ip)*kappa*ybar(it)
    enddo

  end subroutine


  !#############################################################################
  ! FUNCTION solve_household
  !
  ! determines the solution to the household optimization problem
  !#############################################################################
  subroutine solve_household(ij_in, it_in)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it_in, ij_in

    !##### OTHER VARIABLES #####################################################
    integer :: iq, ia, ik, ix, ip, iw, ie, is, ij, it, iq_p, ip_p, io_p
    integer :: ij_max

    if (TT>1) ij_max = min(TT-it_in+1, JJ)
    ! ij_max = JJ

    ! solve household problem recursively
    do ij = ij_max, 1, -1

      it = year(it_in, ij_in, ij)

      ! solve for oldest cohort
      if (ij == JJ) then

        omega_x_t(:, :, :, :, :, :, :, :, ij, it) = 0d0
        omega_k_t(:, :, :, :, :, :, :, :, ij, it) = 0d0

        do iq_p = 0, NQ
            S(:, iq_p, :, :, :, :, :, :, ij, it) = mu_b*max(Q(iq_p), 1d-13)**egam/egam
        enddo ! iq_p

        !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads) shared(ij, it)
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NX
              do ia = 0, NA

                call solve_consumption(0, ia, 0, ix, ip, 1, 1, is, ij, it)

                ! copy decisions
                Q_plus(ia, :, ix, ip, :, :, is, ij, it) = Q_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                a_plus(ia, :, ix, ip, :, :, is, ij, it) = a_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                x_plus(ia, :, ix, ip, :, :, is, ij, it) = x_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                p_plus(ia, :, ix, ip, :, :, is, ij, it) = p_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                k_plus(ia, :, ix, ip, :, :, is, ij, it) = k_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                c(ia, :, ix, ip, :, :, is, ij, it) = c_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                l(ia, :, ix, ip, :, :, is, ij, it) = l_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                inctax(ia, :, ix, ip, :, :, is, ij, it) = inctax_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                captax(ia, :, ix, ip, :, :, is, ij, it) = captax_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                penben(ia, :, ix, ip, :, :, is, ij, it) = penben_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                pencon(ia, :, ix, ip, :, :, is, ij, it) = pencon_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                V(ia, :, ix, ip, :, :, is, ij, it) = V_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)

              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
        !$omp end parallel do

      ! solve for retirement age
      elseif (ij >= JR) then

        !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads) shared(ij, it)
        do is = 1, NS
          do ip_p = 0, NP
            do ix = 0, NX
              do iq_p = 0, NQ

                ! next period retiree
                call solve_retiree(iq_p, 0, ix, ip_p, 1, 1, is, ij, it)

                omega_x_t(:, iq_p, :, ix, ip_p, :, :, is, ij, it) = omega_x_t(0, iq_p, 0, ix, ip_p, 1, 1, is, ij, it)
                omega_k_t(:, iq_p, :, ix, ip_p, :, :, is, ij, it) = omega_k_t(0, iq_p, 0, ix, ip_p, 1, 1, is, ij, it)
                S(:, iq_p, :, ix, ip_p, :, :, is, ij, it) = S(0, iq_p, 0, ix, ip_p, 1, 1, is, ij, it)

              enddo ! iq_p
            enddo ! ix
          enddo ! ip_p
        enddo ! is
        !$omp end parallel do

        !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads) shared(ij, it)
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NX
              do ia = 0, NA

                ! next period worker
                call solve_consumption(0, ia, 0, ix, ip, 1, 1, is, ij, it)

                ! copy decisions
                Q_plus(ia, :, ix, ip, :, :, is, ij, it) = Q_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                a_plus(ia, :, ix, ip, :, :, is, ij, it) = a_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                x_plus(ia, :, ix, ip, :, :, is, ij, it) = x_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                p_plus(ia, :, ix, ip, :, :, is, ij, it) = p_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                k_plus(ia, :, ix, ip, :, :, is, ij, it) = k_plus_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                c(ia, :, ix, ip, :, :, is, ij, it) = c_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                l(ia, :, ix, ip, :, :, is, ij, it) = l_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                inctax(ia, :, ix, ip, :, :, is, ij, it) = inctax_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                captax(ia, :, ix, ip, :, :, is, ij, it) = captax_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                penben(ia, :, ix, ip, :, :, is, ij, it) = penben_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                pencon(ia, :, ix, ip, :, :, is, ij, it) = pencon_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)
                V(ia, :, ix, ip, :, :, is, ij, it) = V_t(0, ia, 0, ix, ip, 1, 1, is, ij, it)

              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
        !$omp end parallel do

      ! solve for working age
      elseif (ij >= 2) then

        !$omp parallel do collapse(4) schedule(dynamic) num_threads(numthreads) shared(ij, it)
        do is = 1, NS
          do ie = 1, NE
            do iw = 1, NW
              do ip_p = 0, NP
                do ix = 0, NX
                  do ik = 0, NK
                    do iq_p = 0, NQ

                      ! next period worker
                      call solve_worker(iq_p, ik, ix, ip_p, iw, ie, is, ij, it)

                      ! next period entrepreneur
                      call solve_entrepreneur(iq_p, ik, ix, ip_p, iw, ie, is, ij, it)

                    enddo ! iq_p
                  enddo ! ik
                enddo ! ix
              enddo ! ip_p
            enddo ! iw
          enddo ! ie
        enddo ! is
        !$omp end parallel do

        !$omp parallel do collapse(4) schedule(dynamic) num_threads(numthreads) shared(ij, it)
        do is = 1, NS
          do ie = 1, NE
            do iw = 1, NW
              do ip = 0, NP
                do ix = 0, NX
                  do ik = 0, NK
                    do ia = 0, NA

                      ! next period worker
                      call solve_consumption(0, ia, ik, ix, ip, iw, ie, is, ij, it)

                      ! next period entrpreneur
                      if (ij < JR-1) call solve_consumption(1, ia, ik, ix, ip, iw, ie, is, ij, it)

                      ! decision on whether to be homeowner or renter next period
                      io_p = 0
                      if (ij < JR-1 .and. V_t(1, ia, ik, ix, ip, iw, ie, is, ij, it) > V_t(0, ia, ik, ix, ip, iw, ie, is, ij, it)) io_p = 1

                      ! copy decisions
                      Q_plus(ia, ik, ix, ip, iw, ie, is, ij, it) = Q_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      a_plus(ia, ik, ix, ip, iw, ie, is, ij, it) = a_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      x_plus(ia, ik, ix, ip, iw, ie, is, ij, it) = x_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      p_plus(ia, ik, ix, ip, iw, ie, is, ij, it) = p_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      k_plus(ia, ik, ix, ip, iw, ie, is, ij, it) = k_plus_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      c(ia, ik, ix, ip, iw, ie, is, ij, it) = c_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      l(ia, ik, ix, ip, iw, ie, is, ij, it) = l_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      inctax(ia, ik, ix, ip, iw, ie, is, ij, it) = inctax_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      captax(ia, ik, ix, ip, iw, ie, is, ij, it) = captax_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      penben(ia, ik, ix, ip, iw, ie, is, ij, it) = penben_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      pencon(ia, ik, ix, ip, iw, ie, is, ij, it) = pencon_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)
                      V(ia, ik, ix, ip, iw, ie, is, ij, it) = V_t(io_p, ia, ik, ix, ip, iw, ie, is, ij, it)

                    enddo ! ia
                  enddo ! ik
                enddo ! ix
              enddo ! ip
            enddo ! iw
          enddo ! ie
        enddo ! is
        !$omp end parallel do

      ! solve for youngest cohort
      elseif (ij == 1) then

        !$omp parallel do collapse(4) schedule(dynamic) num_threads(numthreads) shared(ij, it)
        do is = 1, NS
          do ie = 1, NE
            do iw = 1, NW
              do ip_p = 0, NP
                do iq_p = 0, NQ

                  ! next period worker
                  call solve_worker(iq_p, 0, 0, ip_p, iw, ie, is, ij, it)

                  ! copy decisions
                  omega_x_t(0, iq_p, :, :, ip_p, iw, ie, is, ij, it) = omega_x_t(0, iq_p, 0, 0, ip_p, iw, ie, is, ij, it)
                  omega_k_t(0, iq_p, :, :, ip_p, iw, ie, is, ij, it) = omega_k_t(0, iq_p, 0, 0, ip_p, iw, ie, is, ij, it)
                  S(0, iq_p, :, :, ip_p, iw, ie, is, ij, it) = S(0, iq_p, 0, 0, ip_p, iw, ie, is, ij, it)

                  ! next period entrepreneur
                  call solve_entrepreneur(iq_p, 0, 0, ip_p, iw, ie, is, ij, it)

                  ! copy decisions
                  omega_x_t(1, iq_p, :, :, ip_p, iw, ie, is, ij, it) = omega_x_t(1, iq_p, 0, 0, ip_p, iw, ie, is, ij, it)
                  omega_k_t(1, iq_p, :, :, ip_p, iw, ie, is, ij, it) = omega_k_t(1, iq_p, 0, 0, ip_p, iw, ie, is, ij, it)
                  S(1, iq_p, :, :, ip_p, iw, ie, is, ij, it) = S(1, iq_p, 0, 0, ip_p, iw, ie, is, ij, it)

                enddo ! iq_p
              enddo ! ip_p
            enddo ! iw
          enddo ! ie
        enddo ! is
        !$omp end parallel do

        !$omp parallel do collapse(3) schedule(dynamic) num_threads(numthreads) shared(ij, it)
        ! solve the consumption savings problem
        do is = 1, NS
          do ie = 1, NE
            do iw = 1, NW

              ! next period worker
              call solve_consumption(0, 0, 0, 0, 0, iw, ie, is, ij, it)

              ! next period entrpreneur
              if(ij < JR-1) call solve_consumption(1, 0, 0, 0, 0, iw, ie, is, ij, it)

              ! decision on whether to be homeowner or renter next period
              io_p = 0
              if(ij < JR-1 .and. V_t(1, 0, 0, 0, 0, iw, ie, is, ij, it) > V_t(0, 0, 0, 0, 0, iw, ie, is, ij, it)) io_p = 1

              ! copy decisions
              Q_plus(:, :, :, :, iw, ie, is, ij, it) = Q_plus_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              a_plus(:, :, :, :, iw, ie, is, ij, it) = a_plus_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              x_plus(:, :, :, :, iw, ie, is, ij, it) = x_plus_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              p_plus(:, :, :, :, iw, ie, is, ij, it) = p_plus_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              k_plus(:, :, :, :, iw, ie, is, ij, it) = k_plus_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              c(:, :, :, :, iw, ie, is, ij, it) = c_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              l(:, :, :, :, iw, ie, is, ij, it) = l_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              inctax(:, :, :, :, iw, ie, is, ij, it) = inctax_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              captax(:, :, :, :, iw, ie, is, ij, it) = captax_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              penben(:, :, :, :, iw, ie, is, ij, it) = penben_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              pencon(:, :, :, :, iw, ie, is, ij, it) = pencon_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)
              V(:, :, :, :, iw, ie, is, ij, it) = V_t(io_p, 0, 0, 0, 0, iw, ie, is, ij, it)

            enddo ! iw
          enddo ! ie
        enddo ! is
        !$omp end parallel do

      endif

      call interpolate(ij, it)
      !write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

    enddo ! ij

  end subroutine


  !#############################################################################
  ! SUBROUTINE interpolate
  !
  ! calculates the expected valuefunction of cohort ij
  !#############################################################################
  subroutine interpolate(ij, it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: ij, it

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

                  EV(ia, ik, ix, ip, iw, ie, is, ij, it) = 0d0
                  do ie_p = 1, NE
                    do iw_p = 1, NW
                      EV(ia, ik, ix, ip, iw, ie, is, ij, it) = EV(ia, ik, ix, ip, iw, ie, is, ij, it) &
                        + pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*V(ia, ik, ix, ip, iw, ie, is, ij, it)
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
  subroutine get_distribution(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES #####################################################
    integer :: ia, ik, ix, ip, iw, ie, is, ij, itm, iw_p, ie_p
    integer :: iql, iqr, ial, iar, ikl, ikr, ixl, ixr, ipl, ipr
    real*8 :: varphi_q, varphi_a, varphi_k, varphi_x, varphi_p

    ! get last year
    itm = year(it, 2, 1)

    m(:, :, :, :, :, :, :, :, it) = 0d0
    m_Q(:, :, :, :, :, :, :, :, it) = 0d0

    do is = 1, NS
      do iw = 1, NW
        do ie = 1, NE
            m(0, 0, 0, 0, iw, ie, is, 1, it) = dist_eta(iw, is)*dist_theta(ie, is)*dist_skill(is)
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
                    !if (m(ia, ik, ix, ip, iw, ie, is, ij-1, itm) <= 0d0) cycle

                    ! derive interpolation weights
                    call linint_Grow(Q_plus(ia, ik, ix, ip, iw, ie, is, ij-1, itm), Q_l, Q_u, Q_grow, NQ, iql, iqr, varphi_q)
                    call linint_Grow(a_plus(ia, ik, ix, ip, iw, ie, is, ij-1, itm), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                    if (NK > 0) then
                      call linint_Grow(k_plus(ia, ik, ix, ip, iw, ie, is, ij-1, itm), k_l, k_u, k_grow, NK-1, ikl, ikr, varphi_k)
                    else
                      ikl = 0; ikr = 0; varphi_k = 1d0
                    endif
                    if (NX > 0) then
                      call linint_Grow(x_plus(ia, ik, ix, ip, iw, ie, is, ij-1, itm), x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)
                    else
                      ixl = 0; ixr = 0; varphi_x = 1d0
                    endif
                    call linint_Equi(p_plus(ia, ik, ix, ip, iw, ie, is, ij-1, itm), p_l, p_u, NP, ipl, ipr, varphi_p)

                    ! restrict values to grid just in case
                    iql = min(iql, NQ)
                    iqr = min(iqr, NQ)
                    varphi_q = max(min(varphi_q, 1d0),0d0)

                    ! restrict values to grid just in case
                    ial = min(ial, NA)
                    iar = min(iar, NA)
                    varphi_a = max(min(varphi_a, 1d0),0d0)

                    ! restrict values to grid just in case
                    if (k_plus(ia, ik, ix, ip, iw, ie, is, ij-1, itm) >= k_min) then
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

                        m(ial, ikl, ixl, ipl, iw_p, ie_p, is, ij, it) = m(ial, ikl, ixl, ipl, iw_p, ie_p, is, ij, it) + &
                              varphi_a*varphi_k*varphi_x*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(ial, ikl, ixl, ipr, iw_p, ie_p, is, ij, it) = m(ial, ikl, ixl, ipr, iw_p, ie_p, is, ij, it) + &
                              varphi_a*varphi_k*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(ial, ikl, ixr, ipl, iw_p, ie_p, is, ij, it) = m(ial, ikl, ixr, ipl, iw_p, ie_p, is, ij, it) + &
                              varphi_a*varphi_k*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(ial, ikl, ixr, ipr, iw_p, ie_p, is, ij, it) = m(ial, ikl, ixr, ipr, iw_p, ie_p, is, ij, it) + &
                              varphi_a*varphi_k*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(ial, ikr, ixl, ipl, iw_p, ie_p, is, ij, it) = m(ial, ikr, ixl, ipl, iw_p, ie_p, is, ij, it) + &
                              varphi_a*(1d0-varphi_k)*varphi_x*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(ial, ikr, ixl, ipr, iw_p, ie_p, is, ij, it) = m(ial, ikr, ixl, ipr, iw_p, ie_p, is, ij, it) + &
                              varphi_a*(1d0-varphi_k)*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(ial, ikr, ixr, ipl, iw_p, ie_p, is, ij, it) = m(ial, ikr, ixr, ipl, iw_p, ie_p, is, ij, it) + &
                              varphi_a*(1d0-varphi_k)*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(ial, ikr, ixr, ipr, iw_p, ie_p, is, ij, it) = m(ial, ikr, ixr, ipr, iw_p, ie_p, is, ij, it) + &
                              varphi_a*(1d0-varphi_k)*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(iar, ikl, ixl, ipl, iw_p, ie_p, is, ij, it) = m(iar, ikl, ixl, ipl, iw_p, ie_p, is, ij, it) + &
                              (1d0-varphi_a)*varphi_k*varphi_x*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(iar, ikl, ixl, ipr, iw_p, ie_p, is, ij, it) = m(iar, ikl, ixl, ipr, iw_p, ie_p, is, ij, it) + &
                              (1d0-varphi_a)*varphi_k*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(iar, ikl, ixr, ipl, iw_p, ie_p, is, ij, it) = m(iar, ikl, ixr, ipl, iw_p, ie_p, is, ij, it) + &
                              (1d0-varphi_a)*varphi_k*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(iar, ikl, ixr, ipr, iw_p, ie_p, is, ij, it) = m(iar, ikl, ixr, ipr, iw_p, ie_p, is, ij, it) + &
                              (1d0-varphi_a)*varphi_k*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(iar, ikr, ixl, ipl, iw_p, ie_p, is, ij, it) = m(iar, ikr, ixl, ipl, iw_p, ie_p, is, ij, it) + &
                              (1d0-varphi_a)*(1d0-varphi_k)*varphi_x*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(iar, ikr, ixl, ipr, iw_p, ie_p, is, ij, it) = m(iar, ikr, ixl, ipr, iw_p, ie_p, is, ij, it) + &
                              (1d0-varphi_a)*(1d0-varphi_k)*varphi_x*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(iar, ikr, ixr, ipl, iw_p, ie_p, is, ij, it) = m(iar, ikr, ixr, ipl, iw_p, ie_p, is, ij, it) + &
                              (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*varphi_p*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(iar, ikr, ixr, ipr, iw_p, ie_p, is, ij, it) = m(iar, ikr, ixr, ipr, iw_p, ie_p, is, ij, it) + &
                              (1d0-varphi_a)*(1d0-varphi_k)*(1d0-varphi_x)*(1d0-varphi_p)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)

                      enddo
                    enddo

                    m_Q(iql, ik, ix, ipl, iw, ie, is, ij, it) = m_Q(iql, ik, ix, ipl, iw, ie, is, ij, it) + &
                                varphi_q*varphi_p*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)
                    m_Q(iql, ik, ix, ipr, iw, ie, is, ij, it) = m_Q(iql, ik, ix, ipr, iw, ie, is, ij, it) + &
                                varphi_q*(1d0-varphi_p)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)
                    m_Q(iqr, ik, ix, ipl, iw, ie, is, ij, it) = m_Q(iqr, ik, ix, ipl, iw, ie, is, ij, it) + &
                                (1d0-varphi_q)*varphi_p*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)
                    m_Q(iqr, ik, ix, ipr, iw, ie, is, ij, it) = m_Q(iqr, ik, ix, ipr, iw, ie, is, ij, it) + &
                                (1d0-varphi_q)*(1d0-varphi_p)*m(ia, ik, ix, ip, iw, ie, is, ij-1, itm)

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
  subroutine aggregation(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES #####################################################
    integer :: ia, ik, ix, ip, iw, ie, is, ij, itm, itp
    real*8 :: LC_old

    ! get last and next year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    ! copy labor supply
    LC_old = LC(it)

    ! reset macroeconomic aggregates in each iteration step
    AA(it) = 0d0; AX(it) = 0d0; BQ(it) = 0d0; PBEN(it) = 0d0; PCON(it) = 0d0
    KK(it) = 0d0; KE(it) = 0d0; LC(it) = 0d0;
    YY(it) = 0d0; YC(it) = 0d0; YE(it) = 0d0; CC(it) = 0d0;  II(it) = 0d0; TC(it) = 0d0
    TAc(it) = 0d0; TAr(it) = 0d0; TAw(it) = 0d0; TAk(it) = 0d0
    BQS(:, it) = 0d0;

    do ij = 1, JJ

      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ik = 0, NK
                  do ia = 0, NA

                    ! skip if there is no household
                    !if (m(ia, ik, ix, ip, iw, ie, is, ij, it) <= 0d0 .and. m(ia, ik, ix, ip, iw, ie, is, ij, itm) <= 0d0) cycle

                    ! AA(it) = AA(it) + (a_plus(ia, ik, ix, ip, iw, ie, is, ij, itm)-xi*k_plus(ia, ik, ix, ip, iw, ie, is, ij, itm))*psi(is, ij+1)*m(ia, ik, ix, ip, iw, ie, is, ij, itm)/(1d0+n_p)
                    AA(it) = AA(it) + a_plus(ia, ik, ix, ip, iw, ie, is, ij, itm)*m(ia, ik, ix, ip, iw, ie, is, ij, itm)
                    AX(it) = AX(it) + x(ix)/psi(is, ij)*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    CC(it) = CC(it) + c(ia, ik, ix, ip, iw, ie, is, ij, it)*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    ! BQS(is, it) = BQS(is, it) + (a_plus(ia, ik, ix, ip, iw, ie, is, ij, itm)+(1d0-xi)*k_plus(ia, ik, ix, ip, iw, ie, is, ij, itm))*(1d0-psi(is, ij+1))*m(ia, ik, ix, ip, iw, ie, is, ij, itm)
                    BQS(is, it) = BQS(is, it) + (1d0+r(it))*a_plus(ia, ik, ix, ip, iw, ie, is, ij, itm)*(1d0-psi(is, ij+1))*m(ia, ik, ix, ip, iw, ie, is, ij, itm)
                    KE(it) = KE(it) + k(ik)*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    TC(it) = TC(it) + tr(k(ik), k_plus(ia, ik, ix, ip, iw, ie, is, ij, it))*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    TAc(it) = TAc(it) + tauc(it)*c(ia, ik, ix, ip, iw, ie, is, ij, it)*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    TAw(it) = TAw(it) + inctax(ia, ik, ix, ip, iw, ie, is, ij, it)*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    TAr(it) = TAr(it) + captax(ia, ik, ix, ip, iw, ie, is, ij, it)*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    PBEN(it) = PBEN(it) + penben(ia, ik, ix, ip, iw, ie, is, ij, it)*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    PCON(it) = PCON(it) + pencon(ia, ik, ix, ip, iw, ie, is, ij, it)*m(ia, ik, ix, ip, iw, ie, is, ij, it)

                    if(ik == 0) then
                      LC(it) = LC(it) + eff(is, ij)*eta(iw, is)*l(ia, ik, ix, ip, iw, ie, is, ij, it)*m(ia, ik, ix, ip, iw, ie, is, ij, it)
                    else
                      YE(it) = YE(it) + theta(ie, is)*k(ik)**nu1*(eff(is, ij)*l(ia, ik, ix, ip, iw, ie, is, ij, it))**nu2*m(ia, ik, ix, ip, iw, ie, is, ij, it)
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
    ybar(it) = w(it)*LC(it)/sum(m(:, :, :, :, :, :, :, 1:JR-1, it))

    ! compute stock of capital
    KC(it) = damp*(AA(it) + AX(it) - BB(it)) + (1d0-damp)*KC(it)
    KK(it) = KC(it) + KE(it)

    ! update work supply
    LC(it) = damp*LC(it) + (1d0-damp)*LC_old

    ! compute total bequests
    BQ(it) = sum(BQS(:, it))

    ! commpute investment
    II(it) = (1d0+n_p)*KK(itp) - (1d0-delta_k)*KK(it)

    ! compute output
    YC(it) = Omega*KC(it)**alpha*LC(it)**(1d0-alpha)
    YY(it) = YC(it) + YE(it)

    ! compute corporate tax incom
    TAk(it) = tauk*(YC(it)-delta_k*KC(it)-w(it)*LC(it))

  end subroutine


  !#############################################################################
  ! SUBROUTINE government
  !
  ! calculates government parameters
  !#############################################################################
  subroutine government(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES ##############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES #####################################################
    integer :: itp
    real*8 :: expend

    ! get next year
    itp = year(it, 1, 2)

    ! computes government expenditures
    GG(it) = gy*YY(0)
    BB(it) = by*YY(0)
    expend = GG(it) + (1d0+r(it))*BB(it) - (1d0+n_p)*BB(itp)

    ! calculates consumption tax rate
    tauc(it) = (expend - TAk(it) - TAw(it) - TAr(it))/CC(it)

    ! get budget balancing pension contribution rate
    taup(it) = PBEN(it)/PCON(it)

    ! compute gap on goods market
    DIFF(it) = YY(it)-CC(it)-II(it)-TC(it)-GG(it)

  end subroutine


  ! subroutine for writing output
  subroutine output(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it

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

    if (it == 0) then

      write(*,'(/, a, /)')     '******* CALIBRATION *******'
      write(*,'(a, 3f10.4)')   '- life_exp:            ', life_exp
      write(*,'(a, f10.4)')    '  + (total):           ', sum(life_exp*dist_skill)
      write(*,'(a, f10.4, /)') '- dep. ratio:          ', sum(m(:, :, :, :, :, :, :, JR:JJ, it))/sum(m(:, :, :, :, :, :, :, 1:JR-1, it))*100d0
      write(*,'(a, 3f10.4)')   '- fraction of ent. (%):', sum(m(:, 1:NK, :, :, :, :, 1, 1:JR-1, it))/sum(m(:, :, :, :, :, :, 1, 1:JR-1, it))*100d0, sum(m(:, 1:NK, :, :, :, :, 2, 1:JR-1, it))/sum(m(:, :, :, :, :, :, 2, 1:JR-1, it))*100d0, sum(m(:, 1:NK, :, :, :, :, 3, 1:JR-1, it))/sum(m(:, :, :, :, :, :, 3, 1:JR-1, it))*100d0
      write(*,'(a, f10.4)')    '  + (total):           ', sum(m(:, 1:NK, :, :, :, :, :, 1:JR-1, it))/sum(m(:, :, :, :, :, :, :, 1:JR-1, it))*100d0
      write(*,'(a, f10.4)')    '- avg. lab. supply (h):', sum(l(:, :, :, :, :, :, :, 1:JR-1, it)*m(:, :, :, :, :, :, :, 1:JR-1, it))/sum(m(:, :, :, :, :, :, :, 1:JR-1, it))
      write(*,'(a, f10.4)')    '  + corp. sector:      ', sum(l(:, 0, :, :, :, :, :, 1:JR-1, it)*m(:, 0, :, :, :, :, :, 1:JR-1, it))/sum(m(:, 0, :, :, :, :, :, 1:JR-1, it))
      write(*,'(a, f10.4, /)') '  + non-corp. sector:  ', sum(l(:, 1:NK, :, :, :, :, :, 1:JR-1, it)*m(:, 1:NK, :, :, :, :, :, 1:JR-1, it))/max(sum(m(:, 1:NK, :, :, :, :, :, 1:JR-1, it)), 1d-4)
      write(*,'(a, f10.4)')    '- pen. ben. (%):       ', PBEN(it)/YY(it)*100d0
      write(*,'(a, f10.4)')    '- pen. con. rate (%):  ', taup(it)*100d0
      write(*,'(a, f10.4)')    '- tax rev. (%):        ', (TAc(it)+TAw(it)+TAr(it)+TAk(it))/YY(it)*100d0
      write(*,'(a, f10.4)')    '  + cons. tax (%):     ', TAc(it)/(TAc(it)+TAw(it)+TAr(it)+TAk(it))*100d0
      write(*,'(a, f10.4)')    '  + inc. tax (%):      ', TAw(it)/(TAc(it)+TAw(it)+TAr(it)+TAk(it))*100d0
      write(*,'(a, f10.4)')    '  + cap. tax (%):      ', TAr(it)/(TAc(it)+TAw(it)+TAr(it)+TAk(it))*100d0
      write(*,'(a, f10.4)')    '  + corp. tax (%):     ', TAk(it)/(TAc(it)+TAw(it)+TAr(it)+TAk(it))*100d0
      write(*,'(a, f10.4)')    '- cons. tax rate (%):  ', tauc(it)*100d0
      write(*,'(a, f10.4)')    '- cap.-output ratio:   ', 5d0*KK(it)/YY(it)
      write(*,'(a, f10.4)')    '  + corp. sector:      ', 5d0*KC(it)/YC(it)
      write(*,'(a, f10.4, /)') '  + non-corp. sector:  ', 5d0*KE(it)/max(YE(it), 1d-4)
      write(*,'(a, f10.4)')    '- int. rate p.a. (%):  ', ((1d0+r(it))**0.2d0-1d0)*100d0
      write(*,'(a, f10.4)')    '- bequests (%):        ', BQ(it)/YY(it)*100d0
      write(*,*)

    endif

  end subroutine

end program
