program main

  ! modules
  use globals
  use omp_lib
  use clock
  use toolbox

  implicit none

  integer, parameter :: numthreads = 56

  ! allocate arrays
  if(allocated(aplus))deallocate(aplus)
  if(allocated(xplus))deallocate(xplus)
  if(allocated(pplus))deallocate(pplus)
  if(allocated(c))deallocate(c)
  if(allocated(l))deallocate(l)
  if(allocated(k))deallocate(k)
  if(allocated(mx))deallocate(mx)
  if(allocated(oplus))deallocate(oplus)
  if(allocated(pencon))deallocate(pencon)
  if(allocated(inctax))deallocate(inctax)
  if(allocated(captax))deallocate(captax)
  if(allocated(m))deallocate(m)
  if(allocated(VV))deallocate(VV)
  if(allocated(EV))deallocate(EV)
  if(allocated(v))deallocate(v)
  allocate(aplus(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(xplus(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(pplus(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(c(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(l(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(mx(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(k(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(oplus(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(pencon(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(inctax(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(captax(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(m(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(VV(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(EV(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))
  allocate(v(0:1, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ, 0:TT))

  ! set compensating payments to zero
  v = 0d0

  ! household preference parameters
  gamma  =  0.500d0
  ! invert gamma
  gamma = 1d0/gamma
  sigma  =  0.320d0
  phi1   = -11.600d0
  phi2   =  11.600d0
  ! invert phi2
  phi2 = 1d0/phi2
  sigmaq =  1.500d0
  beta   =  0.988d0
  ! convert variables into per period values
  beta = beta**5d0

  ! production parameters
  alpha = 0.36d0
  delta = 0.065d0
  nu = 1d0!0.88d0
  l_bar = .47d0
  ! convert variables into per period values
  delta = 1d0 - (1d0-delta)**5d0

  ! demographic parameters
  n_p   = 0.007d0
  ! convert variables into per period values
  n_p = (1d0+n_p)**5-1d0

  ! tax and transfers
  tauc   = 0.10d0
  taur   = 0.25d0
  tauy   = 0.15d0
  kappa  = 0.53d0
  sscc   = 2d0
  mu     = 1d0
  lambda = 0d0
  phi    = 0d0
  taup   = 0.10d0
  gy     = 0.19d0
  by     = 0.70d0
  ! convert variables into per period values
  by = by/5d0

  ! size of the asset grid
  a_l    = 0d0
  a_u    = 16384d0
  a_grow = 1.8d0

  ! size of the annuitiy grid
  x_l    = 0d0
  x_u    = 1024d0
  x_grow = 1.8d0

  ! size of the pension claim grid
  p_l  = 0d0
  p_u  = 2d0

  ! simulation parameters
  damp  = 0.60d0
  tol   = 1d-3
  itermax = 200

  ! compute gini
  gini_on = .true.

  ! set switches
  if (NO == 0) ent = .false.
  if (NX == 0) ann = .false.

  ! calculate initial equilibrium
  call get_SteadyState()

  stop

  ! set reform parameters
  !pen_debt = .true.
  !smopec = .true.
  phi(1:TT) = 1d0
  !lambda(1:TT) = 1d0
  !mu(1:TT) = 0d0
  !labor = .false.
  !ent = .false.

  ! calculate transition path without lsra
  lsra_on = .false.
  call get_transition()

  ! calculate transition path with lsra
  lsra_on = .true.
  call get_transition()

  ! close files
  close(21)
  close(22)
  close(23)

contains


  !##############################################################################
  ! SUBROUTINE get_SteadyState
  !
  ! Computes the initial steady state of the economy
  !##############################################################################
  subroutine get_SteadyState()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: iter, iamax(JJ), ixmax(JJ)

    ! initialize remaining variables
    call initialize()

    ! start timer
    call tick(calc)

    ! iterate until value function converges
    do iter = 1, itermax

      !call tick(calc)

      ! get factor and other prices
      call get_prices(0)

      ! solve the household problem
      call solve_household(1, 0)

      ! calculate the distribution of households over state space
      call get_distribution(0)

      ! aggregate individual decisions
      call aggregation(0)

      ! determine the government parameters
      call government(0)

      ! check the grid
      call check_grid(iamax, ixmax, 0)

      write(*,'(i4,6f8.2,2i7,f14.8)')iter, (/5d0*KK(0), CC(0), II(0)/)/YY(0)*100d0, &
        ((1d0+r(0))**0.2d0-1d0)*100d0, w(0), sum(pop_e(:, 0))/(sum(pop_w(:, 0))+sum(pop_e(:, 0)))*100d0, maxval(iamax), maxval(ixmax), DIFF(0)/YY(0)*100d0

      if(abs(DIFF(0)/YY(0))*100d0 < tol)then

        call tock(calc)
        call output(0)
        return
      endif

      !call tock(calc)

    enddo

    call tock(calc)
    call output(0)

    write(*,*)'No Convergence'

  end subroutine


  !##############################################################################
  ! SUBROUTINE get_transition
  !
  ! Computes the transition path of the economy
  !##############################################################################
  subroutine get_transition()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: iter, ij, it, itmax
    logical :: check

    ! initialize remaining variables
    if(.not. lsra_on)then
      call initialize_trn()
    else
      write(*,'(/a/)')'ITER    COMP_OLD  EFFICIENCY        DIFF'
    endif

    ! start timer
    call tick(calc)

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

      ! calculate lsra transfers if needed
      if (lsra_on) call LSRA()

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
      if(.not. lsra_on)then
        write(*,'(i4,6f8.2,f12.6)')iter, (/5d0*KK(TT), CC(TT), II(TT)/)/YY(TT)*100d0, &
          ((1d0+r(TT))**0.2d0-1d0)*100d0, w(TT), sum(pop_e(:, TT))/(sum(pop_w(:, TT))+sum(pop_e(:, TT)))*100d0, DIFF(itmax)/YY(itmax)*100d0
        check = abs(DIFF(itmax)/YY(itmax))*100d0 < tol .and. iter > 0
      else
        write(*,'(i4,2f12.6,f14.8)')iter, lsra_comp/lsra_all*100d0, &
          (Vstar**(1d0/(1d0-gamma))-1d0)*100d0,DIFF(itmax)/YY(itmax)*100d0
          check = abs(DIFF(itmax)/YY(itmax))*100d0 < tol .and. iter > 0 .and. lsra_comp/lsra_all > 0.99999d0
      endif

      ! check for convergence
      if(check)then
        call tock(calc)
        do it = 1, TT
          if (.not. lsra_on) call output(it)
        enddo
        call output_summary()
        return
      endif
    enddo

    call tock(calc)
    do it = 1, TT
      if (.not. lsra_on) call output(it)
    enddo
    call output_summary()

    write(*,*)'No Convergence'

  end subroutine


  !##############################################################################
  ! SUBROUTINE initialize
  !
  ! Initializes the remaining model parameters and variables
  !##############################################################################
  subroutine initialize()

    implicit none

    !##### OTHER VARIABLES ######################################################
    integer :: ia, ix, ip, iw, ie, is, ij
    real*8 :: adj

    write(*,'(/a/)')'INITIAL EQUILIBRIUM'
    write(*,'(a)')'ITER     K/Y     C/Y     I/Y       r       w     ent  iamax  ixmax          DIFF'

    ! initialize asset grid
    a = grid_Cons_Grow(a_l, a_u, a_grow, NA)

    ! initialize annuity grid
    if (ann) then
      x = grid_Cons_Grow(x_l, x_u, x_grow, NX)
    else
      x(0) = 0d0
    endif

    ! initialize pension claim grid
    p = grid_Cons_Equi(p_l, p_u, NP)

    ! get initial guess for savings decision
    do ij = 1, JJ
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ia = 0, NA
                  aplus(0, ia, ix, ip, iw, ie, is, ij, 0) = max(a(ia)/2d0, a(1))
                  aplus(1, ia, ix, ip, iw, ie, is, ij, 0) = max(a(ia)/2d0, a(1))
                enddo ! ia
              enddo ! ix
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is
    enddo ! ij

    ! initial guess for investment decision
    k(:, :, :, :, :, :, :, :, 0) = 1d-4

    ! initial guess for labor decision
    l(:, :, :, :, :, :, :, :, 0) = 0.33d0

    ! initial guess for annuity investment
    mx(:, :, :, :, :, :, :, :, 0) = 0.33d0

    ! distribution of skill classes
    !dist_skill(:) = (/0.2600d0, 0.5500d0, 0.1900d0/)
    dist_skill(:) = (/0.1520d0, 0.5547d0, 0.2933d0/)

    !* check for error *!
    if (sum(dist_skill(:)) /= 1d0) call error('initialize', 'distribution of skill classes in not equal to one')

		! initialize survival probabilities for middle skilled
    open(301, file='sp.dat')
    do ij = 1, JJ+1
      read(301,'(f13.8)')psi(2, ij)
    enddo
    close(301)

    ! compute survival probabilities for high/low skilled
    psi(:, 1) = psi(2, 1)
    psi(:, JJ+1) = 0d0
    adj = 23d0
    do ij = 2, JJ
      psi(1, ij) = psi(2, ij) - exp(0.33d0*(dble(ij-1)-adj))
      psi(3, ij) = psi(2, ij) + exp(0.33d0*(dble(ij-1)-adj))
    enddo

    ! set up population structure
    rpop(:, 1, 0) = dist_skill(:)
    do ij = 2, JJ
      rpop(:, ij, 0) = psi(:, ij)*rpop(:, ij-1, 0)/(1d0+n_p)
    enddo

    pop = 0d0
    do is = 1, NS
      pop(:, 0) = pop(:, 0) + rpop(is, :, 0)
    enddo

    ! set distribution of bequests
    Gama(1) = 0.0d0*pop(1, 0)
    Gama(2) = 0.0d0*pop(2, 0)
    Gama(3) = 1.0d0*pop(3, 0)
    Gama(4) = 1.2d0*pop(4, 0)
    Gama(5) = 1.4d0*pop(5, 0)
    Gama(6) = 1.6d0*pop(6, 0)
    Gama(7) = 1.8d0*pop(7, 0)
    Gama(8) = 1.8d0*pop(8, 0)
    Gama(9) = 1.6d0*pop(9, 0)
    !Gama(1:JR-1) = 1d0
    Gama(JR:JJ) = 0d0
    Gama = Gama/sum(Gama)

    ! initialize age earnings process
    open(302, file='eff.dat')
    do ij = 1, JJ
      read(302,'(3f12.8)')eff(ij, :)
    enddo
    close(302)
    eff(JE:, :) = 0d0

    call discretize_AR(0.95666d0**5d0, 0.0d0, sigma5(0.95666d0, 0.02321d0), eta(:, 1), pi_eta(:, :, 1), dist_eta(:, 1))
    eta(:, 1) = exp(eta(:, 1))/sum(dist_eta(:, 1)*exp(eta(:, 1)))

    call discretize_AR(0.95687d0**5d0, 0.0d0, sigma5(0.95687d0, 0.02812d0), eta(:, 2), pi_eta(:, :, 2), dist_eta(:, 2))
    eta(:, 2) = exp(eta(:, 2))/sum(dist_eta(:, 2)*exp(eta(:, 2)))

    call discretize_AR(0.95828d0**5d0, 0.0d0, sigma5(0.95828d0, 0.03538d0), eta(:, 3), pi_eta(:, :, 3), dist_eta(:, 3))
    eta(:, 3) = exp(eta(:, 3))/sum(dist_eta(:, 3)*exp(eta(:, 3)))

    ! initialize entrepreneurial ability
    call discretize_AR(0.93d0**5d0, -0.35d0, sigma5(0.93d0, 0.03d0), theta(:, 1), pi_theta(:, :, 1), dist_theta(:, 1))
    theta(:, 1) = exp(theta(:, 1))!/sum(dist_theta(:, 1)*exp(theta(:, 1)))

    call discretize_AR(0.94d0**5d0, -0.30d0, sigma5(0.94d0, 0.03d0), theta(:, 2), pi_theta(:, :, 2), dist_theta(:, 2))
    theta(:, 2) = exp(theta(:, 2))!/sum(dist_theta(:, 2)*exp(theta(:, 2)))

    call discretize_AR(0.95d0**5d0, -0.00d0, sigma5(0.95d0, 0.03d0), theta(:, 3), pi_theta(:, :, 3), dist_theta(:, 3))
    theta(:, 3) = exp(theta(:, 3))!/sum(dist_theta(:, 3)*exp(theta(:, 3)))

!    theta(:, 1)       = (/0.000d0, 0.290d0, 1.000d0, 1.710d0/)*1.880d0
!    theta(:, 2)       = theta(:, 1)
!    theta(:, 3)       = theta(:, 1)
!    dist_theta(:, 1)  = (/0.554d0, 0.283d0, 0.099d0, 0.064d0/)
!    dist_theta(:, 2)  = dist_theta(:, 1)
!    dist_theta(:, 3)  = dist_theta(:, 1)
!    pi_theta(1, :, 1) = (/0.780d0, 0.220d0, 0.000d0, 0.000d0/)
!    pi_theta(2, : ,1) = (/0.430d0, 0.420d0, 0.150d0, 0.000d0/)
!    pi_theta(3, :, 1) = (/0.000d0, 0.430d0, 0.420d0, 0.150d0/)
!    pi_theta(4, :, 1) = (/0.000d0, 0.000d0, 0.220d0, 0.780d0/)
!    pi_theta(:, :, 2) = pi_theta(:, :, 1)
!    pi_theta(:, :, 3) = pi_theta(:, :, 1)

    ! initial guesses for macro variables
    inc_bar(0) = 0.61d0
    BQ(0) = 0.50d0
    KC(0) = 4.93d0
    LC(0) = 5.25d0

    ! open files
    open(21, file='output.out')
    open(22, file='summary.out')
    open(23, file='output_mikrozensus2012.out')

  end subroutine


  !##############################################################################
  ! SUBROUTINE initialize_trn
  !
  ! Initializes transitional variables
  !##############################################################################
  subroutine initialize_trn()

    implicit none

    !##### OTHER VARIABLES ######################################################
    integer :: it

    write(*,'(/a/)')'TRANSITION PATH'
    write(*,'(a)')'ITER     K/Y     C/Y     I/Y       r       w     ent        DIFF'

    do it = 1, TT

      r(it) = r(0)
      w(it) = w(0)
      inc_bar(it) = inc_bar(0)
      psix(JJ, it) = psix(JJ, 0)
      pinv(it) = pinv(0)

      KK(it) = KK(0)
      KC(it) = KC(0)
      KE(it) = KE(0)
      AA(it) = AA(0)
      LC(it) = LC(0)
      HH(it) = HH(0)
      YY(it) = YY(0)
      YC(it) = YC(0)
      YE(it) = YE(0)
      CC(it) = CC(0)
      II(it) = II(0)
      GG(it) = GG(0)
      NEX(it) = NEX(0)
      BB(it) = BB(0)
      BF(it) = BF(0)
      BQ(it) = BQ(0)
      TAc(it) = TAc(0)
      TAr(it) = TAr(0)
      TAw(it) = TAw(0)
      TAy(it) = TAy(0)

      tauc(it) = tauc(0)
      taup(it) = taup(0)
      pen(:, :, it) = pen(:, :, 0)
      PP(it) = PP(0)
      PE(it) = PE(0)
      PRE(it) = PRE(0)
      PC(it) = PC(0)

      pop_w(:, it) = pop_w(:, 0)
      pop_e(:, it) = pop_e(:, 0)
      pop_r(:, it) = pop_r(:, 0)
      pop_re(:, it) = pop_re(:, 0)
      bqs(:, it) = bqs(:, 0)
      x_coh(:, it) = x_coh(:, 0)
      bx_coh(:, it) = bx_coh(:, 0)

      rpop(:, :, it) = rpop(:, :, 0)
      pop(:, it) = pop(:, 0)
      beq(:, :, it) = beq(:, :, 0)

      vv_coh(:, it) = vv_coh(:, 0)

      aplus(:, :, :, :, :, :, :, :, it) = aplus(:, :, :, :, :, :, :, :, 0)
      xplus(:, :, :, :, :, :, :, :, it) = xplus(:, :, :, :, :, :, :, :, 0)
      pplus(:, :, :, :, :, :, :, :, it) = pplus(:, :, :, :, :, :, :, :, 0)
      c(:, :, :, :, :, :, :, :, it) = c(:, :, :, :, :, :, :, :, 0)
      l(:, :, :, :, :, :, :, :, it) = l(:, :, :, :, :, :, :, :, 0)
      k(:, :, :, :, :, :, :, :, it) = k(:, :, :, :, :, :, :, :, 0)
      mx(:, :, :, :, :, :, :, :, it) = mx(:, :, :, :, :, :, :, :, 0)
      oplus(:, :, :, :, :, :, :, :, it) = oplus(:, :, :, :, :, :, :, :, 0)
      pencon(:, :, :, :, :, :, :, :, it) = pencon(:, :, :, :, :, :, :, :, 0)
      inctax(:, :, :, :, :, :, :, :, it) = inctax(:, :, :, :, :, :, :, :, 0)
      captax(:, :, :, :, :, :, :, :, it) = captax(:, :, :, :, :, :, :, :, 0)
      VV(:, :, :, :, :, :, :, :, it) = VV(:, :, :, :, :, :, :, :, 0)
      EV(:, :, :, :, :, :, :, :, it) = EV(:, :, :, :, :, :, :, :, 0)
      m(:, :, :, :, :, :, :, :, it) = m(:, :, :, :, :, :, :, :, 0)

    enddo

  end subroutine


  !##############################################################################
  ! SUBROUTINE get_prices
  !
  ! Determines factor prices, indiviudal bequests and pensions at time it
  !##############################################################################
  subroutine get_prices(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES ####################################################
    integer :: ij
    real*8 :: dueink

    ! calculate inverse goods price
    pinv(it) = 1d0/(1d0+tauc(it))

    ! calculate factor prices
    if (.not. smopec) then
      r(it) = (1d0-tauy(it))*(alpha*(KC(it)/LC(it))**(alpha-1d0)-delta)
      w(it) = (1d0-alpha)*(KC(it)/LC(it))**alpha
    endif

    ! calculate individual bequests
    beq(1, :, it) = Gama(:)*bqs(1, it)/rpop(1, :, it)
    beq(2, :, it) = Gama(:)*bqs(2, it)/rpop(2, :, it)
    beq(3, :, it) = Gama(:)*bqs(3, it)/rpop(3, :, it)

    ! calculate individual pensions
    pen(:, :, it) = 0d0
    do ij = JR, JJ
      pen(:, ij, it) = p(:)*kappa(it)*inc_bar(it)
    enddo

    ! determine the income tax system
    dueink = inc_bar(0)

    r1 = 0.286d0*dueink*2d0
    r2 = 0.456d0*dueink*2d0
    r3 = 1.786d0*dueink*2d0

    b1 = (t2-t1)/(r2-r1)
    b2 = (t3-t2)/(r3-r2)

    ! calculate interests of annuities
    do ij = 1, JJ
      psix(ij, it) = 1d0
      if (x_coh(ij, it) > 0d0) psix(ij, it) = (x_coh(ij, it) + bx_coh(ij, it))/x_coh(ij, it)
    enddo


  end subroutine


  !##############################################################################
  ! SUBROUTINE solve_household
  !
  ! Determines the solution to the household optimization problem
  !##############################################################################
  subroutine solve_household(ij_in, it_in)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it_in, ij_in

    !##### OTHER VARIABLES ####################################################
    integer :: ia, ix, ip, iw, ie, is, ij, it
    real*8 :: xy(3), fret, limit

    do ij = JJ, 1, -1

      !call tick(calc)
      !write(*,*)'Optimize for age: ', ij

      it = year(it_in, ij_in, ij)

      ! set up communication variables
      ij_com = ij
      it_com = it

      if (ij >= JE) then

        ! set up communication variables
        ie_com = 1
        iw_com = 1
        io_com = 0

        !$omp parallel do copyin(io_com, iw_com, ie_com, ij_com, it_com) collapse(3) schedule(dynamic, 1) private(xy, fret) num_threads(numthreads)
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NX
              do ia = 0, NA

                ! set up communication variables
                is_com = is
                ip_com = ip
                ix_com = ix
                ia_com = ia

                ! get initial guess for the individual choices
                xy(1) = max(aplus(0, ia, ix, ip, 1, 1, is, ij, it), 1d-4)

                call fminsearch(xy(1), fret, a_l, a_u, valuefunc_r)

                ! copy decisions
                aplus(:, ia, ix, ip,  :,  :, is, ij, it) = xy(1)
                xplus(:, ia, ix, ip,  :,  :, is, ij, it) = xplus_com
                pplus(:, ia, ix, ip,  :,  :, is, ij, it) = p(ip)
                c(:, ia, ix, ip,  :,  :, is, ij, it) = max(c_com, 1d-10)
                l(:, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                k(:, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                mx(:, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                oplus(:, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                pencon(:, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                inctax(:, ia, ix, ip,  :,  :, is, ij, it) = inctax_com
                captax(:, ia, ix, ip,  :,  :, is, ij, it) = captax_com
                VV(:, ia, ix, ip,  :,  :, is, ij, it) = -fret

              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
        !$omp end parallel do

      elseif (ij >= JR) then

        if (ent) then

          ! set up communication variables
          iw_com = 1
          io_com = 1

          !$omp parallel do copyin(io_com, iw_com, ij_com, it_com) collapse(4) schedule(dynamic, 1) private(xy, fret) num_threads(numthreads)
          do is = 1, NS
            do ie = 1, NE
              do ip = 0, NP
                do ix = 0, NX
                  do ia = 0, NA

                    ! set up communication variables
                    is_com = is
                    ie_com = ie
                    ip_com = ip
                    ix_com = ix
                    ia_com = ia

                    ! get initial guess for the individual choices
                    xy(1) = max(aplus(1, ia, ix, ip, 1, ie, is, ij, it), 1d-4)
                    xy(2) = max(k(1, ia, ix, ip, 1, ie, is, ij, it), 1d-4)
                    xy(3) = max(mx(1, ia, ix, ip, 1, ie, is, ij, it), 1d-4)

                    limit = max(1.5d0*a(ia), 1d-4)

                    call fminsearch(xy, fret, (/a_l, 0d0, x_l/), (/a_u, limit, x_u/), valuefunc_e)

                    ! copy decisions
                    aplus(1, ia, ix, ip, :, ie, is, ij, it) = xy(1)
                    xplus(1, ia, ix, ip, :, ie, is, ij, it) = xplus_com
                    k(1, ia, ix, ip, :, ie, is, ij, it) = k_com
                    pplus(1, ia, ix, ip, :, ie, is, ij, it) = pplus_com
                    c(1, ia, ix, ip, :, ie, is, ij, it) = max(c_com, 1d-10)
                    l(1, ia, ix, ip, :, ie, is, ij, it) = l_com
                    mx(1, ia, ix, ip, :, ie, is, ij, it) = mx_com
                    oplus(1, ia, ix, ip, :, ie, is, ij, it) = oplus_com
                    pencon(1, ia, ix, ip, :, ie, is, ij, it) = pencon_com
                    inctax(1, ia, ix, ip, :, ie, is, ij, it) = inctax_com
                    captax(1, ia, ix, ip, :, ie, is, ij, it) = captax_com
                    VV(1, ia, ix, ip, :, ie, is, ij, it) = -fret

                  enddo ! ia
                enddo ! ix
              enddo ! ip
            enddo ! ie
          enddo ! is
          !$omp end parallel do

        endif

        ! set up communication variables
        ie_com = 1
        iw_com = 1
        io_com = 0

        !$omp parallel do copyin(io_com, iw_com, ie_com, ij_com, it_com) collapse(3) schedule(dynamic, 1) private(xy, fret) num_threads(numthreads)
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NX
              do ia = 0, NA

                ! set up communication variables
                is_com = is
                ip_com = ip
                ix_com = ix
                ia_com = ia

                ! get initial guess for the individual choices
                xy(1) = max(aplus(0, ia, ix, ip, 1, 1, is, ij, it), 1d-4)

                call fminsearch(xy(1), fret, a_l, a_u, valuefunc_r)

                ! copy decisions
                aplus(0, ia, ix, ip,  :,  :, is, ij, it) = xy(1)
                xplus(0, ia, ix, ip,  :,  :, is, ij, it) = xplus_com
                pplus(0, ia, ix, ip,  :,  :, is, ij, it) = p(ip)
                c(0, ia, ix, ip,  :,  :, is, ij, it) = max(c_com, 1d-10)
                l(0, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                k(0, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                mx(0, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                oplus(0, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                pencon(0, ia, ix, ip,  :,  :, is, ij, it) = 0d0
                inctax(0, ia, ix, ip,  :,  :, is, ij, it) = inctax_com
                captax(0, ia, ix, ip,  :,  :, is, ij, it) = captax_com
                VV(0, ia, ix, ip,  :,  :, is, ij, it) = -fret

              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
        !$omp end parallel do

      elseif (ij >= 2) then

        !$omp parallel copyin(ij_com, it_com) private(xy, fret, limit) num_threads(numthreads)

        if (ent) then

          ! set up communication variables
          io_com = 1

          !$omp do collapse(4) schedule(dynamic, 1)
          do is = 1, NS
            do ie = 1, NE
              do iw = 1, NW
                do ip = 0, NP
                  do ix = 0, NX
                    do ia = 0, NA

                      ! set up communication variables
                      is_com = is
                      ie_com = ie
                      iw_com = iw
                      ip_com = ip
                      ix_com = ix
                      ia_com = ia

                      ! get initial guess for the individual choices
                      xy(1) = max(aplus(1, ia, ix, ip, iw, ie, is, ij, it), 1d-4)
                      xy(2) = max(k(1, ia, ix, ip, iw, ie, is, ij, it), 1d-4)
                      xy(3) = max(mx(1, ia, ix, ip, iw, ie, is, ij, it), 1d-4)

                      limit = max(1.5d0*a(ia), 1d-4)

                      call fminsearch(xy, fret, (/a_l, 0d0, x_l/), (/a_u, limit, x_u/), valuefunc_e)

                      ! copy decisions
                      aplus(1, ia, ix, ip, iw, ie, is, ij, it) = xy(1)
                      xplus(1, ia, ix, ip, iw, ie, is, ij, it) = xplus_com
                      k(1, ia, ix, ip, iw, ie, is, ij, it) = k_com
                      pplus(1, ia, ix, ip, iw, ie, is, ij, it) = pplus_com
                      c(1, ia, ix, ip, iw, ie, is, ij, it) = max(c_com, 1d-10)
                      l(1, ia, ix, ip, iw, ie, is, ij, it) = l_com
                      mx(1, ia, ix, ip, iw, ie, is, ij, it) = mx_com

                      if (mx_com > 0.001d0) write(*,*) ij, ix, ia, mx_com
                      oplus(1, ia, ix, ip, iw, ie, is, ij, it) = oplus_com
                      pencon(1, ia, ix, ip, iw, ie, is, ij, it) = pencon_com
                      inctax(1, ia, ix, ip, iw, ie, is, ij, it) = inctax_com
                      captax(1, ia, ix, ip, iw, ie, is, ij, it) = captax_com
                      VV(1, ia, ix, ip, iw, ie, is, ij, it) = -fret

                    enddo ! ia
                  enddo ! ix
                enddo ! ip
              enddo ! iw
            enddo ! ie
          enddo ! is
          !$omp end do nowait

        endif

        ! set up communication variables
        io_com = 0

        !$omp do collapse(4) schedule(dynamic, 1)
        do is = 1, NS
          do ie = 1, NE
            do iw = 1, NW
              do ip = 0, NP
                do ix = 0, NX
                  do ia = 0, NA

                    ! set up communication variables
                    is_com = is
                    ie_com = ie
                    iw_com = iw
                    ip_com = ip
                    ix_com = ix
                    ia_com = ia

                    ! get initial guess for the individual choices
                    xy(1) = max(aplus(0, ia, ix, ip, iw, ie, is, ij, it), 1d-4)
                    xy(2) = max(l(0, ia, ix, ip, iw, ie, is, ij, it), 1d-4)

                    call fminsearch(xy(:2), fret, (/a_l, 0d0/), (/a_u, 1d0/), valuefunc_w)

                    ! copy decisions
                    aplus(0, ia, ix, ip, iw, ie, is, ij, it) = xy(1)
                    xplus(0, ia, ix, ip, iw, ie, is, ij, it) = xplus_com
                    k(0, ia, ix, ip, iw, ie, is, ij, it) = k_com
                    pplus(0, ia, ix, ip, iw, ie, is, ij, it) = pplus_com
                    c(0, ia, ix, ip, iw, ie, is, ij, it) = max(c_com, 1d-10)
                    l(0, ia, ix, ip, iw, ie, is, ij, it) = l_com
                    mx(0, ia, ix, ip, iw, ie, is, ij, it) = mx_com
                    oplus(0, ia, ix, ip, iw, ie, is, ij, it) = oplus_com
                    pencon(0, ia, ix, ip, iw, ie, is, ij, it) = pencon_com
                    inctax(0, ia, ix, ip, iw, ie, is, ij, it) = inctax_com
                    captax(0, ia, ix, ip, iw, ie, is, ij, it) = captax_com
                    VV(0, ia, ix, ip, iw, ie, is, ij, it) = -fret

                  enddo ! ia
                enddo ! ix
              enddo ! ip
            enddo ! iw
          enddo ! ie
        enddo ! is
        !$omp end do
        !$omp end parallel

      elseif (ij == 1) then

        ! set up communication variables
        ip_com = 0
        ix_com = 0
        ia_com = 0
        io_com = 0

	      !$omp parallel do copyin(io_com, ia_com, ix_com, ip_com, ij_com, it_com) collapse(3) schedule(dynamic, 1) private(xy, fret) num_threads(numthreads)
        do is = 1, NS
          do ie = 1, NE
            do iw = 1, NW

              ! set up communication variables
              is_com = is
              ie_com = ie
              iw_com = iw

              ! get initial guess for the individual choices
              xy(1) = max(aplus(0, 0, 0, 0, iw, ie, is, ij, it), 1d-4)
              xy(2) = max(l(0, 0, 0, 0, iw, ie, is, ij, it), 1d-4)

              call fminsearch(xy(:2), fret, (/a_l, 0d0/), (/a_u, 1d0/), valuefunc_w)

              ! copy decisions
              aplus(:, :, :, :, iw, ie, is, ij, it) = xy(1)
              xplus(:, :, :, :, iw, ie, is, ij, it) = xplus_com
              pplus(:, :, :, :, iw, ie, is, ij, it) = pplus_com
              c(:, :, :, :, iw, ie, is, ij, it) = max(c_com, 1d-10)
              l(:, :, :, :, iw, ie, is, ij, it) = l_com
              mx(:, :, :, :, iw, ie, is, ij, it) = mx_com
              k(:, :, :, :, iw, ie, is, ij, it) = 0d0
              oplus(:, :, :, :, iw, ie, is, ij, it) = oplus_com
              pencon(:, :, :, :, iw, ie, is, ij, it) = pencon_com
              inctax(:, :, :, :, iw, ie, is, ij, it) = inctax_com
              captax(:, :, :, :, iw, ie, is, ij, it) = captax_com
              VV(:, :, :, :, iw, ie, is, ij, it) = -fret

            enddo ! iw
          enddo ! ie
        enddo ! is
	      !$omp end parallel do

      endif

      ! interpolate individual value function
      call interpolate(ij, it)

      !write(*,*)'Done!'
      !call tock(calc)

    enddo

  end subroutine


  !##############################################################################
  ! SUBROUTINE interpolate
  !
  ! Calculates the expected valuefunction of cohort ij at time it
  !##############################################################################
  subroutine interpolate(ij, it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: ij, it

    !##### OTHER VARIABLES ####################################################
    integer :: ia, ix, ip, iw, ie, is, iw_p, ie_p

    !$omp parallel do collapse(4) schedule(dynamic,1) private(iw_p, ie_p) num_threads(numthreads)
    do is = 1, NS
      do ie = 1, NE
        do iw = 1, NW
          do ip = 0, NP
            do ix = 0, NX
              do ia = 0, NA

                EV(:, ia, ix, ip, iw, ie, is, ij, it) = 0d0
                do ie_p = 1, NE
                  do iw_p = 1, NW
                    EV(0, ia, ix, ip, iw, ie, is, ij, it) = EV(0, ia, ix, ip, iw, ie, is, ij, it) &
                      +pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*VV(0, ia, ip, ix, iw_p, ie_p, is, ij, it)
                    EV(1, ia, ix, ip, iw, ie, is, ij, it) = EV(1, ia, ix, ip, iw, ie, is, ij, it) &
                      +pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*VV(1, ia, ip, ix, iw_p, ie_p, is, ij, it)
                  enddo ! iw_p
                enddo ! ie_p

                EV(:, ia, ix, ip, iw, ie, is, ij, it) &
                  = ((1d0-gamma)*EV(:, ia, ix, ip, iw, ie, is, ij, it))**(1d0/(1d0-gamma))

              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! iw
      enddo ! ie
    enddo ! is
    !$omp end parallel do

  end subroutine


  !##############################################################################
  ! SUBROUTINE get_distribution
  !
  ! Determines the invariant distribution of households at time it
  !##############################################################################
  subroutine get_distribution(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES ####################################################
    integer :: io, ia, ix, ip, iw, ie, is, ij, &
               io_p, iw_p, ie_p, ial, iar, ixr, ixl, ipl, ipr, itm
    real*8 :: varpsi, varchi, varphi

    ! get yesterdays year
    itm = year(it, 2, 1)

    ! set distribution to zero
    m(:, :, :, :, :, :, :, :, it) = 0d0

    ! get initial distribution at age 1
    do is = 1, NS
      do ie = 1, NE
        do iw = 1, NW
          m(0, 0, 0, 0, iw, ie, is, 1, it) = dist_theta(ie, is)*dist_eta(iw, is)*dist_skill(is)
        enddo ! iw
      enddo ! ie
    enddo ! is

    !write(*,*) sum(m(:, :, :, :, :, :, :, 1, it))

    ! successively compute distribution over ages
    do ij = 2, JJ

      ! iterate over yesterdays gridpoints

      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ia = 0, NA
                  do io = 0, NO

                    ! interpolate yesterday's savings decision
                    call linint_Grow(aplus(io, ia, ix, ip, iw, ie, is, ij-1, itm), &
                             a_l, a_u, a_grow, NA, ial, iar, varphi)

                    ! interpolate yesterday's annuitized assets
                    if (ann) then
                      call linint_Grow(xplus(io, ia, ix, ip, iw, ie, is, ij-1, itm), &
                             x_l, x_u, x_grow, NX, ixl, ixr, varchi)
                    else
                      ixl = 0
                      ixr = 0
                      varchi = 1d0
                    endif

                    ! interpolate today's pension claims
                    call linint_Equi(pplus(io, ia, ix, ip, iw, ie, is, ij-1, itm), &
                                       p_l, p_u, NP, ipl, ipr, varpsi)

                    ! this year's occupation
                    io_p = int(oplus(io, ia, ix, ip, iw, ie, is, ij-1, itm))

                    ! redistribute households
                    do ie_p = 1, NE
                      do iw_p = 1, NW
                        m(io_p, ial, ixl, ipl, iw_p, ie_p, is, ij, it) = &
                           m(io_p, ial, ixl, ipl, iw_p, ie_p, is, ij, it) &
                           +varphi*varchi*varpsi*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)&
                           *psi(is, ij)*m(io, ia, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(io_p, ial, ixl, ipr, iw_p, ie_p, is, ij, it) = &
                           m(io_p, ial, ixl, ipr, iw_p, ie_p, is, ij, it) &
                           +varphi*varchi*(1d0-varpsi)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                           *psi(is, ij)*m(io, ia, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(io_p, ial, ixr, ipl, iw_p, ie_p, is, ij, it) = &
                           m(io_p, ial, ixr, ipl, iw_p, ie_p, is, ij, it) &
                           +varphi*(1d0-varchi)*varpsi*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)&
                           *psi(is, ij)*m(io, ia, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(io_p, ial, ixr, ipr, iw_p, ie_p, is, ij, it) = &
                           m(io_p, ial, ixr, ipr, iw_p, ie_p, is, ij, it) &
                           +varphi*(1d0-varchi)*(1d0-varpsi)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                           *psi(is, ij)*m(io, ia, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(io_p, iar, ixl, ipl, iw_p, ie_p, is, ij, it) = &
                           m(io_p, iar, ixl, ipl, iw_p, ie_p, is, ij, it) &
                           +(1d0-varphi)*varchi*varpsi*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                           *psi(is, ij)*m(io, ia, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(io_p, iar, ixl, ipr, iw_p, ie_p, is, ij, it) = &
                           m(io_p, iar, ixl, ipr, iw_p, ie_p, is, ij, it) &
                           +(1d0-varphi)*varchi*(1d0-varpsi)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                           *psi(is, ij)*m(io, ia, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(io_p, iar, ixr, ipl, iw_p, ie_p, is, ij, it) = &
                           m(io_p, iar, ixr, ipl, iw_p, ie_p, is, ij, it) &
                           +(1d0-varphi)*(1d0-varchi)*varpsi*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                           *psi(is, ij)*m(io, ia, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                        m(io_p, iar, ixr, ipr, iw_p, ie_p, is, ij, it) = &
                           m(io_p, iar, ixr, ipr, iw_p, ie_p, is, ij, it) &
                           +(1d0-varphi)*(1d0-varchi)*(1d0-varpsi)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                           *psi(is, ij)*m(io, ia, ix, ip, iw, ie, is, ij-1, itm)/(1d0+n_p)
                      enddo ! iw_p
                    enddo ! ie_p

                  enddo ! ia
                enddo ! ix
              enddo ! ip
            enddo ! iw
          enddo ! ie
        enddo ! is
      enddo ! io

    enddo ! ij

  end subroutine


  !##############################################################################
  ! SUBROUTINE aggregation
  !
  ! Calculates aggregate quantities at time it
  !##############################################################################
  subroutine aggregation(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES ####################################################
    integer :: io, ia, ix, ip, iw, ie, is, ij, ial, iar, ixl, ixr, itm, itp
    real*8 :: LC_old, varchi, varphi

    !write(*,*)'Calculate Aggregation:'
    !call tick(calc)

    LC_old = LC(it)

    ! get tomorrow's year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    ! calculate aggregates
    AA(it)  = 0d0
    CC(it)  = 0d0
    LC(it)  = 0d0
    HH(it)  = 0d0
    KE(it)  = 0d0
    YE(it)  = 0d0
    PE(it)  = 0d0
    PRE(it) = 0d0
    PC(it)  = 0d0
    BQ(it)  = 0d0
    PP(it)  = 0d0
    AX(it)  = 0d0
    XB(it)  = 0d0
    TAc(it)   = 0d0
    TAr(it)   = 0d0
    TAw(it)   = 0d0
    TAy(it)   = 0d0
    bqs(:, it) = 0d0
    pop_w(:, it) = 0d0
    pop_r(:, it) = 0d0
    pop_re(:, it) = 0d0
    pop_e(:, it) = 0d0
    x_coh(:, it) = 0d0
    bx_coh(:, it) = 0d0
    vv_coh(:, it) = 0d0

    do ij = 1, JJ
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ix = 0, NX
                do ia = 0, NA
                  do io = 0, NO

                    call linint_Grow(aplus(io, ia, ix, ip, iw, ie, is, ij, itm), a_l, a_u, a_grow, NA, ial, iar, varphi)
                    if (ann) then
                      call linint_Grow(xplus(io, ia, ix, ip, iw, ie, is, ij, itm), x_l, x_u, x_grow, NX, ixl, ixr, varchi)
                    else
                      ixl = 0
                      ixr = 0
                      varchi = 1d0
                    endif

                    AA(it) = AA(it) + (varphi*a(ial) + (1d0-varphi)*a(iar) + varchi*x(ixl) + (1d0-varchi)*x(ixr)) &
                              *m(io, ia, ix, ip, iw, ie, is, ij, itm)/(1d0+n_p)
                    x_coh(ij, it) = x_coh(ij, it) + x(ix) &
                              *m(io, ia, ix, ip, iw, ie, is, ij, it)
                    bx_coh(ij, it) = bx_coh(ij, it) + x(ix) &
                             *m(io, ia, ix, ip, iw, ie, is, ij, it)/psi(is, ij)*(1d0-psi(is, ij))
                    CC(it) = CC(it) + c(io, ia, ix, ip, iw, ie, is, ij, it) &
                              *m(io, ia, ix, ip, iw, ie, is, ij, it)
                    if (io == 0 .and. ij < JR) then
                      LC(it) = LC(it) + eff(ij, is)*eta(iw, is)*l(io, ia, ix, ip, iw, ie, is, ij, it) &
                                *m(io, ia, ix, ip, iw, ie, is, ij, it)
                      HH(it) = HH(it) + l(io, ia, ix, ip, iw, ie, is, ij, it) &
                                *m(io, ia, ix, ip, iw, ie, is, ij, it)
                      PC(it) = PC(it) + min(w(it)*eff(ij, is)*eta(iw, is)*l(io, ia, ix, ip, iw, ie, is, ij, it), sscc(it)*inc_bar(it)) &
                                *m(io, ia, ix, ip, iw, ie, is, ij, it)
                    else
                      KE(it) = KE(it) + k(io, ia, ix, ip, iw, ie, is, ij, it)*m(io, ia, ix, ip, iw, ie, is, ij, it)
                      YE(it) = YE(it) + theta(ie, is)*(k(io, ia, ix, ip, iw, ie, is, ij, it)**alpha*(eff(ij, is)*l_bar)**(1d0-alpha))**nu &
                                *m(io, ia, ix, ip, iw, ie, is, ij, it)
                      if (ij < JR) PC(it) = PC(it) + phi(it)*min(profent(k(io, ia, ix, ip, iw, ie, is, ij, it), ij, ia, is, ie, it), sscc(it)*inc_bar(it)) &
                                *m(io, ia, ix, ip, iw, ie, is, ij, it)
                      PE(it) = PE(it) + profent(k(io, ia, ix, ip, iw, ie, is, ij, it), ij, ia, is, ie, it) &
                                *m(io, ia, ix, ip, iw, ie, is, ij, it)
                      if (ij >= JR) then
                        PRE(it) = PRE(it) + profent(k(io, ia, ix, ip, iw, ie, is, ij, it), ij, ia, is, ie, it) &
                                  *m(io, ia, ix, ip, iw, ie, is, ij, it)
                      endif
                    endif
                    PP(it) = PP(it) + pen(ip, ij, it)*m(io, ia, ix, ip, iw, ie, is, ij, it)
                    bqs(is, it) = bqs(is, it) + (1d0+r(it))*(varphi*a(ial) + (1d0-varphi)*a(iar))*(1d0-psi(is, ij+1)) &
                                  *m(io, ia, ix, ip, iw, ie, is, ij, itm)/(1d0+n_p)
                    TAc(it) = TAc(it) + tauc(it)*c(io, ia, ix, ip, iw, ie, is, ij, it) &
                              *m(io, ia, ix, ip, iw, ie, is, ij, it)
                    TAr(it) = TAr(it) + captax(io, ia, ix, ip, iw, ie, is, ij, it)*m(io, ia, ix, ip, iw, ie, is, ij, it)
                    TAw(it) = TAw(it) + inctax(io, ia, ix, ip, iw, ie, is, ij, it)*m(io, ia, ix, ip, iw, ie, is, ij, it)
                    if (io == 1 .and. ij < JR) then
                      pop_e(is, it) = pop_e(is, it) + m(io, ia, ix, ip, iw, ie, is, ij, it)
                    elseif (io == 1 .and. ij >= JR) then
                      pop_re(is, it) = pop_re(is, it) + m(io, ia, ix, ip, iw, ie, is, ij, it)
                    elseif (io == 0 .and. ij < JR) then
                      pop_w(is, it) = pop_w(is, it) + m(io, ia, ix, ip, iw, ie, is, ij, it)
                    else
                      pop_r(is, it) = pop_r(is, it) + m(io, ia, ix, ip, iw, ie, is, ij, it)
                    endif
                    vv_coh(ij, it) = vv_coh(ij, it) + VV(io, ia, ix, ip, iw, ie, is, ij, it) &
                                    *m(io, ia, ix, ip, iw, ie, is, ij, it)

                  enddo ! io
                enddo ! ia
              enddo ! ix
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is

      vv_coh(ij, it) = vv_coh(ij, it)/sum(m(:, :, :, :, :, :, :, ij, it))

    enddo ! ij


    ! damping and other quantities
    LC(it) = damp*LC(it) + (1d0-damp)*LC_old

    if (smopec) then
      KC(it) = damp*(LC(it)*((r(it)/(1d0-tauy(it))+delta)/alpha)**(1d0/(alpha-1d0)))+(1d0-damp)*KC(it)
      BF(it) = AA(it) - KE(it) - KC(it) - BB(it) - BP(it) - BA(it)
      NEX(it) = (n_p-r(it))*BF(it)
    else
      KC(it) = damp*(AA(it) - KE(it) - BB(it) - BP(it) - BA(it))+(1d0-damp)*KC(it)
      BF(it) = 0d0
      NEX(it) = 0d0
    endif

    BQ(it) = sum(bqs(:, it))

    KK(it) = KC(it) + KE(it)

    II(it) = (1d0+n_p)*KK(itp) - (1d0-delta)*KK(it)
    YC(it) = KC(it)**alpha * LC(it)**(1d0-alpha)
    YY(it) = YC(it) + YE(it)

    TAy(it) = TAy(it) + tauy(it)*(YC(it) - delta*KC(it) - w(it)*LC(it))

    inc_bar(it) = (w(it)*LC(it) + PE(it))/(sum(pop_w(:, it)) + sum(pop_e(:, it)))

    ! get difference on goods market
    DIFF = YY(it)-CC(it)-II(it)-GG(it)-NEX(it)

    !write(*,*)'Done!'
    !call tock(calc)

  end subroutine


  !##############################################################################
  ! SUBROUTINE government
  !
  ! Calculates government parameters at time it
  !##############################################################################
  subroutine government(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES ####################################################
    integer :: itp, itt
    real*8 :: expend, PP_bar, PC_bar

    ! last year
    itp = year(it, 1, 2)

    ! set government quantities and pension payments
    GG(it) = gy*YY(0)
    BB(it) = by*YY(0)

    ! calculate government expenditure
    expend = GG(it) + (1d0+r(it))*BB(it) - (1d0+n_p)*BB(itp)

    ! get budget balancing tax rate
    tauc(it) = (expend - (TAr(it) + TAw(it) + TAy(it)))/CC(it)

    ! get budget balancing contribution rate to pension system
    if (.not. pen_debt) taup(it) = PP(it)/PC(it)

    if (pen_debt .and. it == 1) then

      PP_bar = PP(TT)/(r(TT)-n_p)*(1d0+r(TT))
      PC_bar = PC(TT)/(r(TT)-n_p)*(1d0+r(TT))

      do itt = TT-1, 1, -1
        PP_bar = (1d0+n_p)*PP_bar/(1d0+r(itt+1)) + PP(itt)
        PC_bar = (1d0+n_p)*PC_bar/(1d0+r(itt+1)) + PC(itt)
      enddo

      if (PP_bar > 0d0) then
        taup(it) = PP_bar/PC_bar
      else
        taup(it) = taup(0)
      endif

      BP(1:TT) = 0d0

    elseif (pen_debt .and. it > 1) then

      taup(it) = taup(1)
      if (it > 1) BP(it) = ((1d0+r(it-1))*BP(it-1) + PP(it-1) - taup(it-1)*PC(it-1))/(1d0+n_p)
      if (it == TT) BP(it) = (PP(TT) - taup(TT)*PC(TT))/(n_p - r(TT))

    endif

  end subroutine


  !##############################################################################
  ! SUBROUTINE LSRA
  !
  ! Calculates LSRA payments
  !##############################################################################
  subroutine LSRA()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: io, ia, ix, ip, iw, ie, is, ij, it
    real*8 :: VV_today, VV_target, dVV_da, v_tilde
    real*8 :: pv_today, pv_target, pv_trans

    ! initialize variables
    SV(:) = 0d0

    ! initialize counters
    lsra_comp   = 0d0
    lsra_all    = 0d0

    do ij = 2, JJ
      do ia = 0, NA
        do ix = 0, NX
          do ip = 0, NP
            do iw = 1, NW
              do ie = 1, NE
                do is = 1, NS
                  do io = 0, NO

                    ! do not do anything for an agent at retirement without pension and savings
                    if(ij >= JR .and. ia == 0 .and. (pen(ip, ij, 0) <= 1d-10 .or. pen(ip, ij, 1) <= 1d-10))then
                      v(io, ia, ix, ip, iw, ie, is, ij, 1) = 0d0
                      cycle
                    endif

                    ! get today's utility
                    VV_today = VV(io, ia, ix, ip, iw, ie, is, ij, 1)
!                    if (VV_today >= 0d0)write(*,*)'VV_today', VV_today
!                    if (VV_today <= -1d10)write(*,*)'VV_today', VV_today

                    ! get target utility
                    VV_target = VV(io, ia, ix, ip, iw, ie, is, ij, 0)
!                    if (VV_target >= 0d0)write(*,*)'VV_target', VV_target
!                    if (VV_target <= -1d10)write(*,*)'VV_target', VV_target

                    ! get derivative of the value function
                    dVV_da = margu(c(io, ia, ix, ip, iw, ie, is, ij, 1), l(io, ia, ix, ip, iw, ie, is, ij, 1), 1)
!                    if (dVV_da < 0d0)write(*,*)'dVV_da', dVV_da

                    ! calculate change in transfers
                    v_tilde = (VV_target-VV_today)/dVV_da

                    ! check whether individual is already compensated
                    lsra_all = lsra_all + m(io, ia, ix, ip, iw, ie, is, ij, 1)*pop(ij, 1)
                    if(abs((VV_today-VV_target)/VV_target) < tol) &
                      lsra_comp = lsra_comp + m(io, ia, ix, ip, iw, ie, is, ij, 1)*pop(ij, 1)

                    ! calculate total transfer
                    v(io, ia, ix, ip, iw, ie, is, ij, 1) = v(io, ia, ix, ip, iw, ie, is, ij, 1) + damp*damp*v_tilde

                    ! aggregate transfers by cohort
                    SV(1) = SV(1) + v(io, ia, ix, ip, iw, ie, is, ij, 1)*m(io, ia, ix, ip, iw, ie, is, ij, 1)*pop(ij, 1)

                  enddo ! io
                enddo ! is
              enddo ! ie
            enddo ! iw
          enddo ! ip
        enddo ! ix
      enddo ! ia
    enddo ! ij

    ! initialize present value variables
    pv_today = 0d0
    pv_target = 0d0
    pv_trans = 0d0

    ! calculate present value of utility changes (in monetary values)
    do it = TT, 1, -1

      ! get today's ex ante utility
      VV_today = damp*damp*vv_coh(1, it)

      ! get damped target utility
      VV_target = damp*damp*vv_coh(1, 0)

      ! get derivative of expected utility function
      dVV_da = 0d0
      do iw = 1, NW
        do ie = 1, NE
          do is = 1, NS
            dVV_da = dVV_da + margu(c(0, 0, 0, 0, iw, ie, is, 1, it), l(0, 0, 0, 0, iw, ie, is, 1, it), it)*m(0, 0, 0, 0, iw, ie, is, 1, it)
          enddo ! is
        enddo ! ie
      enddo ! iw

      ! calculate present values
      if(it == TT)then
        pv_today  = VV_today/dVV_da  *(1d0+r(it))/(r(it)-n_p)
        pv_target = VV_target/dVV_da   *(1d0+r(it))/(r(it)-n_p)
        pv_trans  = v(0, 0, 0, 0, 1, 1, 1, 1, it)*(1d0+r(it))/(r(it)-n_p)
      else
        pv_today  = pv_today *(1d0+n_p)/(1d0+r(it+1)) + VV_today/dVV_da
        pv_target = pv_target*(1d0+n_p)/(1d0+r(it+1)) + VV_target/dVV_da
        pv_trans  = pv_trans *(1d0+n_p)/(1d0+r(it+1)) + v(0, 0, 0, 0, 1, 1, 1, 1, it)
      endif
    enddo

    ! calculate the constant utility gain/loss for future generations
    Vstar = (pv_today-pv_trans-SV(1))/pv_target

    ! calculate compensation payments for future cohorts
    do it = TT, 1, -1

      ! get today's ex ante utility
      VV_today = damp*damp*vv_coh(1, it)

      ! get target utility
      VV_target = damp*damp*vv_coh(1, 0)*Vstar

      ! get derivative of expected utility function
      dVV_da = 0d0
      do iw = 1, NW
        do ie = 1, NE
          do is = 1, NS
            dVV_da = dVV_da + margu(c(0, 0, 0, 0, iw, ie, is, 1, it), l(0, 0, 0, 0, iw, ie, is, 1, it), it)*m(0, 0, 0, 0, iw, ie, is, 1, it)
          enddo ! is
        enddo ! ie
      enddo ! iw

      ! compute change in transfers (restricted)
      v_tilde = (VV_target-VV_today)/dVV_da

      ! calculate cohort transfer level
      v(0, 0, 0, 0, :, :, :, 1, it) = v(0, 0, 0, 0, :, :, :, 1, it) + v_tilde

      ! aggregate transfers
      SV(it) = SV(it) + v(0, 0, 0, 0, 1, 1, 1, 1, it)*pop(1, it)

    enddo

    ! determine sequence of LSRA debt/savings
    BA(2) = SV(1)/(1d0+n_p)
    do it = 3, TT
      BA(it) = ((1d0+r(it-1))*BA(it-1) + SV(it-1))/(1d0+n_p)
    enddo

  end subroutine


  !##############################################################################
  ! SUBROUTINE output
  !
  ! Writes output of year it to file 21
  !##############################################################################
  subroutine output(it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: it

    !##### OTHER VARIABLES ####################################################
    integer :: io, ia, ix, ip, iw, ie, is, ij, io_p, ixmax(JJ), iamax(JJ)
    real*8 :: c_coh(0:1, JJ, 0:TT), a_coh(0:1, JJ, 0:TT), ax_coh(0:1, JJ, 0:TT), k_coh(JJ, 0:TT)
    real*8 :: inc_coh(0:1, JJ, 0:TT), o_coh(0:1, 0:1, JJ, 0:TT), flc_coh(JJ, 0:TT)
    real*8, allocatable :: wealth(:, :, :, :, :, :, :, :), grossinc(:, :, :, :, :, :, :, :), netinc(:, :, :, :, :, :, :, :)
    real*8 :: life_exp(NS), punb(NS, JJ)

    if(allocated(wealth))deallocate(wealth)
    if(allocated(grossinc))deallocate(grossinc)
    if(allocated(netinc))deallocate(netinc)
    allocate(wealth(0:NO, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ))
    allocate(grossinc(0:NO, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ))
    allocate(netinc(0:NO, 0:NA, 0:NX, 0:NP, NW, NE, NS, JJ))

    life_exp = 0d0
    do is = 1, NS
      punb(is, 1) = psi(is, 1)
      life_exp(is) = life_exp(is) + 22d0*punb(is, 1)*(1d0-psi(is, 2))
      do ij = 2, JJ
        punb(is, ij) = punb(is, ij-1)*psi(is, ij)
        life_exp(is) = life_exp(is) + (22d0 + 5d0*dble(ij-1))*punb(is, ij)*(1d0-psi(is, ij+1))
      enddo ! ij
    enddo ! is

    c_coh(:, :, it) = 0d0
    a_coh(:, :, it) = 0d0
    ax_coh(:, :, it) = 0d0
    k_coh(:, it) = 0d0
    inc_coh(:, :, it) = 0d0
    o_coh(:, :, :, it) = 0d0
    os_coh(:, :, :, :, it) = 0d0
    flc_coh(:, it) = 0d0

    do ij = 1, JJ
      do ia = 0, NA
        do ix = 0, NX
          do ip = 0, NP
            do iw = 1, NW
              do ie = 1, NE
                do is = 1, NS
                  do io = 0, NO

                    ! Cohort average variables
                    c_coh(io, ij, it) = c_coh(io, ij, it) + c(io, ia, ix, ip, iw, ie, is, ij, it) &
                                        *m(io, ia, ix, ip, iw, ie, is, ij, it)
                    a_coh(io, ij, it) = a_coh(io, ij, it) + a(ia)*m(io, ia, ix, ip, iw, ie, is, ij, it)
                    ax_coh(io, ij, it) = ax_coh(io, ij, it) + x(ix)*m(io, ia, ix, ip, iw, ie, is, ij, it)
                    io_p = int(oplus(io, ia, ix, ip, iw, ie, is, ij, it))
                    o_coh(io, io_p, ij, it) = o_coh(io, io_p, ij, it) + m(io, ia, ix, ip, iw, ie, is, ij, it)
                    os_coh(io, io_p, is, ij, it) = os_coh(io, io_p, is, ij, it) &
                                     + m(io, ia, ix, ip, iw, ie, is, ij, it)
                    if (io == 1) then
                      k_coh(ij, it) = k_coh(ij, it) + k(io, ia, ix, ip, iw, ie, is, ij, it) &
                                      *m(io, ia, ix, ip, iw, ie, is, ij, it)
                      inc_coh(io, ij, it) = inc_coh(io, ij, it) + profent(k(io, ia, ix, ip, iw, ie, is, ij, it), ij, ia, is, ie, it) &
                                            *m(io, ia, ix, ip, iw, ie, is, ij, it)
                    else
                      inc_coh(io, ij, it) = inc_coh(io, ij, it) + w(it)*eff(ij, is)*eta(iw, is)*l(io, ia, ix, ip, iw, ie, is, ij, it) &
                                            *m(io, ia, ix, ip, iw, ie, is, ij, it)
                    endif
                    if (aplus(io, ia, ix, ip, iw, ie, is, ij, it) <= 1d-10) then
                      flc_coh(ij, it) = flc_coh(ij, it) + m(io, ia, ix, ip, iw, ie, is, ij, it)
                    endif

                    ! Individual variables
                    if (io == 0) then
                      grossinc(io, ia, ix, ip, iw, ie, is, ij) = a(ia)*r(it) + pen(ip, ij, it) + eff(ij, is)*eta(iw, is)*l(io, ia, ip, ix, iw, ie, is, ij, it)*w(it)
                    else
                      grossinc(io, ia, ix, ip, iw, ie, is, ij) = max(a(ia)-k(1, ia, ip, ix, iw, ie, is, ij, it), 0d0)*r(it) + pen(ip, ij, it) + profent(k(io, ia, ip, ix, iw, ie, is, ij, it), ij, ia, is, ie, it)
                    endif
                    netinc(io, ia, ix, ip, iw, ie, is, ij) = grossinc(io, ia, ix, ip, iw, ie, is, ij) - captax(io, ia, ix, ip, iw, ie, is, ij, it) - inctax(io, ia, ix, ip, iw, ie, is, ij, it) - pencon(io, ia, ix, ip, iw, ie, is, ij, it)
                    wealth(io, ia, ix, ip, iw, ie, is, ij) = a(ia)

                  enddo ! io
                enddo ! is
              enddo ! ie
            enddo ! iw
          enddo ! ip
        enddo ! ix
      enddo ! ia
    enddo ! ij

    do ij = 1, JJ
      c_coh(0, ij, it) = c_coh(0, ij, it)/sum(m(0, :, :, :, :, :, :, ij, it))
      c_coh(1, ij, it) = c_coh(1, ij, it)/sum(m(1, :, :, :, :, :, :, ij, it))
      a_coh(0, ij, it) = a_coh(0, ij, it)/sum(m(0, :, :, :, :, :, :, ij, it))
      a_coh(1, ij, it) = a_coh(1, ij, it)/sum(m(1, :, :, :, :, :, :, ij, it))
      ax_coh(0, ij, it) = ax_coh(0, ij, it)/sum(m(0, :, :, :, :, :, :, ij, it))
      ax_coh(1, ij, it) = ax_coh(1, ij, it)/sum(m(1, :, :, :, :, :, :, ij, it))
      inc_coh(0, ij, it) = inc_coh(0, ij, it)/sum(m(0, :, :, :, :, :, :, ij, it))
      inc_coh(1, ij, it) = inc_coh(1, ij, it)/sum(m(1, :, :, :, :, :, :, ij, it))
      k_coh(ij, it) = k_coh(ij, it)/sum(m(1, :, :, :, :, :, :, ij, it))
      flc_coh(ij, it) = flc_coh(ij, it)/sum(m(:, :, :, :, :, :, :, ij, it))
      do is = 1, NS
        os_coh(0, :, is, ij, it) = os_coh(0, :, is, ij, it)/sum(m(:, :, :, :, :, :, is, ij, it))
        os_coh(1, :, is, ij, it) = os_coh(1, :, is, ij, it)/sum(m(:, :, :, :, :, :, is, ij, it))
      end do ! is
    enddo ! ij

    ! Output
    write(21,'(a, i3/)')'EQUILIBRIUM YEAR ', it
    write(21,'(a)')'CAPITAL       KK      KC      KE      AA       r    p.a.'
    write(21,'(8x,6f8.2)')KK(it), KC(it), KE(it), AA(it), r(it), ((1d0+r(it))**(1d0/5d0)-1d0)*100d0
    write(21,'(a,4f8.2/)')'(in %)  ',(/KK(it), KC(it), KE(it), AA(it)/)/YY(it)*500d0

    write(21,'(a)')'LABOR         LC       w     inc   l_bar'
    write(21,'(8x,4f8.2/)')LC(it), w(it), inc_bar(it), HH(it)/sum(pop_w(:, it))

    write(21,'(a)')'GOODS         YY      YC      YE      CC      II      NX    DIFF'
    write(21,'(8x,6f8.2,f8.3)')YY(it), YC(it), YE(it), CC(it), II(it), NEX(it), diff(it)
    write(21,'(a,6f8.2,f8.3/)')'(in %)  ',(/YY(it), YC(it), YE(it), CC(it), II(it), NEX(it), diff(it)/)/YY(it)*100d0

    write(21,'(a)')'GOV         TAUC    TAUR    TAUW    TAUY   TOTAL      GG      BB      BF'
    write(21,'(8x,8f8.2)')TAc(it), TAr(it), TAw(it), TAy(it), TAc(it)+TAr(it)+TAw(it)+TAy(it), GG(it), BB(it), BF(it)
    write(21,'(a,8f8.2)')'(in %)  ',(/TAc(it), TAr(it), TAw(it), TAy(it)/)/(TAc(it)+TAr(it)+TAw(it)+TAy(it))*100d0, &
               (/TAc(it)+TAr(it)+TAw(it)+TAy(it), GG(it), BB(it)*5d0, BF(it)*5d0/)/YY(it)*100d0
    write(21,'(a,4f8.2/)')'(rate)  ',(/tauc(it), taur(it), TAw(it)/(w(it)*LC(it)+ PE(it) + PRE(it)), tauy(it)/)*100d0

    write(21,'(a)')'PENS        TAUP     PEN      PP      PC      BQ'
    write(21,'(8x,5f8.2)')taup(it)*PC(it), PP(it)/sum(pop(JR:, it)), PP(it), PC(it), BQ(it)
    write(21,'(a,5f8.2/)')'(in %)  ',(/taup(it), kappa(it), PP(it)/YY(it), PC(it)/YY(it), BQ(it)/YY(it)/)*100d0

    write(21,'(a)')'INCOME     TOTAL     WOR     ENT ENT/WOR     y_w     y_e y_e/y_w'
    write(21, '(8x,7f8.2)')w(it)*LC(it) + PE(it) + PRE(it), w(it)*LC(it), PE(it) + PRE(it), (PE(it)+PRE(it))/(w(it)*LC(it)), w(it)*LC(it)/sum(pop_w(:, it)), (PE(it) + PRE(it))/(sum(pop_e(:, it) + pop_re(:, it))), (PE(it) + PRE(it))/sum(pop_e(:, it) + pop_re(:, it))/w(it)*LC(it)/sum(pop_w(:, it))
    write(21, '(a,4f8.2/)')'(in %)  ',(/w(it)*LC(it) + PE(it) + PRE(it), w(it)*LC(it), PE(it) + PRE(it)/)/(w(it)*LC(it) + PE(it) + PRE(it))*100d0, (PE(it)+PRE(it))/(w(it)*LC(it))*100d0

    write(21,'(a)')'POP        TOTAL     65-     65+ 65+/65-'
    write(21,'(8x,3f8.2)')sum(pop(:, it)), sum(pop(:JR-1, it)), sum(pop(JR:, it))
    write(21,'(a,4f8.2/)')'(in %)  ',(/sum(pop(:, it)), sum(pop(:JR-1, it)), sum(pop(JR:, it))/)/sum(pop(:, it))*100d0, &
                sum(pop(JR:, it))/sum(pop(:JR-1, it))*100d0

    write(21,'(a)')'           WORFO     WOR     ENT    WOR1    ENT1    WOR2    ENT2    WOR3    ENT3'
    write(21,'(8x,9f8.2)')sum(pop_w(:, it))+sum(pop_e(:, it)), sum(pop_w(:, it)), sum(pop_e(:, it)), &
                pop_w(1, it), pop_e(1, it), pop_w(2, it), pop_e(2, it), pop_w(3, it), pop_e(3, it)
    write(21, '(a, 9f8.2/)')'(in %)  ',(/sum(pop_w(:, it)+pop_e(:, it)), sum(pop_w(:, it)), sum(pop_e(:, it))/)/(sum(pop_w(:, it)+pop_e(:, it)))*100d0, &
                (/pop_w(1, it), pop_e(1, it)/)/(pop_w(1, it)+pop_e(1, it))*100d0, &
                (/pop_w(2, it), pop_e(2, it)/)/(pop_w(2, it)+pop_e(2, it))*100d0, &
                (/pop_w(3, it), pop_e(3, it)/)/(pop_w(3, it)+pop_e(3, it))*100d0

    write(21,'(a)')'LIFE       j_bar  j_bar1  j_bar2  j_bar3'
    write(21,'(8x,4f8.2/)')sum(life_exp*dist_skill), life_exp(:)

    ! check for the maximium grid point used

    call check_grid(iamax, ixmax, it)

    write(21, '(a,a)')' IJ   CONSw   CONSe   ASSw     ASSe   ASSxw   ASSxe    INCw    INCe    INVe     ENT    ENTs1  ENTs2   ENTs3', &
        '    ENTn    WORn     FLC      VALUE  IAMAX  IXMAX'
    write(21,'(a)')'------------------------------------------------------------------------------------------------------------------------------------------------------------'
    do ij = 1, JJ
      write(21, '(i3, 16f8.3, f11.3, 2i7)')ij, c_coh(0, ij, it), c_coh(1, ij, it), a_coh(0, ij, it), a_coh(1, ij, it), ax_coh(0, ij, it), ax_coh(1, ij, it), inc_coh(0, ij, it), inc_coh(1, ij, it), &
          k_coh(ij, it), sum(o_coh(1, :, ij, it)), sum(os_coh(1, :, 1, ij, it)), sum(os_coh(1, :, 2, ij, it)), &
          sum(os_coh(1, :, 3, ij, it)), o_coh(0, 1, ij, it), o_coh(1, 0, ij, it), flc_coh(ij, it), vv_coh(ij, it), iamax(ij), ixmax(ij)
      if (ij == JR-1) write(21,'(a)')'------------------------------------------------------------------------------------------------------------------------------------------------------------'
    enddo

    if (gini_on .and. (it == 0 .or. TT)) then

      write(21,'(a)')'------------------------------------------------------------------------------------------------------------------------------------------------------------'
      write(21,'(a/)')' '

      write(21, '(a)')'WEALTH    GINI     1%     5%    10%    25%    50%    75%    90%    95%    99%'
      write(21,'(4x, f10.3, 9f7.2/)')gini(reshape(wealth(:, :, :, :, :, :, :, :), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), reshape(m(:, :, :, :, :, :, :, :, it), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/))), &
                                    percentiles(reshape(wealth(:, :, :, :, :, :, :, :), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), reshape(m(:, :, :, :, :, :, :, :, it), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))

      write(21, '(a)')'GROSSINC  GINI     1%     5%    10%    25%    50%    75%    90%    95%    99%'
      write(21,'(4x, f10.3, 9f7.2/)')gini(reshape(grossinc(:, :, :, :, :, :, :, :), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), reshape(m(:, :, :, :, :, :, :, :, it), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/))), &
                                     percentiles(reshape(grossinc(:, :, :, :, :, :, :, :), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), reshape(m(:, :, :, :, :, :, :, :, it), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))

      write(21, '(a)')'NETINC    GINI     1%     5%    10%    25%    50%    75%    90%    95%    99%'
      write(21,'(4x, f10.3, 9f7.2/)')gini(reshape(netinc(:, :, :, :, :, :, :, :), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), reshape(m(:, :, :, :, :, :, :, :, it), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/))), &
                                     percentiles(reshape(netinc(:, :, :, :, :, :, :, :), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), reshape(m(:, :, :, :, :, :, :, :, it), (/2*NS*NE*NW*(NP+1)*(NX+1)*(NA+1)*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))

    endif

    write(21,'(a)')'------------------------------------------------------------------------------------------------------------------------------------------------------------'
    write(21,'(a/)')'------------------------------------------------------------------------------------------------------------------------------------------------------------'



    !##### Output Mikrozensus 2012 ############################################

    write(23,'(/a)')'Shares of entrepreneurs'
    do ij = 1, JJ
      write(23,'(i3, 4f8.2)') ij, (/sum(os_coh(1, :, 1, ij, it)), &
                                    sum(os_coh(1, :, 2, ij, it)), &
                                    sum(os_coh(1, :, 3, ij, it)), &
                                    sum(os_coh(1, :, :, ij, it)*reshape((/rpop(:, ij, it), rpop(:, ij, it)/), (/2, NS/)))/pop(ij, it)/)*100d0
    enddo

    write(23,'(/a)')'Entry rates'
    do ij = 1, JJ
      write(23,'(i3, 4f8.2)') ij, (/os_coh(0, 1, 1, ij, it), &
                                    os_coh(0, 1, 2, ij, it), &
                                    os_coh(0, 1, 3, ij, it), &
                                    sum(os_coh(0, 1, :, ij, it)*rpop(:, ij, it))/pop(ij, it)/)*100d0
    enddo

    write(23,'(/a)')'Exit rates'
    do ij = 1, JJ
      write(23,'(i3, 4f8.2)') ij, (/os_coh(1, 0, 1, ij, it), &
                                   os_coh(1, 0, 2, ij, it), &
                                   os_coh(1, 0, 3, ij, it), &
                                   sum(os_coh(1, 0, :, ij, it)*rpop(:, ij, it))/pop(ij, it)/)*100d0
    enddo

    write(23,'(/a)')'Mean netinc (worker)'
    do ij = 1, JJ
      write(23,'(i3, 4f12.6)') ij, sum(reshape(netinc(0, :, :, :, :, :, 1, ij), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/))*reshape(m(0, :, :, :, :, :, 1, ij, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/)))/sum(m(0, :, :, :, :, :, 1, ij, it)), &
                                   sum(reshape(netinc(0, :, :, :, :, :, 2, ij), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/))*reshape(m(0, :, :, :, :, :, 2, ij, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/)))/sum(m(0, :, :, :, :, :, 2, ij, it)), &
                                   sum(reshape(netinc(0, :, :, :, :, :, 3, ij), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/))*reshape(m(0, :, :, :, :, :, 3, ij, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/)))/sum(m(0, :, :, :, :, :, 3, ij, it)), &
                                   sum(reshape(netinc(0, :, :, :, :, :, :, ij), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS/))*reshape(m(0, :, :, :, :, :, :, ij, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS/)))/sum(m(0, :, :, :, :, :, :, ij, it))
    enddo

    write(23,'(/a)')'Mean netinc (entrepreneur)'
    do ij = 1, JJ
      write(23,'(i3, 4f12.6)') ij, sum(reshape(netinc(1, :, :, :, :, :, 1, ij), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/))*reshape(m(1, :, :, :, :, :, 1, ij, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/)))/sum(m(1, :, :, :, :, :, 1, ij, it)), &
                                   sum(reshape(netinc(1, :, :, :, :, :, 2, ij), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/))*reshape(m(1, :, :, :, :, :, 2, ij, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/)))/sum(m(1, :, :, :, :, :, 2, ij, it)), &
                                   sum(reshape(netinc(1, :, :, :, :, :, 3, ij), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/))*reshape(m(1, :, :, :, :, :, 3, ij, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE/)))/sum(m(1, :, :, :, :, :, 3, ij, it)), &
                                   sum(reshape(netinc(1, :, :, :, :, :, :, ij), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS/))*reshape(m(1, :, :, :, :, :, :, ij, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS/)))/sum(m(1, :, :, :, :, :, :, ij, it))
    enddo

    write(23,'(/a)')'GINI (netinc)'
    write(23,'(4f8.4)') gini(reshape(netinc(0, :, :, :, :, :, 1, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(0, :, :, :, :, :, 1, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(0, :, :, :, :, :, 2, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(0, :, :, :, :, :, 2, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(0, :, :, :, :, :, 3, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(0, :, :, :, :, :, 3, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(0, :, :, :, :, :, :, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), reshape(m(0, :, :, :, :, :, :, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)))
    write(23,'(4f8.4)') gini(reshape(netinc(1, :, :, :, :, :, 1, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(1, :, :, :, :, :, 1, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(1, :, :, :, :, :, 2, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(1, :, :, :, :, :, 2, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(1, :, :, :, :, :, 3, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(1, :, :, :, :, :, 3, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(1, :, :, :, :, :, :, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), reshape(m(1, :, :, :, :, :, :, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)))
    write(23,'(4f8.4)') gini(reshape(netinc(:, :, :, :, :, :, 1, :), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(:, :, :, :, :, :, 1, :, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(:, :, :, :, :, :, 2, :), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(:, :, :, :, :, :, 2, :, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(:, :, :, :, :, :, 3, :), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(:, :, :, :, :, :, 3, :, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/))), &
                        gini(reshape(netinc(:, :, :, :, :, :, :, :), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), reshape(m(:, :, :, :, :, :, :, :, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)))

    write(23,'(/a)')'Percentiles (netinc)'
    write(23,'(9f8.4)') percentiles(reshape(netinc(:, :, :, :, :, :, :, :), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), reshape(m(:, :, :, :, :, :, :, :, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(0, :, :, :, :, :, :, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), reshape(m(0, :, :, :, :, :, :, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(1, :, :, :, :, :, :, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), reshape(m(1, :, :, :, :, :, :, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(:, :, :, :, :, :, 1, :), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(:, :, :, :, :, :, 1, :, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(:, :, :, :, :, :, 2, :), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(:, :, :, :, :, :, 2, :, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(:, :, :, :, :, :, 3, :), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(:, :, :, :, :, :, 3, :, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(0, :, :, :, :, :, 1, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(0, :, :, :, :, :, 1, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(1, :, :, :, :, :, 1, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(1, :, :, :, :, :, 1, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(0, :, :, :, :, :, 2, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(0, :, :, :, :, :, 2, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(1, :, :, :, :, :, 2, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(1, :, :, :, :, :, 2, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(0, :, :, :, :, :, 3, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(0, :, :, :, :, :, 3, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))
    write(23,'(9f8.4)') percentiles(reshape(netinc(1, :, :, :, :, :, 3, :), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), reshape(m(1, :, :, :, :, :, 3, :, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*JJ/)), (/0.01d0, 0.05d0, 0.10d0, 0.25d0, 0.50d0, 0.75d0, 0.90d0, 0.95d0, 0.99d0/))

    write(23,'(/a)')'Mean (worktime)'
    write(23,'(4f8.4)') sum(reshape(l(0, :, :, :, :, :, 1, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(0, :, :, :, :, :, 1, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(0, :, :, :, :, :, 1, :JR-1, it)), &
                        sum(reshape(l(0, :, :, :, :, :, 2, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(0, :, :, :, :, :, 2, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(0, :, :, :, :, :, 2, :JR-1, it)), &
                        sum(reshape(l(0, :, :, :, :, :, 3, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(0, :, :, :, :, :, 3, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(0, :, :, :, :, :, 3, :JR-1, it)), &
                        sum(reshape(l(0, :, :, :, :, :, :, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*(JR-1)/))*reshape(m(0, :, :, :, :, :, :, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*(JR-1)/)))/sum(m(0, :, :, :, :, :, :, :JR-1, it))
    write(23,'(4f8.4)') sum(reshape(l(1, :, :, :, :, :, 1, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(1, :, :, :, :, :, 1, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(1, :, :, :, :, :, 1, :JR-1, it)), &
                        sum(reshape(l(1, :, :, :, :, :, 2, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(1, :, :, :, :, :, 2, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(1, :, :, :, :, :, 2, :JR-1, it)), &
                        sum(reshape(l(1, :, :, :, :, :, 3, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(1, :, :, :, :, :, 3, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(1, :, :, :, :, :, 3, :JR-1, it)), &
                        sum(reshape(l(1, :, :, :, :, :, :, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*(JR-1)/))*reshape(m(1, :, :, :, :, :, :, :JR-1, it), (/(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*(JR-1)/)))/sum(m(1, :, :, :, :, :, :, :JR-1, it))
    write(23,'(4f8.4)') sum(reshape(l(:, :, :, :, :, :, 1, :JR-1, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(:, :, :, :, :, :, 1, :JR-1, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(:, :, :, :, :, :, 1, :JR-1, it)), &
                        sum(reshape(l(:, :, :, :, :, :, 2, :JR-1, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(:, :, :, :, :, :, 2, :JR-1, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(:, :, :, :, :, :, 2, :JR-1, it)), &
                        sum(reshape(l(:, :, :, :, :, :, 3, :JR-1, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/))*reshape(m(:, :, :, :, :, :, 3, :JR-1, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*(JR-1)/)))/sum(m(:, :, :, :, :, :, 3, :JR-1, it)), &
                        sum(reshape(l(:, :, :, :, :, :, :, :JR-1, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*(JR-1)/))*reshape(m(:, :, :, :, :, :, :, :JR-1, it), (/2*(NA+1)*(NX+1)*(NP+1)*NW*NE*NS*(JR-1)/)))/sum(m(:, :, :, :, :, :, :, :JR-1, it))

  end subroutine


  !##############################################################################
  ! SUBROUTINE output_summary
  !
  ! Writes output summary of to file 22
  !##############################################################################
  subroutine output_summary()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: io, ia, ix, ip, iw, ie, is, ij, it
    real*8 :: HEV(-(JJ-2):TT), HEV_help, mas(-(JJ-2):0), HEVs(0:1, NS, -(JJ-2):0), mass(0:1, NS, (-JJ-2):0)

    ! aggregate ex post welfare changes of current generations
    HEV = 0d0
    mas = 0d0
    HEVs = 0d0
    mass = 0d0

    do io = 0, NO
      do ij = JJ, 2, -1
        do ia = 0, NA
          do ix = 0, NX
            do ip = 0, NP
              do is = 1, NS
                do iw = 1, NW
                  do ie = 1, NE
                    if(ij >= JR .and. ia == 0 .and. (pen(ip, ij, 0) <= 1d-10 .or. pen(ip, ij, 1) <= 1d-10))then
                      cycle
                    endif
                    HEV_help = ((VV(io, ia, ix, ip, iw, ie, is, ij, 1)/max(VV(io, ia, ix, ip, iw, ie, is, ij, 0), -1d10))**(1d0/(1d0-gamma))-1d0)*100d0
                    HEV(-(ij-2)) = HEV(-(ij-2)) + HEV_help*m(io, ia, ix, ip, iw, ie, is, ij, 1)
                    mas(-(ij-2)) = mas(-(ij-2)) + m(io, ia, ix, ip, iw, ie, is, ij, 1)
                    HEVs(io, is, -(ij-2)) = HEVs(io, is, -(ij-2)) + HEV_help*m(io, ia, ix, ip, iw, ie, is, ij, 1)
                    mass(io, is, -(ij-2)) = mass(io, is, -(ij-2)) + m(io, ia, ix, ip, iw, ie, is, ij, 1)
                  enddo ! ie
                enddo ! iw
              enddo ! is
            enddo ! ip
          enddo ! ix
        enddo ! ia
      enddo ! ij
    enddo ! io

    do io = 0, NO
      do is = 1, NS
        HEVs(io, is, -(JJ-2):0) = HEVs(io, is, -(JJ-2):0)/mass(io, is, -(JJ-2):0)
      enddo ! is
    enddo ! io
    HEV(-(JJ-2):0) = HEV(-(JJ-2):0)/mas

    ! calculate ex ante welfare of future generations
    do it = 1, TT
      HEV(it) = ((vv_coh(1, it)/vv_coh(1, 0))**(1d0/(1d0-gamma))-1d0)*100d0
    enddo

    ! headline
    write(22,'(/a,a)')'          A      KK      KC      KE      LC     l_bar       r       w     inc       C       I',  &
      '      YY      YC      YE    tauc    taup      PP      BQ     ent    ent1    ent2    ent3     HEV      DIFF'

    ! current generations
    do ij = -(JJ-2), 0
      write(22,'(i3,88x,3f8.2,20x,3f8.2,20x,f8.2)')ij, HEVs(0, :, ij), HEVs(1, :, ij), HEV(ij)
    enddo

    write(22,'(a)')'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'

    ! future generations
    do it = 1, TT
      write(22,'(i3,23f8.2,f10.5)')it,(/AA(it)/AA(0)-1d0, KK(it)/KK(0)-1d0, KC(it)/KC(0)-1d0, KE(it)/KE(0)-1d0, LC(it)/LC(0)-1d0, &
        HH(it)/sum(pop_w(:, it))/(HH(0)/sum(pop_w(:, 0)))-1d0, &
        (1d0+r(it))**0.2d0-(1d0+r(0))**0.2d0, w(it)/w(0)-1d0, inc_bar(it)/inc_bar(0)-1d0, CC(it)/CC(0)-1d0, &
        II(it)/II(0)-1d0, YY(it)/YY(0)-1d0, YC(it)/YC(0)-1d0, YE(it)/YE(0)-1d0, tauc(it)-tauc(0), taup(it)-taup(0), PP(it)/YY(it)-PP(0)/YY(0), &
        BQ(it)/BQ(0)-1d0, sum(pop_e(:, it))/sum(pop(:JR-1, it))-sum(pop_e(:, 0))/sum(pop(:JR-1, 0)), &
        pop_e(1, it)/sum(pop(:JR-1, it))-pop_e(1, 0)/sum(pop(:JR-1, 0)), &
        pop_e(2, it)/sum(pop(:JR-1, it))-pop_e(2, 0)/sum(pop(:JR-1, 0)), &
        pop_e(3, it)/sum(pop(:JR-1, it))-pop_e(3, 0)/sum(pop(:JR-1, 0))/)*100d0, &
        HEV(it), DIFF(it)/YY(it)*100d0
    enddo

    if(lsra_on)write(22, '(/a,f12.6)')'EFFICIENCY GAIN: ', (Vstar**(1d0/(1d0-gamma))-1d0)*100d0

  end subroutine


end program
