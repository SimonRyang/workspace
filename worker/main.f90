program main

  ! modules
  use omp_lib
	use accel_lib
  use globals
  use clock
  use toolbox

  implicit none

  integer, parameter :: numthreads = 28
  integer :: it

  ! set up how many kernels you want to use
  call omp_set_num_threads(numthreads)

  ! household preference parameters
  gamma = 0.500d0
  sigma = 0.318d0
  beta  = 0.998d0
  ! convert variables into per period values
  beta = beta**5d0

  ! production parameters
  alpha = 0.36d0
  delta = 0.06d0
  ! convert variables into per period values
  delta = 1d0 - (1d0-delta)**5d0

  ! demographic parameters
  n_p   = 0.0064d0
  ! convert variables into per period values
  n_p = (1d0+n_p)**5-1d0

  ! tax and transfers
  taur = 0.25d0
  tauy = 0.15d0
  tauw = 0.15d0
  kappa = 0.55d0
  mu = 1d0
  lambda = 0d0
  gy = 0.19d0
  by = 0.80d0
  ! convert variables into per period values
  by = by/5d0

  ! size of the asset grid
  a_l  = 0d0
  a_u  = 64d0
  a_grow = 0.5d0

  ! size of the pension claim grid
  p_l  = 0d0
  p_u  = 2d0
  p_grow = 0.05d0

  ! simulation parameters
  damp  = 0.50d0
  tol   = 1d-5
  itermax = 200

  ! turn gini calculation on
  gini_on = .true.

  ! calculate initial equilibrium
  call get_SteadyState()

	stop

  ! set reform parameters
  !pen_debt = .true.
  !smopec = .true.
  !phi(1:TT) = 1d0
  !lambda(1:TT) = 1d0
  !mu(1:TT) = 0d0

  !taur(1:TT) = 0.20d0

  ! calculate transition path without lsra
  lsra_on = .false.
  call get_transition()

  stop

  ! calculate transition path with lsra
  lsra_on = .true.
  call get_transition()

  ! close files
  close(21)
  close(22)

contains


  !##############################################################################
  ! SUBROUTINE get_SteadyState
  !
  ! Computes the initial steady state of the economy
  !##############################################################################
  subroutine get_SteadyState()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: iter

    ! initialize remaining variables
    call initialize()

    ! start timer
    call tick(calc)

    ! iterate until value function converges
    do iter = 1, itermax

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

      write(*,'(i4,5f8.2,f12.6)')iter, (/5d0*KK(0), CC(0), II(0)/)/YY(0)*100d0, &
        ((1d0+r(0))**0.2d0-1d0)*100d0, w(0), DIFF(0)/YY(0)*100d0

      if(abs(DIFF(0)/YY(0))*100d0 < tol)then

        call tock(calc)
        call output(0)
        return
      endif

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
      write(*,'(a)')'ITER  COMP_OLD  EFFICIENCY    DIFF'
    endif

    ! start timer
    call tic()

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
        write(*,'(i4,5f8.2,f12.6)')iter, (/5d0*KK(TT), CC(TT), II(TT)/)/YY(TT)*100d0, &
          ((1d0+r(TT))**0.2d0-1d0)*100d0, w(TT), DIFF(itmax)/YY(itmax)*100d0
        check = abs(DIFF(itmax)/YY(itmax))*100d0 < tol .and. iter > 1
      else
        write(*,'(i4,3f12.5)')iter, lsra_comp/lsra_all*100d0, &
          (Vstar**(1d0/(1d0-1d0/gamma))-1d0)*100d0,DIFF(itmax)/YY(itmax)*100d0
          check = abs(DIFF(itmax)/YY(itmax))*100d0 < tol .and. iter > 1 .and. lsra_comp/lsra_all > 0.99999d0
      endif

      ! check for convergence
      if(check)then
        call toc
        do it = 1, TT
          if (.not. lsra_on) call output(it)
        enddo
        call output_summary()
        return
      endif
    enddo

    call toc
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

    !##### OTHER VARIABLES ####################################################
    integer :: ix, ip, is, iw, ij
    real*8 :: adj

    write(*,'(/a/)')'INITIAL EQUILIBRIUM'
    write(*,'(a)')'ITER   K/Y   C/Y   I/Y     r     w    DIFF'

    ! initialize asset grid
    a = grid_Cons_Grow(a_l, a_u, a_grow, NA)

    ! initialize pension claim grid
    if (p_grow > 0d0) p = grid_Cons_Grow(p_l, p_u, p_grow, NP)
    if (p_grow <= 0d0) p = grid_Cons_Equi(p_l, p_u, NP)

    ! get initial guess for savings decision
    do ij = 1, JJ
      do iw = 1, NW
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NA
              aplus(:, ix, ip, is, iw, ij, 0) = max(a(:)/2d0, a(1))
            enddo
          enddo
        enddo
      enddo
    enddo

    ! initial guess for labor decision
    l(:, :, :, :, :, :, 0) = 0.33d0

    ! initial guess for annuity investment
    mx(:, :, :, :, :, :, 0) = 0.33d0

    ! distribution of skill classes
    dist_skill = (/0.263d0, 0.545d0, 0.192d0/)

		! initialize survival probabilities for middle skilled
    open(301, file='sp.dat')
    do ij = 1, JJ+1
      read(301,'(f13.8)')psi(ij, 2)
    enddo
    close(301)

    ! compute survival probabilities for high/low skilled
    psi(1, :) = psi(1, 2)
    psi(JJ+1, :) = 0d0
    adj = 22d0
    do ij = 2, JJ
      psi(ij, 1) = psi(ij, 2) - exp(0.33d0*(dble(ij-1)-adj))
      psi(ij, 3) = psi(ij, 2) + exp(0.33d0*(dble(ij-1)-adj))
    enddo

    !psi(10:, 1) = psi(10:, 2)
    !psi(10:, 3) = psi(10:, 2)

    ! set up population structure
    rpop(1, :, 0) = 1d0
    do ij = 2, JJ
      rpop(ij, :, 0) = psi(ij, :)*rpop(ij-1, :, 0)/(1d0+n_p)
    enddo

    pop = 0d0
    do is = 1, NS
      pop(:, 0) = pop(:, 0) + rpop(:, is, 0)*dist_skill(is)
    enddo

    ! calculate average survival probabilities
    psi_avg(1) = 1d0
    do ij = 1, JJ
      psi_avg(ij+1) = sum(rpop(ij, :, 0)*dist_skill(:)*psi(ij+1, :))/pop(ij, 0)
    enddo

    ! set distribution of bequests
    Gama(1:JR-1) = pop(1:JR-1, 0)/sum(pop(1:JR-1, 0))
    Gama(JR:JJ) = 0d0

    ! initialize age earnings process
    open(302, file='eff.dat')
    do ij = 1, JJ
      read(302,'(3f12.8)')eff(ij, :)
    enddo
    close(302)
    eff(JR:, :) = 0d0

    call discretize_AR(0.95666d0**5d0, 0.0d0, sigma5(0.95666d0, 0.02321d0), eta(1, :), pi_eta(1, :, :), dist_eta(1, :))
    eta(1, :) = exp(eta(1, :))/sum(dist_eta(1,:)*exp(eta(1, :)))

    call discretize_AR(0.95687d0**5d0, 0.0d0, sigma5(0.95687d0, 0.02812d0), eta(2, :), pi_eta(2, :, :), dist_eta(2, :))
    eta(2, :) = exp(eta(2, :))/sum(dist_eta(2,:)*exp(eta(2, :)))

    call discretize_AR(0.95828d0**5d0, 0.0d0, sigma5(0.95828d0, 0.03538d0), eta(3, :), pi_eta(3, :, :), dist_eta(3, :))
    eta(3, :) = exp(eta(3, :))/sum(dist_eta(3,:)*exp(eta(3, :)))

    ! initial guesses for macro variables
    tauc(0) = 0.05d0
    taup(0)  = 0.20d0
    inc_bar(0) = 0.40d0
    BQ(0) = 0.10d0
    YY(0) = 5.00d0
    BB(0) = by*YY(0)
    GG(0) = gy*YY(0)
    KC(0) = 3.80d0
    LC(0) = 5.50d0

    ! open files
    open(21, file='output.out')
    open(22, file='summary.out')

  end subroutine


  !##############################################################################
  ! SUBROUTINE initialize_trn
  !
  ! Initializes transitional variables
  !##############################################################################
  subroutine initialize_trn()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: ij, is, it

    write(*,'(/a/)')'TRANSITION PATH'
    write(*,'(a)')'ITER   K/Y   C/Y   I/Y     r     w    DIFF'

    do it = 1, TT

      taup(it) = taup(0)
      tauc(it) = tauc(0)

      TAc(it) = TAc(0)
      TAr(it) = TAr(0)
      TAw(it) = TAw(0)
      TAy(it) = TAy(0)

      pop_w(:, it) = pop_w(:, 0)
      pop_r(:, it) = pop_r(:, 0)

      r(it) = r(0)
      w(it) = w(0)
      inc_bar(it) = inc_bar(0)
      px(:, it) = px(:, 0)
      beq(:, it) = beq(:, 0)
      pen(:, :, it) = pen(:, :, 0)

      KK(it) = KK(0)
      KC(it) = KC(0)
      AA(it) = AA(0)
      BB(it) = BB(0)
      HH(it) = HH(0)
      LC(it) = LC(0)
      YY(it) = YY(0)
      YC(it) = YC(0)
      CC(it) = CC(0)
      II(it) = II(0)
      GG(it) = GG(0)
      BQ(it) = BQ(0)
      PP(it) = PP(0)
      PC(it) = PC(0)
      BF(it) = BF(0)
      NX(it) = NX(0)

      rpop(:, :, it) = rpop(:, :, 0)
      pop(:, it) = pop(:, 0)

      c_coh(:, it) = c_coh(:, 0)
      a_coh(:, it) = a_coh(:, 0)
      x_coh(:, it) = x_coh(:, 0)
      bx_coh(:, it) = bx_coh(:, 0)
      inc_coh(:, it) = inc_coh(:, 0)
      flc_coh(:, it) = flc_coh(:, 0)
      vv_coh(:, it) = vv_coh(:, 0)

      aplus(:, :, :,  :, :, :, it) = aplus(:, :, :, :, :, :, 0)
      xplus(:, :, :, :, :, :, it) = xplus(:, :, :, :, :, :, 0)
      mx(:, :, :, :, :, :, it) = mx(:, :, :, :, :, :, 0)
      pplus(:, :, :, :, :, :, it) = pplus(:, :, :, :, :, :, 0)
      c(:, :, :, :, :, :, it) = c(:, :, :, :, :, :, 0)
      l(:, :, :, :, :, :, it) = l(:, :, :, :, :, :, 0)
      m(:, :, :, :, :, :, it) = m(:, :, :, :, :, :, 0)
      pencon(:, :, :, :, :, :, it) = pencon(:, :, :, :, :, :, 0)
      inctax(:, :, :, :, :, :, it) = inctax(:, :, :, :, :, :, 0)
      captax(:, :, :, :, :, :, it) = captax(:, :, :, :, :, :, 0)
      VV(:, :, :, :, :, :, it) = VV(:, :, :, :, :, :, 0)
      EV(:, :, :, :, :, :, it) = EV(:, :, :, :, :, :, 0)

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

    ! calculate factor prices
    if (.not. smopec) then
      r(it) = (1d0-tauy(it))*(alpha*(KC(it)/LC(it))**(alpha-1d0)-delta)
      w(it) = (1d0-alpha)*(KC(it)/LC(it))**alpha
    endif

    ! calculate individual bequests
    beq(:, it) = Gama*BQ(it)/pop(:, it)

    ! calculate individual pensions
    pen(:, :, it) = 0d0
    do ij = JR, JJ
      pen(ij, :, it) = p(:)*kappa(it)*inc_bar(it)
    enddo

    ! calculate interests of annuities
    do ij = 1, JJ
      px(ij, it) = 1d0
      if (x_coh(ij, it) > 0d0) px(ij, it) = x_coh(ij, it)/(x_coh(ij, it) + bx_coh(ij, it))
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
    integer ::ia, ix, ip, is, iw, ij, it
		real*8 :: x(3), fret

		!$acc routine(fminsearch)

    do ij = JJ, 1, -1

      it = year(it_in, ij_in, ij)

      if (ij >= JR) then
				
				ij_com = ij
				it_com = it
				iw_com = 1

  	  	!$omp parallel do copyin(ij_com, it_com, iw_com) collapse(3) schedule(dynamic,1) private(x, fret) num_threads(numthreads)
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NA
              do ia = 0, NA

							! set up communication variables
							ia_com = ia
							ix_com = ix
							ip_com = ip
							is_com = is

							! get initial guess for the individual choices
							x(1) = max(aplus(ia, ix, ip, is, 1, ij, it), 1d-4)
							x(2) = max(l(ia, ix, ip, is, 1, ij, it), 1d-4)
							x(3) = max(mx(ia, ix, ip, is, 1, ij, it), 1d-4)

							call fminsearch(x, fret, (/a_l, 0d0, a_l/), (/a_u, 1d0, a_u/), valuefunc_w)

							! copy decisions
							aplus(ia, ix, ip, is, :, ij, it) = x(1)
							xplus(ia, ix, ip, is, :, ij, it) = xplus_com
							pplus(ia, ix, ip, is, :, ij, it) = pplus_com
							c(ia, ix, ip, is, :, ij, it) = max(c_com, 1d-16)
							l(ia, ix, ip, is, :, ij, it) = l_com
							mx(ia, ix, ip, is, :, ij, it) = mx_com
							pencon(ia, ix, ip, is, :, ij, it) = pencon_com
							inctax(ia, ix, ip, is, :, ij, it) = inctax_com
							captax(ia, ix, ip, is, :, ij, it) = captax_com
							VV(ia, ix, ip, is, :, ij, it) = -fret

              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
				!$omp end parallel do

      elseif (ij >= 2) then

				ij_com = ij
				it_com = it

				!$acc parallel loop
  	  	!!$omp parallel do copyin(ij_com, it_com) collapse(3) schedule(dynamic,1) private(x, fret) num_threads(numthreads)
        do iw = 1, NW
          do is = 1, NS
            do ip = 0, NP
              do ix = 0, NA
                do ia = 0, NA

									! set up communication variables
									ia_com = ia
									ix_com = ix
									ip_com = ip
									is_com = is
									iw_com = iw
									ij_com = ij
									it_com = it

									! get initial guess for the individual choices
									x(1) = max(aplus(ia, ix, ip, is, iw, ij, it), 1d-4)
									x(2) = max(l(ia, ix, ip, is, iw, ij, it), 1d-4)
									x(3) = max(mx(ia, ix, ip, is, iw, ij, it), 1d-4)

									call fminsearch(x, fret, (/a_l, 0d0, a_l/), (/a_u, 1d0, a_u/), valuefunc_w)

									! copy decisions
									aplus(ia, ix, ip, is, iw, ij, it) = x(1)
									xplus(ia, ix, ip, is, iw, ij, it) = xplus_com
									pplus(ia, ix, ip, is, iw, ij, it) = pplus_com
									c(ia, ix, ip, is, iw, ij, it) = max(c_com, 1d-16)
									l(ia, ix, ip, is, iw, ij, it) = l_com
									mx(ia, ix, ip, is, iw, ij, it) = mx_com
									pencon(ia, ix, ip, is, iw, ij, it) = pencon_com
									inctax(ia, ix, ip, is, iw, ij, it) = inctax_com
									captax(ia, ix, ip, is, iw, ij, it) = captax_com
									VV(ia, ix, ip, is, iw, ij, it) = -fret

                enddo ! ia
              enddo ! ix
            enddo ! ip
          enddo ! is
        enddo ! iw
				!!$omp end parallel do
				!$acc end parallel loop

      elseif (ij == 1) then

				ij_com = ij
				it_com = it
				ia_com = 0
				ix_com = 0
				ip_com = 0

  	  	!$omp parallel do copyin(ij_com, it_com, ia_com, ix_com, ip_com) collapse(2) schedule(dynamic,1) private(x, fret) num_threads(numthreads)
        do iw = 1, NW
          do is = 1, NS

						! set up communication variables
						is_com = is
						iw_com = iw

						! get initial guess for the individual choices
						x(1) = max(aplus(0, 0, 0, is, iw, ij, it), 1d-4)
						x(2) = max(l(0, 0, 0, is, iw, ij, it), 1d-4)
						x(3) = max(mx(0, 0, 0, is, iw, ij, it), 1d-4)

						call fminsearch(x, fret, (/a_l, 0d0, a_l/), (/a_u, 1d0, a_u/), valuefunc_w)

            ! copy decisions
            aplus(:, :, :, is, iw, ij, it) = x(1)
            xplus(:, :, :, is, iw, ij, it) = xplus_com
            pplus(:, :, :, is, iw, ij, it) = pplus_com
            c(:, :, :, is, iw, ij, it) = max(c_com, 1d-16)
            l(:, :, :, is, iw, ij, it) = l_com
            mx(:, :, :, is, iw, ij, it) = mx_com
            pencon(:, :, :, is, iw, ij, it) = pencon_com
            inctax(:, :, :, is, iw, ij, it) = inctax_com
            captax(:, :, :, is, iw, ij, it) = captax_com
            VV(:, :, :, is, iw, ij, it) = -fret

          enddo ! is
        enddo ! iw
				!$omp end parallel do

      endif

      ! interpolate individual value function
      call interpolate(ij, it)

    enddo

  end subroutine



  !##############################################################################
  ! SUBROUTINE get_decision
  !
  ! Calculates the optimal decision of an individual
  !##############################################################################
  subroutine get_decision(ia, ix, ip, is, iw, ij, it)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: ia, ix, ip, is, iw, ij, it

    !##### OTHER VARIABLES ####################################################
    real*8 :: x(3), fret

    ! set up communication variables
    ia_com = ia
    ix_com = ix
    ip_com = ip
    is_com = is
    iw_com = iw
    ij_com = ij
    it_com = it

    ! get initial guess for the individual choices
    x(1) = max(aplus(ia, ix, ip, is, iw, ij, it), 1d-4)
    x(2) = max(l(ia, ix, ip, is, iw, ij, it), 1d-4)
    x(3) = max(mx(ia, ix, ip, is, iw, ij, it), 1d-4)

    call fminsearch(x, fret, (/a_l, 0d0, a_l/), (/a_u, 1d0, a_u/), valuefunc_w)

    ! copy decisions
    aplus(ia, ix, ip, is, iw, ij, it) = x(1)
    xplus(ia, ix, ip, is, iw, ij, it) = xplus_com
    pplus(ia, ix, ip, is, iw, ij, it) = pplus_com
    c(ia, ix, ip, is, iw, ij, it) = max(c_com, 1d-16)
    l(ia, ix, ip, is, iw, ij, it) = l_com
    mx(ia, ix, ip, is, iw, ij, it) = mx_com
    pencon(ia, ix, ip, is, iw, ij, it) = pencon_com
    inctax(ia, ix, ip, is, iw, ij, it) = inctax_com
    captax(ia, ix, ip, is, iw, ij, it) = captax_com
    VV(ia, ix, ip, is, iw, ij, it) = -fret

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
    integer :: ia, ix, ip, is, iw

    do iw = 1, NW
      do is = 1, NS
        do ip = 0, NP
          do ix = 0, NA
            do ia = 0, NA

              EV(ia, ix, ip, is, iw, ij, it) &
                = (1d0-1d0/gamma)*sum(pi_eta(is, iw, :)*VV(ia, ix, ip, is, :, ij, it))**(1d0/(1d0-1d0/gamma))

            enddo ! ia
          enddo ! ix
        enddo ! ip
      enddo ! is
    enddo ! iw

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
    integer :: ia, ix, ip, is, iw, ij, ial, iar, ixl, ixr, ipl, ipr, itm
    real*8 :: varphi, varchi, varpsi

    ! get yesterdays year
    itm = year(it, 2, 1)

    ! set distribution to zero
    m(:, :, :, :, :, :, it) = 0d0

    ! get initial distribution at age 1
    do iw = 1, NW
      do is = 1, NS
        m(0, 0, 0, is, iw, 1, it) = dist_skill(is)*dist_eta(is, iw)
      enddo ! is
    enddo ! iw

    ! successively compute distribution over ages
    do ij = 2, JJ

      ! iterate over yesterdays gridpoints
      do iw = 1, NW
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NA
              do ia = 0, NA

                ! interpolate yesterday's savings decision
                call linint_Grow(aplus(ia, ix, ip, is, iw, ij-1, itm), &
                          a_l, a_u, a_grow, NA, ial, iar, varphi)
                call linint_Grow(xplus(ia, ix, ip, is, iw, ij-1, itm), &
                          a_l, a_u, a_grow, NA, ixl, ixr, varchi)

                ! interpolate today's pension claims
                if (p_grow > 0d0) call linint_Grow(pplus(ia, ix, ip, is, iw, ij-1, itm), &
                                  p_l, p_u, p_grow, NP, ipl, ipr, varpsi)
                if (p_grow <= 0d0) call linint_Equi(pplus(ia, ix, ip, is, iw, ij-1, itm), &
                                  p_l, p_u, NP, ipl, ipr, varpsi)

                ! redistribute households
                m(ial, ixl, ipl, is, :, ij, it) = &
                  m(ial, ixl, ipl, is, :, ij, it) &
                  +varphi*varchi*varpsi*psi(ij, is)*pi_eta(is, iw, :) &
                  *m(ia, ix, ip, is, iw, ij-1, itm)
                m(ial, ixl, ipr, is, :, ij, it) = &
                  m(ial, ixl, ipr, is, :, ij, it) &
                  +varphi*varchi*(1d0-varpsi)*psi(ij, is)*pi_eta(is, iw, :) &
                  *m(ia, ix, ip, is, iw, ij-1, itm)
                m(ial, ixr, ipl, is, :, ij, it) = &
                  m(ial, ixr, ipl, is, :, ij, it) &
                  +varphi*(1d0-varchi)*varpsi*psi(ij, is)*pi_eta(is, iw, :) &
                  *m(ia, ix, ip, is, iw, ij-1, itm)
                m(ial, ixr, ipr, is, :, ij, it) = &
                  m(ial, ixr, ipr, is, :, ij, it) &
                  +varphi*(1d0-varchi)*(1d0-varpsi)*psi(ij, is)*pi_eta(is, iw, :) &
                  *m(ia, ix, ip, is, iw, ij-1, itm)
                m(iar, ixl, ipl, is, :, ij, it) = &
                  m(iar, ixl, ipl, is, :, ij, it) &
                  +(1d0-varphi)*varchi*varpsi*psi(ij, is)*pi_eta(is, iw, :) &
                  *m(ia, ix, ip, is, iw, ij-1, itm)
                m(iar, ixl, ipr, is, :, ij, it) = &
                  m(iar, ixl, ipr, is, :, ij, it) &
                  +(1d0-varphi)*varchi*(1d0-varpsi)*psi(ij, is)*pi_eta(is, iw, :) &
                  *m(ia, ix, ip, is, iw, ij-1, itm)
                m(iar, ixr, ipl, is, :, ij, it) = &
                  m(iar, ixr, ipl, is, :, ij, it) &
                  +(1d0-varphi)*(1d0-varchi)*varpsi*psi(ij, is)*pi_eta(is, iw, :) &
                  *m(ia, ix, ip, is, iw, ij-1, itm)
                m(iar, ixr, ipr, is, :, ij, it) = &
                  m(iar, ixr, ipr, is, :, ij, it) &
                  +(1d0-varphi)*(1d0-varchi)*(1d0-varpsi)*psi(ij, is)*pi_eta(is, iw, :) &
                  *m(ia, ix, ip, is, iw, ij-1, itm)

              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
      enddo ! iw

      m(:, :, :, :, :, ij, it) = m(:, :, :, :, :, ij, it)/sum(m(:, :, :, :, :, ij, it))

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
    integer :: ij, ia, ix, ip, is, iw, itp
    real*8 :: LC_old

    LC_old = LC(it)

    ! get tomorrow's year
    itp = year(it, 1, 2)

    ! calculate aggregates
    AA(it)  = 0d0
    CC(it)  = 0d0
    LC(it)  = 0d0
    HH(it)  = 0d0
    PC(it)  = 0d0
    BQ(it)  = 0d0
    PP(it)  = 0d0
    TAc(it)   = 0d0
    TAr(it)   = 0d0
    TAw(it)   = 0d0
    TAy(it)   = 0d0
    pop_w(:, it) = 0d0
    pop_r(:, it) = 0d0
    x_coh(:, it) = 0d0
    bx_coh(:, it) = 0d0
    vv_coh(:, it) = 0d0

    do ij = 1, JJ
      do iw = 1, NW
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NA
              do ia = 0, NA

                AA(it) = AA(it) + a(ia)*m(ia, ix, ip, is, iw, ij, it)/psi(ij, is)*pop(ij, it)
                AA(it) = AA(it) + a(ix)*m(ia, ix, ip, is, iw, ij, it)/psi(ij, is)*pop(ij, it)
                x_coh(ij, it) = x_coh(ij, it) + a(ix)*m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                bx_coh(ij, it) = bx_coh(ij, it) + a(ix)*(1d0-psi(ij, is)) &
                          *m(ia, ix, ip, is, iw, ij, it)/psi(ij, is)*pop(ij, it)
                CC(it) = CC(it) + c(ia, ix, ip, is, iw, ij, it) &
                          *m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                LC(it) = LC(it) + eff(ij, is)*eta(is, iw)*l(ia, ix, ip, is, iw, ij, it) &
                          *m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                HH(it) = HH(it) + l(ia, ix, ip, is, iw, ij, it) &
                          *m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                PC(it) = PC(it) + min(w(it)*eff(ij, is)*eta(is, iw)*l(ia, ix, ip, is, iw, ij, it), 2d0*inc_bar(it)) &
                          *m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                PP(it) = PP(it) + pen(ij, ip, it)*m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                BQ(it) = BQ(it) + (1d0+r(it))*a(ia)*(1d0-psi(ij, is)) &
                          *m(ia, ix, ip, is, iw, ij, it)/psi(ij, is)*pop(ij, it)
                TAc(it) = TAc(it) + tauc(it)*c(ia, ix, ip, is, iw, ij, it) &
                          *m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                TAr(it) = TAr(it) + captax(ia, ix, ip, is, iw, ij, it)*m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                TAw(it) = TAw(it) + inctax(ia, ix, ip, is, iw, ij, it)*m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                if (ij < JR) then
                  pop_w(is, it) = pop_w(is, it) + m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                else
                  pop_r(is, it) = pop_r(is, it) + m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                endif
                vv_coh(ij, it) = vv_coh(ij, it) + VV(ia, ix, ip, is, iw, ij, it) &
                                *m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)


              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
      enddo ! iw
    enddo ! ij

    ! damping and other quantities
    LC(it) = damp*LC(it) + (1d0-damp)*LC_old

    if (smopec) then
      KC(it) = damp*(LC(it)*((r(it)/(1d0-tauy(it))+delta)/alpha)**(1d0/(alpha-1d0)))+(1d0-damp)*KC(it)
      BF(it) = AA(it) - KC(it) - BB(it) - BP(it) - BA(it)
      NX(it) = (n_p-r(it))*BF(it)
    else
      KC(it) = damp*(AA(it) - BB(it) - BP(it) - BA(it))+(1d0-damp)*KC(it)
      BF(it) = 0d0
      NX(it) = 0d0
    endif

    KK(it) = KC(it)

    II(it) = (1d0+n_p)*KK(itp) - (1d0-delta)*KK(it)
    YC(it) = KC(it)**alpha * LC(it)**(1d0-alpha)
    YY(it) = YC(it)

    TAy(it) = TAy(it) + tauy(it)*(YC(it) - delta*KC(it) - w(it)*LC(it))

    inc_bar(it) = w(it)*LC(it)/sum(pop_w(:, it))

    ! get difference on goods market
    DIFF = YY(it)-CC(it)-II(it)-GG(it)-NX(it)

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
      !write(*,*)taup(0), PP_bar, PC_bar, taup(1)

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
    integer :: ij, ia, ix, ip, is, iw, it
    real*8 :: VV_today, VV_target, dVV_da, v_tilde
    real*8 :: pv_today, pv_target, pv_trans

    ! initialize variables
    SV(:) = 0d0

    ! initialize counters
    lsra_comp   = 0d0
    lsra_all    = 0d0

    do ij = 2, JJ
      do iw = 1, NW
        do is = 1, NS
          do ip = 0, NP
            do ia = 0, NA

              ! do not do anything for an agent at retirement without pension and savings
              if(ij >= JR .and. ia == 0 .and. (pen(ij, ip, 0) <= 1d-10 .or. pen(ij, ip, 1) <= 1d-10))then
                v(ia, ix, ip, is, iw, ij, 1) = 0d0
                cycle
              endif

              ! get today's utility
              VV_today = VV(ia, ix, ip, is, iw, ij, 1)
              if (VV_today >= 0d0)write(*,*)'VV_today', VV_today
              if (VV_today <= -1d10)write(*,*)'VV_today', VV_today

              ! get target utility
              VV_target = VV(ij, ia, ix, ip, is, iw, 0)
              if (VV_target >= 0d0)write(*,*)'VV_target', VV_target
              if (VV_target <= -1d10)write(*,*)'VV_target', VV_target

              ! get derivative of the value function
              dVV_da = margu(c(ia, ix, ip, is, iw, ij, 1),l(ia, ix, ip, is, iw, ij, 1), 1)
              if (dVV_da < 0d0)write(*,*)'dVV_da', dVV_da

              ! calculate change in transfers
              v_tilde = (VV_target-VV_today)/dVV_da

              ! check whether individual is already compensated
              lsra_all = lsra_all + m(ia, ix, ip, is, iw, ij, 1)*pop(ij, 1)
              if(abs((VV_today-VV_target)/VV_target) < tol) &
                lsra_comp = lsra_comp + m(ia, ix, ip, is, iw, ij, 1)*pop(ij, 1)

              ! calculate total transfer
              v(ia, ix, ip, is, iw, ij, 1) = v(ia, ix, ip, is, iw, ij, 1) + damp*damp*v_tilde

              ! aggregate transfers by cohort
              SV(1) = SV(1) + v(ia, ix, ip, is, iw, ij, 1)*m(ia, ix, ip, is, iw, ij, 1)*pop(ij, 1)

            enddo ! is
          enddo ! ip
        enddo ! is
      enddo ! iw
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
      do is = 1, NS
        do iw = 1, NW
          dVV_da = dVV_da + margu(c(0, 0, 0, is, iw, 1, it), l(0, 0, 0, is, iw, 1, it), it)*m(0, 0, 0, is, iw, 1, it)
        enddo
      enddo

      ! calculate present values
      if(it == TT)then
        pv_today  = VV_today/dVV_da  *(1d0+r(it))/(r(it)-n_p)
        pv_target = VV_target/dVV_da   *(1d0+r(it))/(r(it)-n_p)
        pv_trans  = v(0, 0, 0, 1, 1, 1, it)*(1d0+r(it))/(r(it)-n_p)
      else
        pv_today  = pv_today *(1d0+n_p)/(1d0+r(it+1)) + VV_today/dVV_da
        pv_target = pv_target*(1d0+n_p)/(1d0+r(it+1)) + VV_target/dVV_da
        pv_trans  = pv_trans *(1d0+n_p)/(1d0+r(it+1)) + v(0, 0, 0, is, iw, 1, it)
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
      do is = 1, NS
        do iw = 1, NW
          dVV_da = dVV_da + margu(c(0, 0, 0, is, iw, 1, it), l(0, 0, 0, is, iw, 1, it), it)*m(0, 0, 0, is, iw, 1, it)
        enddo ! iw
      enddo ! is

      ! compute change in transfers (restricted)
      v_tilde = (VV_target-VV_today)/dVV_da

      ! calculate cohort transfer level
      v(0, 0, 0, :, :, 1, it) = v(0, 0, 0, :, :, 1, it) + v_tilde

      ! aggregate transfers
      SV(it) = SV(it) + v(0, 0, 0, is, iw, 1, it)

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
    integer :: ij, ia, ix, ip, is, iw, iamax(JJ), ixmax(JJ)
    real*8 :: life_exp(NS), punb(JJ, NS)

    life_exp = 0d0
    do is = 1, NS
      punb(1, is) = psi(1, is)
      life_exp(is) = life_exp(is) + 2.5d0*dble(1)*punb(1,is)*(1d0-psi(2,is))
      do ij = 2, JJ
        punb(ij, is) = punb(ij-1, is)*psi(ij, is)
        life_exp(is) = life_exp(is) + (2.5d0 + 5d0*dble(ij-1))*punb(ij, is)*(1d0-psi(ij+1, is))
      enddo ! ij
    enddo ! is

    c_coh(:, it) = 0d0
    a_coh(:, it) = 0d0
    x_coh(:, it) = 0d0
    inc_coh(:, it) = 0d0
    flc_coh(:, it) = 0d0

    do ij = 1, JJ
      do ia = 0, NA
        do ix = 0, NA
          do ip = 0, NP
            do is = 1, NS
              do iw = 1, NW
                c_coh(ij, it) = c_coh(ij, it) + c(ia, ix, ip, is, iw, ij, it) &
                                    *m(ia, ix, ip, is, iw, ij, it)
                a_coh(ij, it) = a_coh(ij, it) + a(ia)*m(ia, ix, ip, is, iw, ij, it)
                x_coh(ij, it) = x_coh(ij, it) + a(ix)*m(ia, ix, ip, is, iw, ij, it)
                inc_coh(ij, it) = inc_coh(ij, it) + (w(it)*eff(ij, is)*eta(is, iw)*l(ia, ix, ip, is, iw, ij, it) &
                                   + pen(ij, ip, it))*m(ia, ix, ip, is, iw, ij, it)
                if (aplus(ia, ix, ip, is, iw, ij, it) <= 1d-10) then
                  flc_coh(ij, it) = flc_coh(ij, it) + m(ia, ix, ip, is, iw, ij, it)
                endif
              enddo ! iw
            enddo ! is
          enddo ! ip
        enddo ! ix
      enddo ! ia
    enddo ! ij

    do ij = 1, JJ
      c_coh(ij, it) = c_coh(ij, it)
      a_coh(ij, it) = a_coh(ij, it)
      inc_coh(ij, it) = inc_coh(ij, it)
    enddo ! ij

    ! Output
    write(21,'(a, i3/)')'EQUILIBRIUM YEAR ', it
    write(21,'(a)')'CAPITAL     KK    KC    AA     r  p.a.'
    write(21,'(8x,5f8.2)')KK(it), KC(it), AA(it), r(it), ((1d0+r(it))**(1d0/5d0)-1d0)*100d0
    write(21,'(a,3f8.2/)')'(in %)  ',(/KK(it), KC(it), AA(it)/)/YY(it)*500d0

    write(21,'(a)')'LABOR     LC     w   inc   l_bar'
    write(21,'(8x,4f8.2/)')LC(it), w(it), inc_bar(it), HH(it)/sum(pop_w(:, it))

    write(21,'(a)')'GOODS     YY    YC    CC    II    NX  DIFF'
    write(21,'(8x,5f8.2,f8.3)')YY(it), YC(it), CC(it), II(it), NX(it), diff(it)
    write(21,'(a,5f8.2,f8.3/)')'(in %)  ',(/YY(it), YC(it), CC(it), II(it), NX(it), diff(it)/)/YY(it)*100d0

    write(21,'(a)')'GOV     TAUC  TAUR  TAUW  TAUY   TOTAL    GG    BB    BF'
    write(21,'(8x,8f8.2)')TAc(it), TAr(it), TAw(it), TAy(it), TAc(it)+TAr(it)+TAw(it)+TAy(it), GG(it), BB(it), BF(it)
    write(21,'(a,8f8.2)')'(in %)  ',(/TAc(it), TAr(it), TAw(it), TAy(it), TAc(it)+TAr(it)+TAw(it)+TAy(it), &
               GG(it), BB(it)*5d0, BF(it)*5d0/)/YY(it)*100d0
    write(21,'(a,4f8.2/)')'(rate)  ',(/tauc(it), taur(it), tauw(it), tauy(it)/)*100d0

    write(21,'(a)')'PENS    TAUP   PEN    PP    PC    BQ'
    write(21,'(8x,5f8.2)')taup(it)*PC(it), PP(it)/sum(pop(JR:, it)), PP(it), PC(it), BQ(it)
    write(21,'(a,5f8.2/)')'(in %)  ',(/taup(it), kappa(it), PP(it)/YY(it), PC(it)/YY(it), BQ(it)/YY(it)/)*100d0

    write(21,'(a)')'POP    TOTAL   65-   65+ 65+/65-'
    write(21,'(8x,3f8.2)')sum(pop(:, it)), sum(pop(:JR-1, it)), sum(pop(JR:, it))
    write(21,'(a,4f8.2/)')'(in %)  ',(/sum(pop(:, it)), sum(pop(:JR-1, it)), sum(pop(JR:, it))/)/sum(pop(:, it))*100d0, &
                sum(pop(JR:, it))/sum(pop(:JR-1, it))*100d0

    write(21,'(a)')'LIFE     j_bar  j_bar1  j_bar2  j_bar3'
    write(21,'(8x,4f8.2/)')sum(life_exp*dist_skill), life_exp(:)

    ! check for the maximium grid point used

    call check_grid(iamax, ixmax, it)

    write(21, '(a,a)')' IJ   CONSw  ASSw  ASSx  INCw   FLC    VALUE  IAMAX  IXMAX'
    do ij = 1, JJ
      write(21, '(i3, 5f8.3, f11.3, 2i7)')ij, c_coh(ij, it), a_coh(ij, it), x_coh(ij, it), inc_coh(ij, it), &
                         flc_coh(ij, it), vv_coh(ij, it), iamax(ij), ixmax(ij)

    enddo

    write(21,'(a/)')' '

    if (gini_on .and. (it == 0 .or. it == TT)) then

      call compute_gini(it, 0)
      call compute_gini(it, 1)

      write(21, '(a)')'WEALTH  GINI   1%   5%  10%  20%  40%  60%  FLC'
      write(21,'(4x, f10.3, 8f7.2)')gini(it, 0), percentiles(:, it, 0)*100d0, sum(pop(:, it)*flc_coh(:, it))/sum(pop(:, it))*100d0

      write(21, '(a)')'INCOME  GINI   1%   5%  10%  20%  40%  60%'
      write(21,'(4x, f10.3, 6f7.2/)')gini(it, 1), percentiles(:, it, 1)*100d0

    endif

    write(21,'(a/)')'------------------------------------------------------------------------------------------------------------'

  end subroutine


  !##############################################################################
  ! SUBROUTINE output_summary
  !
  ! Writes output summary of to file 22
  !##############################################################################
  subroutine output_summary()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: ij, ia, ix, ip, is, iw, it
    real*8 :: HEV(-(JJ-2):TT), HEV_help, mas(-(JJ-2):0), HEVs(NS, -(JJ-2):0), mass(NS, (-JJ-2):0)

    ! aggregate ex post welfare changes of current generations
    HEV = 0d0
    mas = 0d0
    HEVs = 0d0
    mass = 0d0

    do ij = JJ, 2, -1
      do ia = 0, NA
        do ip = 0, NP
          do is = 1, NS
            do iw = 1, NW
              if(ij >= JR .and. ia == 0 .and. (pen(ij, ip, 0) <= 1d-10 .or. pen(ij, ip, 1) <= 1d-10))then
                cycle
              endif
              HEV_help = ((VV(ia, ix, ip, is, iw, ij, 1) &
                          /max(VV(ij, ia, ix, ip, is, iw, 0), -1d10))**(1d0/(1d0-1d0/gamma))-1d0)*100d0
              HEV(-(ij-2)) = HEV(-(ij-2)) + HEV_help*m(ia, ix, ip, is, iw, ij, 1)
              mas(-(ij-2)) = mas(-(ij-2)) + m(ia, ix, ip, is, iw, ij, 1)
              HEVs(is, -(ij-2)) = HEVs(is, -(ij-2)) + HEV_help*m(ia, ix, ip, is, iw, ij, 1)
              mass(is, -(ij-2)) = mass(is, -(ij-2)) + m(ia, ix, ip, is, iw, ij, 1)
            enddo ! iw
          enddo ! is
        enddo ! ip
      enddo ! ia
    enddo ! ij

    do is = 1, NS
      HEVs(is, -(JJ-2):0) = HEVs(is, -(JJ-2):0)/mass(is, -(JJ-2):0)
    enddo ! is
    HEV(-(JJ-2):0) = HEV(-(JJ-2):0)/mas

    ! calculate ex ante welfare of future generations
    do it = 1, TT
      HEV(it) = ((vv_coh(1, it)/vv_coh(1, 0))**(1d0/(1d0-1d0/gamma))-1d0)*100d0
    enddo

    ! headline
    write(22,'(/a,a)')'      A    KK    KC    LC   l_bar     r     w   inc     C     I',  &
      '    YY    YC  tauc  taup    PP   HEV    DIFF'

    ! current generations
    do ij = -(JJ-2), 0
      write(22,'(i3,80x,3f8.2,16x,f8.2)')ij,HEVs(:, ij),HEV(ij)
    enddo

    write(22,'(a)')'-------------------------------------------------------------------------------------------------'

    ! future generations
    do it = 1, TT
      write(22,'(i3,22f8.2,f10.5)')it,(/AA(it)/AA(0)-1d0, KK(it)/KK(0)-1d0, KC(it)/KC(0)-1d0, LC(it)/LC(0)-1d0, &
        HH(it)/sum(pop_w(:, it))/(HH(0)/sum(pop_w(:, 0)))-1d0, &
        (1d0+r(it))**0.2d0-(1d0+r(0))**0.2d0, w(it)/w(0)-1d0, inc_bar(it)/inc_bar(0)-1d0, CC(it)/CC(0)-1d0, &
        II(it)/II(0)-1d0, YY(it)/YY(0)-1d0, YC(it)/YC(0)-1d0, tauc(it)-tauc(0), taup(it)-taup(0), &
            PP(it)/YY(it)-PP(0)/YY(0)/)*100d0, HEV(it), DIFF(it)/YY(it)*100d0
    enddo

    if(lsra_on)write(22, '(/a,f12.6)')'EFFICIENCY GAIN: ', (Vstar**(1d0/(1d0-1d0/gamma))-1d0)*100d0

  end subroutine

end program
