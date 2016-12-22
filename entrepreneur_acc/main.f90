program main

  ! modules
  use globals
  use omp_lib
  use clock
  use toolbox

  implicit none

  integer, parameter :: numthreads = 14

  ! allocate arrays
  if(allocated(aplus))deallocate(aplus)
  if(allocated(pplus))deallocate(pplus)
  if(allocated(c))deallocate(c)
  if(allocated(l))deallocate(l)
  if(allocated(k))deallocate(k)
  if(allocated(oplus))deallocate(oplus)
  if(allocated(pencon))deallocate(pencon)
  if(allocated(inctax))deallocate(inctax)
  if(allocated(captax))deallocate(captax)
  if(allocated(m))deallocate(m)
  if(allocated(VV))deallocate(VV)
  if(allocated(EV))deallocate(EV)
  allocate(aplus(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(pplus(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(c(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(l(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(k(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(oplus(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(pencon(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(inctax(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(captax(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(m(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(VV(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))
  allocate(EV(0:1, NS, NE, NW, 0:NP, 0:NA, JJ))

  ! household preference parameters
  gamma  =  0.500d0
  ! invert gamma
  gamma = 1d0/gamma
  sigma  =  0.318d0
  phi1   = -2.000d0
  phi2   =  0.900d0
  ! invert phi2
  phi2 = 1d0/phi2
  sigmaq =  1.500d0
  beta   =  0.998d0
  ! convert variables into per period values
  beta = beta**5d0

  ! production parameters
  alpha = 0.36d0
  delta = 0.06d0
  nu = 0.88d0
  l_bar = .41d0
  ! convert variables into per period values
  delta = 1d0 - (1d0-delta)**5d0

  ! demographic parameters
  n_p   = 0.0064d0
  ! convert variables into per period values
  n_p = (1d0+n_p)**5-1d0

  ! tax and transfers
  tauc   = 0.10d0
  taur   = 0.25d0
  tauy   = 0.15d0
  kappa  = 0.55d0
  sscc   = 2d0
  mu     = 1d0
  lambda = 0d0
  phi    = 0d0
  taup   = 0.10d0
  gy     = 0.19d0
  by     = 0.80d0
  ! convert variables into per period values
  by = by/5d0

  ! size of the asset grid
  a_l    = 0d0
  a_u    = 128d0
  a_grow = 1.2d0

  ! size of the pension claim grid
  p_l  = 0d0
  p_u  = 2d0

  ! simulation parameters
  damp  = 0.60d0
  tol   = 1d-6
  itermax = 200

  ! compute gini
  gini_on = .false.

  ! set switches
  if (NO == 0) ent = .false.

  ! calculate initial equilibrium
  call get_SteadyState()

  ! set reform parameters
  !pen_debt = .true.
  !smopec = .true.
  !phi = 1d0
  !lambda = 1d0
  !mu = 0d0
  !labor = .false.
  !ent = .false.

  ! close files
  close(21)

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

      !call tick(calc)

      ! get factor and other prices
      call get_prices

      ! solve the household problem
      call solve_household(1)

      ! calculate the distribution of households over state space
      call get_distribution

      ! aggregate individual decisions
      call aggregation

      ! determine the government parameters
      call government

      write(*,'(i4,6f8.2,f14.8)')iter, (/5d0*KK, CC, II/)/YY*100d0, &
        ((1d0+r)**0.2d0-1d0)*100d0, w, sum(pop_e(:))/(sum(pop_w(:))+sum(pop_e(:)))*100d0, DIFF/YY*100d0

      if(abs(DIFF/YY)*100d0 < tol)then

        call tock(calc)
        call output
        return
      endif

      !call tock(calc)

    enddo

    call tock(calc)
    call output

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
    integer :: is, ie, iw, ip, ia, ij
    real*8 :: adj

    write(*,'(/a/)')'INITIAL EQUILIBRIUM'
    write(*,'(a)')'ITER     K/Y     C/Y     I/Y       r       w     ent          DIFF'

    ! initialize asset grid
    a = grid_Cons_Grow(a_l, a_u, a_grow, NA)

    ! initialize pension claim grid
    p = grid_Cons_Equi(p_l, p_u, NP)

    ! get initial guess for savings decision
    do ij = 1, JJ
      do ia = 0, NA
        do ip = 0, NP
          do iw = 1, NW
            do ie = 1, NE
              do is = 1, NS
                aplus(0, is, ie, iw, ip, ia, ij) = max(a(ia)/2d0, a(1))
                aplus(1, is, ie, iw, ip, ia, ij) = max(a(ia)/2d0, a(1))
              enddo ! is
            enddo ! ie
          enddo ! iw
        enddo ! ip
      enddo ! ia
    enddo ! ij

    ! initial guess for investment decision
    k(:, :, :, :, :, :, :) = 1d-4

    ! initial guess for labor decision
    l(:, :, :, :, :, :, :) = 0.33d0

    ! distribution of skill classes
    dist_skill = (/0.263d0, 0.545d0, 0.192d0/)

		! initialize survival probabilities for middle skilled
    open(301, file='sp.dat')
    do ij = 1, JJ+1
      read(301,'(f13.8)')psi(2, ij)
    enddo
    close(301)

    ! compute survival probabilities for high/low skilled
    psi(:, 1) = psi(2, 1)
    psi(:, JJ+1) = 0d0
    adj = 22d0
    do ij = 2, JJ
      psi(1, ij) = psi(2, ij) - exp(0.33d0*(dble(ij-1)-adj))
      psi(3, ij) = psi(2, ij) + exp(0.33d0*(dble(ij-1)-adj))
    enddo

    !psi(1, :) = psi(2, :)
    !psi(3, :) = psi(2, :)

    ! set up population structure
    rpop(:, 1) = 1d0
    do ij = 2, JJ
      rpop(:, ij) = psi(:, ij)*rpop(:, ij-1)/(1d0+n_p)
    enddo

    pop = 0d0
    do is = 1, NS
      pop(:) = pop(:) + rpop(is, :)*dist_skill(is)
    enddo

    ! set distribution of bequests
    Gama(1) = 0d0*pop(1)
    Gama(2) = 0d0*pop(2)
    Gama(3) = 1d0*pop(3)
    Gama(4) = 2d0*pop(4)
    Gama(5) = 3d0*pop(5)
    Gama(6) = 4d0*pop(6)
    Gama(7) = 3d0*pop(7)
    Gama(8) = 2d0*pop(8)
    Gama(9) = 1d0*pop(9)
    !Gama(1:JR-1) = 1d0
    Gama(JR:JJ) = 0d0
    Gama = Gama/sum(Gama)

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

    ! initialize entrepreneurial ability
!    theta       = (/0.000d0, 0.290d0, 1.000d0, 1.710d0/)*1.880d0
!    dist_theta    = (/0.554d0, 0.283d0, 0.099d0, 0.064d0/)
!    pi_theta(1,:)   = (/0.780d0, 0.220d0, 0.000d0, 0.000d0/)
!    pi_theta(2,:)   = (/0.430d0, 0.420d0, 0.150d0, 0.000d0/)
!    pi_theta(3,:)   = (/0.000d0, 0.430d0, 0.420d0, 0.150d0/)
!    pi_theta(4,:)   = (/0.000d0, 0.000d0, 0.220d0, 0.780d0/)

    call discretize_AR(0.75666d0**5d0, 0.0d0, sigma5(0.75666d0, 0.3538d0), theta(3, :), pi_theta(3, :, :), dist_theta(3, :))
    theta(3, :) = exp(theta(3, :))/sum(dist_theta(3, :)*exp(theta(3, :)))
    call discretize_AR(0.85687d0**5d0, 0.0d0, sigma5(0.85687d0, 0.2538d0), theta(2, :), pi_theta(2, :, :), dist_theta(2, :))
    theta(2, :) = exp(theta(2, :))/sum(dist_theta(2, :)*exp(theta(2, :)))
    call discretize_AR(0.95828d0**5d0, 0.0d0, sigma5(0.95828d0, 0.1538d0), theta(1, :), pi_theta(1, :, :), dist_theta(1, :))
    theta(1, :) = exp(theta(1, :))/sum(dist_theta(1, :)*exp(theta(1, :)))

    ! initial guesses for macro variables
    inc_bar = 0.58d0
    BQ = 0.39d0
    KC = 5.09d0
    LC = 5.15d0

    ! open files
    open(21, file='output.out')
    open(22, file='summary.out')

  end subroutine


  !##############################################################################
  ! SUBROUTINE get_prices
  !
  ! Determines factor prices, indiviudal bequests and pensions at time it
  !##############################################################################
  subroutine get_prices()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: ij
    real*8 :: dueink

    ! calculate inverse goods price
    pinv = 1d0/(1d0+tauc)

    ! calculate factor prices
    if (.not. smopec) then
      r = (1d0-tauy)*(alpha*(KC/LC)**(alpha-1d0)-delta)
      w = (1d0-alpha)*(KC/LC)**alpha
    endif

    ! calculate individual bequests
    beq(1, :) = Gama*bqs(1)/rpop(1, :)/dist_skill(1)
    beq(2, :) = Gama*bqs(2)/rpop(2, :)/dist_skill(2)
    beq(3, :) = Gama*bqs(3)/rpop(3, :)/dist_skill(3)

    ! calculate individual pensions
    pen(:, :) = 0d0
    do ij = JR, JJ
      pen(:, ij) = p(:)*kappa*inc_bar
    enddo

    ! determine the income tax system
    dueink = inc_bar

    r1 = 0.286d0*dueink*2d0
    r2 = 0.456d0*dueink*2d0
    r3 = 1.786d0*dueink*2d0

    b1 = (t2-t1)/(r2-r1)
    b2 = (t3-t2)/(r3-r2)

  end subroutine


  !##############################################################################
  ! SUBROUTINE solve_household
  !
  ! Determines the solution to the household optimization problem
  !##############################################################################
  subroutine solve_household(ij_in)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: ij_in

    !##### OTHER VARIABLES ####################################################
    integer :: is, ie, iw, ip, ia, ij
    real*8 :: xy(2), fret, limit

    do ij = JJ, 1, -1

      !call tick(calc)
      !write(*,*)'Optimize for age: ', ij

      ! set up communication variables
      ij_com = ij

      if (ij >= JR) then

        ! set up communication variables
        iw_com = 1
        ie_com = 1
        io_com = 0

        !$omp parallel do copyin(ij_com, iw_com, ie_com, io_com) collapse(2) &
        !$omp             schedule(dynamic, 1) private(xy, fret) num_threads(numthreads)
        do ia = 0, NA
          do ip = 0, NP
            do is = 1, NS

              ! set up communication variables
              ia_com = ia
              ip_com = ip
              is_com = is

              ! get initial guess for the individual choices
              xy(1) = max(aplus(0, is, 1, 1, ip, ia, ij), 1d-4)

              call fminsearch(xy(1), fret, a_l, a_u, valuefunc_r)

              ! copy decisions
              aplus(:, is, :, :, ip, ia, ij) = xy(1)
              pplus(:, is, :, :, ip, ia, ij) = p(ip)
              c(:, is, :, :, ip, ia, ij) = max(c_com, 1d-10)
              l(:, is, :, :, ip, ia, ij) = 0d0
              k(:, is, :, :, ip, ia, ij) = 0d0
              oplus(:, is, :, :, ip, ia, ij) = 0d0
              pencon(:, is, :, :, ip, ia, ij) = 0d0
              inctax(:, is, :, :, ip, ia, ij) = inctax_com
              captax(:, is, :, :, ip, ia, ij) = captax_com
              VV(:, is, :, :, ip, ia, ij) = -fret

            enddo ! is
          enddo ! ip
        enddo ! ia
        !$omp end parallel do

      elseif (ij >= 2) then


        !$omp parallel copyin(ij_com) private(xy, fret, limit) num_threads(numthreads)

        if (ent) then

          ! set up communication variables
          io_com = 1

          !$omp do collapse(2) schedule(dynamic, 1)
          do ia = 0, NA
            do ip = 0, NP
              do iw = 1, NW
                do ie = 1, NE
                  do is = 1, NS

                    ! set up communication variables
                    ia_com = ia
                    ip_com = ip
                    is_com = is
                    iw_com = iw
                    ie_com = ie

                    ! get initial guess for the individual choices
                    xy(1) = max(aplus(1, is, ie, iw, ip, ia, ij), 1d-4)
                    xy(2) = max(k(1, is, ie, iw, ip, ia, ij), 1d-4)

                    limit = max(1.5d0*a(ia), 1d-4)

                    call fminsearch(xy, fret, (/a_l, 0d0/), (/a_u, limit/), valuefunc_e)

                    ! copy decisions
                    aplus(1, is, ie, iw, ip, ia, ij) = xy(1)
                    k(1, is, ie, iw, ip, ia, ij) = k_com
                    pplus(1, is, ie, iw, ip, ia, ij) = pplus_com
                    c(1, is, ie, iw, ip, ia, ij) = max(c_com, 1d-10)
                    l(1, is, ie, iw, ip, ia, ij) = l_com
                    oplus(1, is, ie, iw, ip, ia, ij) = oplus_com
                    pencon(1, is, ie, iw, ip, ia, ij) = pencon_com
                    inctax(1, is, ie, iw, ip, ia, ij) = inctax_com
                    captax(1, is, ie, iw, ip, ia, ij) = captax_com
                    VV(1, is, ie, iw, ip, ia, ij) = -fret

                  enddo ! is
                enddo ! ie
              enddo ! iw
            enddo ! ip
          enddo ! ia
          !$omp end do nowait

        endif

        ! set up communication variables
        io_com = 0

        !$omp do collapse(2) schedule(dynamic, 1)
        do ia = 0, NA
          do ip = 0, NP
            do iw = 1, NW
              do ie = 1, NE
                do is = 1, NS

                  ! set up communication variables
                  ia_com = ia
                  ip_com = ip
                  iw_com = iw
                  ie_com = ie
                  is_com = is

                  ! get initial guess for the individual choices
                  xy(1) = max(aplus(0, is, ie, iw, ip, ia, ij), 1d-4)
                  xy(2) = max(l(0, is, ie, iw, ip, ia, ij), 1d-4)

                  call fminsearch(xy(:2), fret, (/a_l, 0d0/), (/a_u, 1d0/), valuefunc_w)

                  ! copy decisions
                  aplus(0, is, ie, iw, ip, ia, ij) = xy(1)
                  k(0, is, ie, iw, ip, ia, ij) = k_com
                  pplus(0, is, ie, iw, ip, ia, ij) = pplus_com
                  c(0, is, ie, iw, ip, ia, ij) = max(c_com, 1d-10)
                  l(0, is, ie, iw, ip, ia, ij) = l_com
                  oplus(0, is, ie, iw, ip, ia, ij) = oplus_com
                  pencon(0, is, ie, iw, ip, ia, ij) = pencon_com
                  inctax(0, is, ie, iw, ip, ia, ij) = inctax_com
                  captax(0, is, ie, iw, ip, ia, ij) = captax_com
                  VV(0, is, ie, iw, ip, ia, ij) = -fret

                enddo ! is
              enddo ! ie
            enddo ! iw
          enddo ! ip
        enddo ! ia
        !$omp end do
        !$omp end parallel

      elseif (ij == 1) then

        ! set up communication variables
        ia_com = 0
        ip_com = 0
        io_com = 0

	      !$omp parallel do copyin(ij_com, ia_com, ip_com, io_com) collapse(3) &
        !$omp             schedule(dynamic, 1) private(xy, fret) num_threads(numthreads)
        do iw = 1, NW
          do ie = 1, NE
            do is = 1, NS

              ! set up communication variables
              iw_com = iw
              ie_com = ie
              is_com = is

              ! get initial guess for the individual choices
              xy(1) = max(aplus(0, is, ie, iw, 0, 0, ij), 1d-4)
              xy(2) = max(l(0, is, ie, iw, 0, 0, ij), 1d-4)

              call fminsearch(xy(:2), fret, (/a_l, 0d0/), (/a_u, 1d0/), valuefunc_w)

              ! copy decisions
              aplus(:, is, ie, iw, :, :, ij) = xy(1)
              pplus(:, is, ie, iw, :, :, ij) = pplus_com
              c(:, is, ie, iw, :, :, ij) = max(c_com, 1d-10)
              l(:, is, ie, iw, :, :, ij) = l_com
              k(:, is, ie, iw, :, :, ij) = 0d0
              oplus(:, is, ie, iw, :, :, ij) = oplus_com
              pencon(:, is, ie, iw, :, :, ij) = pencon_com
              inctax(:, is, ie, iw, :, :, ij) = inctax_com
              captax(:, is, ie, iw, :, :, ij) = captax_com
              VV(:, is, ie, iw, :, :, ij) = -fret

            enddo ! is
          enddo ! ie
        enddo ! iw
	      !$omp end parallel do

      endif

      ! interpolate individual value function
      call interpolate(ij)

      !write(*,*)'Done!'
      !call tock(calc)

    enddo

  end subroutine


  !##############################################################################
  ! SUBROUTINE interpolate
  !
  ! Calculates the expected valuefunction of cohort ij at time it
  !##############################################################################
  subroutine interpolate(ij)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    integer, intent(in) :: ij

    !##### OTHER VARIABLES ####################################################
    integer :: is, ie, iw, ip, ia, ie_p, iw_p

    !!$omp parallel do collapse(2) schedule(dynamic,1) private(iw_p, ie_p) num_threads(numthreads)
    !$acc parallel loop
    do ia = 0, NA
      do ip = 0, NP
        do iw = 1, NW
          do ie = 1, NE
            do is = 1, NS

                EV(:, is, ie, iw, ip, ia, ij) = 0d0
                do iw_p = 1, NW
                  do ie_p = 1, NE
                    EV(0, is, ie, iw, ip, ia, ij) = EV(0, is, ie, iw, ip, ia, ij) &
                      +pi_eta(is, iw, iw_p)*pi_theta(is, ie, ie_p)*VV(0, is, ie_p, iw_p, ip, ia, ij)
                    EV(1, is, ie, iw, ip, ia, ij) = EV(1, is, ie, iw, ip, ia, ij) &
                      +pi_eta(is, iw, iw_p)*pi_theta(is, ie, ie_p)*VV(1, is, ie_p, iw_p, ip, ia, ij)
                  enddo ! ie_p
                enddo ! iw_p

                EV(:, is, ie, iw, ip, ia, ij) &
                  = ((1d0-gamma)*EV(:, is, ie, iw, ip, ia, ij))**(1d0/(1d0-gamma))

            enddo ! is
          enddo ! ie
        enddo ! iw
      enddo ! ip
    enddo ! ia
    !$acc end parallel loop
    !!$omp end parallel do

  end subroutine


  !##############################################################################
  ! SUBROUTINE get_distribution
  !
  ! Determines the invariant distribution of households at time it
  !##############################################################################
  subroutine get_distribution()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: io, is, ie, iw, ip, ia, ij, &
           io_p, ie_p, iw_p, ipl, ipr, ial, iar
    real*8 :: varpsi, varphi

    ! set distribution to zero
    m(:, :, :, :, :, :, :) = 0d0

    ! get initial distribution at age 1
    do iw = 1, NW
      do ie = 1, NE
        do is = 1, NS
          m(0, is, ie, iw, 0, 0, 1) = dist_skill(is)*dist_theta(is, ie)*dist_eta(is, iw)
        enddo ! is
      enddo ! ie
    enddo ! iw

    !write(*,*) sum(m(:, :, :, :, :, :, :, 1))

    ! successively compute distribution over ages
    do ij = 2, JJ

      ! iterate over yesterdays gridpoints

      do ia = 0, NA
        do ip = 0, NP
          do iw = 1, NW
            do ie = 1, NE
              do is = 1, NS
                do io = 0, NO

                  ! interpolate yesterday's savings decision
                  call linint_Grow(aplus(io, is, ie, iw, ip, ia, ij-1), &
                           a_l, a_u, a_grow, NA, ial, iar, varphi)


                  ! interpolate today's pension claims
                  call linint_Equi(pplus(io, is, ie, iw, ip, ia, ij-1), &
                                     p_l, p_u, NP, ipl, ipr, varpsi)

                  ! this year's occupation
                  io_p = int(oplus(io, is, ie, iw, ip, ia, ij-1))

                  ! redistribute households
                  do iw_p = 1, NW
                    do ie_p = 1, NE
                      m(io_p, is, ie_p, iw_p, ipl, ial, ij) = &
                         m(io_p, is, ie_p, iw_p, ipl, ial, ij) &
                         +varphi*varpsi*pi_eta(is, iw, iw_p)*pi_theta(is, ie, ie_p)&
                         *psi(is, ij)*m(io, is, ie, iw, ip, ia, ij-1)
                      m(io_p, is, ie_p, iw_p, ipr, ial, ij) = &
                         m(io_p, is, ie_p, iw_p, ipr, ial, ij) &
                         +varphi*(1d0-varpsi)*pi_eta(is, iw, iw_p)*pi_theta(is, ie, ie_p) &
                         *psi(is, ij)*m(io, is, ie, iw, ip, ia, ij-1)
                      m(io_p, is, ie_p, iw_p, ipl, iar, ij) = &
                         m(io_p, is, ie_p, iw_p, ipl, iar, ij) &
                         +(1d0-varphi)*varpsi*pi_eta(is, iw, iw_p)*pi_theta(is, ie, ie_p) &
                         *psi(is, ij)*m(io, is, ie, iw, ip, ia, ij-1)
                      m(io_p, is, ie_p, iw_p, ipr, iar, ij) = &
                         m(io_p, is, ie_p, iw_p, ipr, iar, ij) &
                         +(1d0-varphi)*(1d0-varpsi)*pi_eta(is, iw, iw_p)*pi_theta(is, ie, ie_p) &
                         *psi(is, ij)*m(io, is, ie, iw, ip, ia, ij-1)
                    enddo ! ie_p
                  enddo ! iw_p

                enddo ! io
              enddo ! is
            enddo ! ie
          enddo ! iw
        enddo ! ip
      enddo ! ia

      m(:, :, :, :, :, :, ij) = m(:, :, :, :, :, :, ij)/sum(m(:, :, :, :, :, :, ij))

      !write(*,*) sum(m(:, :, :, :, :, :, ij))

    enddo ! ij

  end subroutine


  !##############################################################################
  ! SUBROUTINE aggregation
  !
  ! Calculates aggregate quantities at time it
  !##############################################################################
  subroutine aggregation()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: io, is, ie, iw, ip, ia, ij, ial, iar
    real*8 :: LC_old, varchi, varphi

    !write(*,*)'Calculate Aggregation:'
    !call tick(calc)

    LC_old = LC

    ! calculate aggregates
    AA  = 0d0
    CC  = 0d0
    LC  = 0d0
    HH  = 0d0
    KE  = 0d0
    YE  = 0d0
    PE  = 0d0
    PRE   = 0d0
    PC  = 0d0
    BQ  = 0d0
    PP  = 0d0
    TAc   = 0d0
    TAr   = 0d0
    TAw   = 0d0
    TAy   = 0d0
    bqs(:) = 0d0
    pop_w(:) = 0d0
    pop_r(:) = 0d0
    pop_re(:) = 0d0
    pop_e(:) = 0d0
    vv_coh(:) = 0d0

    do ij = 1, JJ
      do ia = 0, NA
        do ip = 0, NP
          do iw = 1, NW
            do ie = 1, NE
              do is = 1, NS
                do io = 0, NO

                  call linint_Grow(aplus(io, is, ie, iw, ip, ia, ij), a_l, a_u, a_grow, NA, ial, iar, varphi)

                  AA = AA + (varphi*a(ial) + (1d0-varphi)*a(iar)) &
                            *m(io, is, ie, iw, ip, ia, ij)*pop(ij)/(1d0+n_p)
                  CC = CC + c(io, is, ie, iw, ip, ia, ij) &
                            *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  if (io == 0 .and. ij < JR) then
                    LC = LC + eff(ij, is)*eta(is, iw)*l(io, is, ie, iw, ip, ia, ij) &
                              *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                    HH = HH + l(io, is, ie, iw, ip, ia, ij) &
                              *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                    PC = PC + min(w*eff(ij, is)*eta(is, iw)*l(io, is, ie, iw, ip, ia, ij), sscc*inc_bar) &
                              *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  else
                    KE = KE + k(io, is, ie, iw, ip, ia, ij)*m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                    YE = YE + theta(is, ie)*(k(io, is, ie, iw, ip, ia, ij)**alpha*(eff(ij, is)*l_bar)**(1d0-alpha))**nu &
                              *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                    if (ij < JR) PC = PC + phi*min(profent(k(io, is, ie, iw, ip, ia, ij), ij, ia, is, ie), sscc*inc_bar) &
                              *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                    PE = PE + profent(k(io, is, ie, iw, ip, ia, ij), ij, ia, is, ie) &
                              *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                    if (ij >= JR) then
                      PRE = PRE + profent(k(io, is, ie, iw, ip, ia, ij), ij, ia, is, ie) &
                                *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                    endif
                  endif
                  PP = PP + pen(ip, ij)*m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  bqs(is) = bqs(is) + (1d0+r)*(varphi*a(ial) + (1d0-varphi)*a(iar))*(1d0-psi(is, ij+1)) &
                                *m(io, is, ie, iw, ip, ia, ij)*pop(ij)/(1d0+n_p)
                  TAc = TAc + tauc*c(io, is, ie, iw, ip, ia, ij) &
                            *m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  TAr = TAr + captax(io, is, ie, iw, ip, ia, ij)*m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  TAw = TAw + inctax(io, is, ie, iw, ip, ia, ij)*m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  if (io == 1 .and. ij < JR) then
                    pop_e(is) = pop_e(is) + m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  elseif (io == 1 .and. ij > JR) then
                    pop_re(is) = pop_re(is) + m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  elseif (io == 0 .and. ij < JR) then
                    pop_w(is) = pop_w(is) + m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  else
                    pop_r(is) = pop_r(is) + m(io, is, ie, iw, ip, ia, ij)*pop(ij)
                  endif
                  vv_coh(ij) = vv_coh(ij) + VV(io, is, ie, iw, ip, ia, ij) &
                                  *m(io, is, ie, iw, ip, ia, ij)*pop(ij)

                enddo ! io
              enddo ! is
            enddo ! ie
          enddo ! iw
        enddo ! ip
      enddo ! ia
    enddo ! ij

    ! damping and other quantities
    LC = damp*LC + (1d0-damp)*LC_old

    if (smopec) then
      KC = damp*(LC*((r/(1d0-tauy)+delta)/alpha)**(1d0/(alpha-1d0)))+(1d0-damp)*KC
      BF = AA - KE - KC - BB - BP - BA
      NEX = (n_p-r)*BF
    else
      KC = damp*(AA - KE - BB - BP - BA)+(1d0-damp)*KC
      BF = 0d0
      NEX = 0d0
    endif

    BQ = sum(bqs(:))

    KK = KC + KE

    II = (1d0+n_p)*KK - (1d0-delta)*KK
    YC = KC**alpha * LC**(1d0-alpha)
    YY = YC + YE

    TAy = TAy + tauy*(YC - delta*KC - w*LC)

    inc_bar = (w*LC + PE)/(sum(pop_w(:)) + sum(pop_e(:)))

    ! get difference on goods market
    DIFF = YY-CC-II-GG-NEX

    !write(*,*)'Done!'
    !call tock(calc)

  end subroutine


  !##############################################################################
  ! SUBROUTINE government
  !
  ! Calculates government parameters at time it
  !##############################################################################
  subroutine government()

    implicit none

    !##### OTHER VARIABLES ####################################################
    real*8 :: expend, PP_bar, PC_bar

    ! set government quantities and pension payments
    GG = gy*YY
    BB = by*YY

    ! calculate government expenditure
    expend = GG + (1d0+r)*BB - (1d0+n_p)*BB

    ! get budget balancing tax rate
    tauc = (expend - (TAr + TAw + TAy))/CC

    ! get budget balancing contribution rate to pension system
    taup = PP/PC

  end subroutine


  !##############################################################################
  ! SUBROUTINE output
  !
  ! Writes output of year it to file 21
  !##############################################################################
  subroutine output()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: io, is, ie, iw, ip, ia, ij, io_p, iamax(JJ)
    real*8 :: c_coh(0:1, JJ), a_coh(0:1, JJ), k_coh(JJ)
    real*8 :: inc_coh(0:1, JJ), o_coh(0:1, 0:1, JJ), os_coh(0:1, 0:1, NS, JJ), flc_coh(JJ)
    real*8 :: life_exp(NS), punb(JJ, NS)

    life_exp = 0d0
    do is = 1, NS
      punb(1, is) = psi(is, 1)
      life_exp(is) = life_exp(is) + 2.5d0*dble(1)*punb(1,is)*(1d0-psi(is, 2))
      do ij = 2, JJ
        punb(ij, is) = punb(ij-1, is)*psi(is, ij)
        life_exp(is) = life_exp(is) + (2.5d0 + 5d0*dble(ij-1))*punb(ij, is)*(1d0-psi(is, ij+1))
      enddo ! ij
    enddo ! is

    c_coh(:, :) = 0d0
    a_coh(:, :) = 0d0
    k_coh(:) = 0d0
    inc_coh(:, :) = 0d0
    o_coh(:, :, :) = 0d0
    os_coh(:, :, :, :) = 0d0
    flc_coh(:) = 0d0

    do ij = 1, JJ
      do ia = 0, NA
        do ip = 0, NP
          do iw = 1, NW
            do ie = 1, NE
              do is = 1, NS
                do io = 0, NO

                  c_coh(io, ij) = c_coh(io, ij) + c(io, is, ie, iw, ip, ia, ij) &
                                      *m(io, is, ie, iw, ip, ia, ij)
                  a_coh(io, ij) = a_coh(io, ij) + a(ia)*m(io, is, ie, iw, ip, ia, ij)
                  io_p = int(oplus(io, is, ie, iw, ip, ia, ij))
                  o_coh(io, io_p, ij) = o_coh(io, io_p, ij) + m(io, is, ie, iw, ip, ia, ij)
                  os_coh(io, io_p, is, ij) = os_coh(io, io_p, is, ij) &
                                   + m(io, is, ie, iw, ip, ia, ij)
                  if (io == 1) then
                    k_coh(ij) = k_coh(ij) + k(io, is, ie, iw, ip, ia, ij) &
                                    *m(io, is, ie, iw, ip, ia, ij)
                    inc_coh(io, ij) = inc_coh(io, ij) + profent(k(io, is, ie, iw, ip, ia, ij), ij, ia, is, ie) &
                                          *m(io, is, ie, iw, ip, ia, ij)
                  else
                    inc_coh(io, ij) = inc_coh(io, ij) + w*eff(ij, is)*eta(is, iw)*l(io, is, ie, iw, ip, ia, ij) &
                                          *m(io, is, ie, iw, ip, ia, ij)
                  endif
                  if (aplus(io, is, ie, iw, ip, ia, ij) <= 1d-10) then
                    flc_coh(ij) = flc_coh(ij) + m(io, is, ie, iw, ip, ia, ij)
                  endif

                enddo ! io
              enddo ! is
            enddo ! ie
          enddo ! iw
        enddo ! ip
      enddo ! ia
    enddo ! ij

    do ij = 1, JJ
      c_coh(0, ij) = c_coh(0, ij)/sum(m(0, :, :, :, :, :, ij))
      c_coh(1, ij) = c_coh(1, ij)/sum(m(1, :, :, :, :, :, ij))
      a_coh(0, ij) = a_coh(0, ij)/sum(m(0, :, :, :, :, :, ij))
      a_coh(1, ij) = a_coh(1, ij)/sum(m(1, :, :, :, :, :, ij))
      inc_coh(0, ij) = inc_coh(0, ij)/sum(m(0, :, :, :, :, :, ij))
      inc_coh(1, ij) = inc_coh(1, ij)/sum(m(1, :, :, :, :, :, ij))
      k_coh(ij) = k_coh(ij)/sum(m(1, :, :, :, :, :, ij))
    enddo ! ij

    ! Output
    write(21,'(a/)')'EQUILIBRIUM YEAR '
    write(21,'(a)')'CAPITAL       KK      KC      KE      AA       r    p.a.'
    write(21,'(8x,6f8.2)')KK, KC, KE, AA, r, ((1d0+r)**(1d0/5d0)-1d0)*100d0
    write(21,'(a,4f8.2/)')'(in %)  ',(/KK, KC, KE, AA/)/YY*500d0

    write(21,'(a)')'LABOR         LC       w     inc   l_bar'
    write(21,'(8x,4f8.2/)')LC, w, inc_bar, HH/sum(pop_w(:))

    write(21,'(a)')'GOODS         YY      YC      YE      CC      II      NX    DIFF'
    write(21,'(8x,6f8.2,f8.3)')YY, YC, YE, CC, II, NEX, diff
    write(21,'(a,6f8.2,f8.3/)')'(in %)  ',(/YY, YC, YE, CC, II, NEX, diff/)/YY*100d0

    write(21,'(a)')'GOV         TAUC    TAUR    TAUW    TAUY   TOTAL      GG      BB      BF'
    write(21,'(8x,8f8.2)')TAc, TAr, TAw, TAy, TAc+TAr+TAw+TAy, GG, BB, BF
    write(21,'(a,8f8.2)')'(in %)  ',(/TAc, TAr, TAw, TAy, TAc+TAr+TAw+TAy, &
               GG, BB*5d0, BF*5d0/)/YY*100d0
    write(21,'(a,4f8.2/)')'(rate)  ',(/tauc, taur, 0d0, tauy/)*100d0

    write(21,'(a)')'PENS        TAUP     PEN      PP      PC      BQ'
    write(21,'(8x,5f8.2)')taup*PC, PP/sum(pop(JR:)), PP, PC, BQ
    write(21,'(a,5f8.2/)')'(in %)  ',(/taup, kappa, PP/YY, PC/YY, BQ/YY/)*100d0

    write(21,'(a)')'INCOME     TOTAL     WOR     ENT ENT/WOR     y_w     y_e y_e/y_w'
    write(21, '(8x,7f8.2)')w*LC + PE + PRE, w*LC, PE + PRE, (PE+PRE)/(w*LC), w*LC/sum(pop_w(:)), &
                           (PE + PRE)/(sum(pop_e(:) + pop_re(:))), (PE + PRE)/sum(pop_e(:) + pop_re(:))/w*LC/sum(pop_w(:))
    write(21, '(a,4f8.2/)')'(in %)  ',(/w*LC + PE + PRE, w*LC, PE + PRE/)/(w*LC + PE + PRE)*100d0, (PE+PRE)/(w*LC)*100d0

    write(21,'(a)')'POP        TOTAL     65-     65+ 65+/65-'
    write(21,'(8x,3f8.2)')sum(pop(:)), sum(pop(:JR-1)), sum(pop(JR:))
    write(21,'(a,4f8.2/)')'(in %)  ',(/sum(pop(:)), sum(pop(:JR-1)), sum(pop(JR:))/)/sum(pop(:))*100d0, &
                sum(pop(JR:))/sum(pop(:JR-1))*100d0

    write(21,'(a)')'           WORFO     WOR     ENT    WOR1    ENT1    WOR2    ENT2    WOR3    ENT3'
    write(21,'(8x,9f8.2)')sum(pop_w(:))+sum(pop_e(:)), sum(pop_w(:)), sum(pop_e(:)), &
                pop_w(1), pop_e(1), pop_w(2), pop_e(2), pop_w(3), pop_e(3)
    write(21, '(a, 9f8.2/)')'(in %)  ',(/sum(pop_w(:)+pop_e(:)), sum(pop_w(:)), sum(pop_e(:)), &
                pop_w(1), pop_e(1), pop_w(2), pop_e(2), pop_w(3), pop_e(3)/)/(sum(pop_w(:)+pop_e(:)))*100d0

    write(21,'(a)')'LIFE       j_bar  j_bar1  j_bar2  j_bar3'
    write(21,'(8x,4f8.2/)')sum(life_exp*dist_skill), life_exp(:)

    ! check for the maximium grid point used

    call check_grid(iamax)

    write(21, '(a,a)')' IJ   CONSw   CONSe   ASSw     ASSe    INCw    INCe    INVe     ENT', &
        '   ENTs1   ENTs2   ENTs3    ENTn    WORn     FLC      VALUE  IAMAX'
    write(21,'(a, a)')'-------------------------------------------------------------------', &
                      '-------------------------------------------------------------------'
    do ij = 1, JJ
      write(21, '(i3, 14f8.3, f11.3, i7)')ij, c_coh(0, ij), c_coh(1, ij), a_coh(0, ij), a_coh(1, ij), &
          inc_coh(0, ij), inc_coh(1, ij), k_coh(ij), sum(o_coh(1, :, ij)), sum(os_coh(1, :, 1, ij))/sum(os_coh(:, :, 1, ij)), &
          sum(os_coh(1, :, 2, ij))/sum(os_coh(:, :, 2, ij)), sum(os_coh(1, :, 3, ij))/sum(os_coh(:, :, 3, ij)), o_coh(0, 1, ij), &
          o_coh(1, 0, ij), flc_coh(ij), vv_coh(ij), iamax(ij)
      if (ij == JR-1) write(21,'(a, a)')'----------------------------------------------------------', &
                      '----------------------------------------------------------------------------'
    enddo

    if (gini_on) then

      write(21,'(a, a)')'-------------------------------------------------------------------', &
                        '-------------------------------------------------------------------'
      write(21,'(a/)')' '

      call gini_wealth
      call gini_income

      write(21, '(a)')'WEALTH  GINI   1%   5%  10%  20%  40%  60%  FLC'
      write(21,'(4x, f10.3, 8f7.2)')gini_w, percentiles_w(:)*100d0, sum(pop(:)*flc_coh(:))/sum(pop(:))*100d0

      write(21, '(a)')'INCOME  GINI   1%   5%  10%  20%  40%  60%'
      write(21,'(4x, f10.3, 6f7.2/)')gini_i, percentiles_i(:)*100d0

    endif

      write(21,'(a, a)')'--------------------------------------------------------------------', &
                        '--------------------------------------------------------------------'
      write(21,'(a, a/)')'-------------------------------------------------------------------', &
                        '--------------------------------------------------------------------'

  end subroutine

end program
