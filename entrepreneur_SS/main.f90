program main

  ! modules
  use globals
  use omp_lib
  use clock
  use toolbox

  implicit none

  integer, parameter :: numthreads = 14
  integer :: ij
  real*8 :: shares_target(JJ, NS), shares_result(JJ, NS), share_target, share_result
  real*8 :: mu_val(5), sigma_val(5), rho_val(5)
  integer :: s3, h3, m3

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
  allocate(aplus(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(pplus(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(c(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(l(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(k(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(oplus(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(pencon(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(inctax(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(captax(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(m(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(VV(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))
  allocate(EV(0:1, 0:NA, 0:NP, NW, NE, NS, JJ))

  ! household preference parameters
  gamma  =  0.500d0
  ! invert gamma
  gamma = 1d0/gamma
  sigma  =   0.320d0
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
  delta = 0.06d0
  nu = 0.88d0
  l_bar = .47d0
  ! convert variables into per period values
  delta = 1d0 - (1d0-delta)**5d0

  ! demographic parameters
  n_p   = 0.007d0
  ! convert variables into per period values
  n_p = (1d0+n_p)**5-1d0

  ! tax and transfers
  taur   = 0.25d0
  tauy   = 0.15d0
  kappa  = 0.53d0
  sscc   = 2d0
  mu     = 1d0
  lambda = 0d0
  phi    = 0d0
  gy     = 0.19d0
  by     = 0.70d0
  ! convert variables into per period values
  by = by/5d0

  ! size of the asset grid
  a_l    = 0d0
  a_u    = 8192d0
  a_grow = 2.0d0

  ! size of the pension claim grid
  p_l  = 0d0
  p_u  = 2d0

  ! simulation parameters
  damp  = 0.60d0
  tol   = 1d-4
  itermax = 20

  ! compute gini
  gini_on = .true.

  ! set switches
  if (NO == 0) ent = .false.

  ! initialize target shares
  open(303, file='shares.dat')
  do ij = 1, JJ
    read(303,'(3f6.2)')shares_target(ij, :)
  enddo
  close(303)
  shares_target(:, :) = shares_target(:, :)

  share_target = 10.4035d0

  mu_val(:) = -(/0.02d0, 0.03d0, 0.05d0, 0.07d0, 0.1d0/)

  sigma_val(:) = (/0.037d0, 0.0365d0, 0.036d0, 0.0355d0, 0.035d0/)

  rho_val(:) = (/0.945d0, 0.94d0, 0.935d0, 0.930d0, 0.925d0/)

  open(307, file='results.out')

  suc = 0.625d0

  do m3 = 3, 3
    do s3 = 2, 2
      do h3 = 5, 5

        write(*,*)m3, s3, h3

        ! calculate initial equilibrium
        call get_SteadyState()

        share_result = sum(pop_e(:))/(sum(pop_w(:)+pop_e(:)))*100d0
        shares_result(:, 1) = (os_coh(1, 0, 1, :)+os_coh(1, 1, 1, :))*100d0
        shares_result(:, 2) = (os_coh(1, 0, 2, :)+os_coh(1, 1, 2, :))*100d0
        shares_result(:, 3) = (os_coh(1, 0, 3, :)+os_coh(1, 1, 3, :))*100d0


        write(307, '(3i3, 8f8.4)')m3, s3, h3, share_result, &
            sqrt(0d0*(share_target-share_result)**2d0 + sum((shares_target(:, 1)-shares_result(:, 1))**2d0) &
                                                      + sum((shares_target(:, 2)-shares_result(:, 2))**2d0) &
                                                      + sum((shares_target(:, 3)-shares_result(:, 3))**2d0)), &
            sqrt(1d0*(share_target-share_result)**2d0 + sum((shares_target(:, 1)-shares_result(:, 1))**2d0) &
                                                      + sum((shares_target(:, 2)-shares_result(:, 2))**2d0) &
                                                      + sum((shares_target(:, 3)-shares_result(:, 3))**2d0)), &
            sqrt(2d0*(share_target-share_result)**2d0 + sum((shares_target(:, 1)-shares_result(:, 1))**2d0) &
                                                      + sum((shares_target(:, 2)-shares_result(:, 2))**2d0) &
                                                      + sum((shares_target(:, 3)-shares_result(:, 3))**2d0)), &
            sqrt(4d0*(share_target-share_result)**2d0 + sum((shares_target(:, 1)-shares_result(:, 1))**2d0) &
                                                      + sum((shares_target(:, 2)-shares_result(:, 2))**2d0) &
                                                      + sum((shares_target(:, 3)-shares_result(:, 3))**2d0)), &
            sqrt(8d0*(share_target-share_result)**2d0 + sum((shares_target(:, 1)-shares_result(:, 1))**2d0) &
                                                      + sum((shares_target(:, 2)-shares_result(:, 2))**2d0) &
                                                      + sum((shares_target(:, 3)-shares_result(:, 3))**2d0)), &
            sqrt(16d0*(share_target-share_result)**2d0 + sum((shares_target(:, 1)-shares_result(:, 1))**2d0) &
                                                      + sum((shares_target(:, 2)-shares_result(:, 2))**2d0) &
                                                      + sum((shares_target(:, 3)-shares_result(:, 3))**2d0)), &
            sqrt(32d0*(share_target-share_result)**2d0 + sum((shares_target(:, 1)-shares_result(:, 1))**2d0) &
                                                      + sum((shares_target(:, 2)-shares_result(:, 2))**2d0) &
                                                      + sum((shares_target(:, 3)-shares_result(:, 3))**2d0))
          write(*,'(16f8.4)')shares_result(:, 1)
          write(*,'(16f8.4)')shares_target(:, 1)

          write(*,'(16f8.4)')shares_result(:, 2)
          write(*,'(16f8.4)')shares_target(:, 2)

          write(*,'(16f8.4)')shares_result(:, 3)
          write(*,'(16f8.4)')shares_target(:, 3)

      end do
    end do
  end do


!    call plot((/(dble(ij), ij=1,JJ)/), c_coh(0, :))
!    call plot((/(dble(ij), ij=1,JJ)/), c_coh(1, :))
!    call execplot
!
!    call plot((/(dble(ij), ij=1,JJ)/), a_coh(0, :))
!    call plot((/(dble(ij), ij=1,JJ)/), a_coh(1, :))
!    call execplot
!
!    call plot((/(dble(ij), ij=1,JJ)/), inc_coh(0, :))
!    call plot((/(dble(ij), ij=1,JJ)/), inc_coh(1, :))
!    call execplot

!    write(*,'(16f8.4)')shares_result(:, 1)
!    write(*,'(16f8.4)')shares_target(:, 1)
!
!    write(*,'(16f8.4)')shares_result(:, 2)
!    write(*,'(16f8.4)')shares_target(:, 2)
!
!    write(*,'(16f8.4)')shares_result(:, 3)
!    write(*,'(16f8.4)')shares_target(:, 3)
!
!    write(*,*)share_result

    call plot((/(dble(ij), ij=1,JJ)/), shares_target(:, 1), color='blue')
    call plot((/(dble(ij), ij=1,JJ)/), shares_target(:, 2), color='red')
    call plot((/(dble(ij), ij=1,JJ)/), shares_target(:, 3), color='green')

    call plot((/(dble(ij), ij=1,jj)/), (os_coh(1, 0, 1, :)+os_coh(1, 1, 1, :))*100d0, color='blue', linewidth=4d0)
    call plot((/(dble(ij), ij=1,jj)/), (os_coh(1, 0, 2, :)+os_coh(1, 1, 2, :))*100d0, color='red', linewidth=4d0)
    call plot((/(dble(ij), ij=1,jj)/), (os_coh(1, 0, 3, :)+os_coh(1, 1, 3, :))*100d0, color='green', linewidth=4d0)
    call execplot()

  ! close files
  close(307)
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
    integer :: iter, iamax(JJ)

    ! initialize remaining variables
    call initialize()

    ! start timer
    call tick(calc)

    ! iterate until value function converges
    do iter = 1, itermax

      !call tick(calc)

      ! get factor and other prices
      call get_prices()

      ! solve the household problem
      call solve_household()

      ! calculate the distribution of households over state space
      call get_distribution()

      ! aggregate individual decisions
      call aggregation()

      ! determine the government parameters
      call government()

      ! check the grid
      call check_grid(iamax)

      write(*,'(i4,6f8.2,i7,f14.8)')iter, (/5d0*KK, CC, II/)/YY*100d0, &
        ((1d0+r)**0.2d0-1d0)*100d0, w, sum(pop_e(:))/(sum(pop_w(:))+sum(pop_e(:)))*100d0, maxval(iamax), DIFF/YY*100d0

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
    integer :: ia, ip, iw, ie, is, ij
    real*8 :: adj

    write(*,'(/a/)')'INITIAL EQUILIBRIUM'
    write(*,'(a)')'ITER     K/Y     C/Y     I/Y       r       w     ent  iamax          DIFF'

    ! initialize asset grid
    a = grid_Cons_Grow(a_l, a_u, a_grow, NA)

    ! initialize pension claim grid
    p = grid_Cons_Equi(p_l, p_u, NP)

    ! get initial guess for savings decision
    do ij = 1, JJ
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ia = 0, NA
                aplus(0, ia, ip, iw, ie, is, ij) = max(a(ia)/2d0, a(1))
                aplus(1, ia, ip, iw, ie, is, ij) = max(a(ia)/2d0, a(1))
              enddo ! ia
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is
    enddo ! ij

    ! initial guess for investment decision
    k(:, :, :, :, :, :, :) = 1d-4

    ! initial guess for labor decision
    l(:, :, :, :, :, :, :) = 0.33d0

    ! distribution of skill classes
    dist_skill = (/0.1520d0, 0.5547d0, 0.2933d0/)

		! initialize survival probabilities for middle skilled
    open(319, file='sp.dat')
    do ij = 1, JJ+1
      read(319,'(f13.8)')psi(2, ij)
    enddo
    close(319)

    ! compute survival probabilities for high/low skilled
    psi(:, 1) = psi(2, 1)
    psi(:, JJ+1) = 0d0
    adj = 23d0
    do ij = 2, JJ
      psi(1, ij) = psi(2, ij) - exp(0.33d0*(dble(ij-1)-adj))
      psi(3, ij) = psi(2, ij) + exp(0.33d0*(dble(ij-1)-adj))
    enddo

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
    Gama(1) = 0.0d0*pop(1)
    Gama(2) = 0.0d0*pop(2)
    Gama(3) = 1.0d0*pop(3)
    Gama(4) = 1.2d0*pop(4)
    Gama(5) = 1.4d0*pop(5)
    Gama(6) = 1.6d0*pop(6)
    Gama(7) = 1.8d0*pop(7)
    Gama(8) = 1.8d0*pop(8)
    Gama(9) = 1.6d0*pop(9)
!    Gama(1:JR-1) = 1d0
    Gama(JR:JJ) = 0d0
    Gama = Gama/sum(Gama)

    ! initialize age earnings process
    open(302, file='eff.dat')
    do ij = 1, JJ
      read(302,'(3f12.8)')eff(ij, :)
    enddo
    close(302)
    eff(JE:, :) = 0d0

    ! initialize productivity shocks
    call discretize_AR(0.95666d0**5d0, 0.0d0, sigma5(0.95666d0, 0.02321d0), eta(:, 1), pi_eta(:, :, 1), dist_eta(:, 1))
    eta(:, 1) = exp(eta(:, 1))/sum(dist_eta(:, 1)*exp(eta(:, 1)))

    call discretize_AR(0.95687d0**5d0, 0.0d0, sigma5(0.95687d0, 0.02812d0), eta(:, 2), pi_eta(:, :, 2), dist_eta(:, 2))
    eta(:, 2) = exp(eta(:, 2))/sum(dist_eta(:, 2)*exp(eta(:, 2)))

    call discretize_AR(0.95828d0**5d0, 0.0d0, sigma5(0.95828d0, 0.03538d0), eta(:, 3), pi_eta(:, :, 3), dist_eta(:, 3))
    eta(:, 3) = exp(eta(:, 3))/sum(dist_eta(:, 3)*exp(eta(:, 3)))

    ! initialize entrepreneurial ability
    call discretize_AR(0.935d0**5d0, -0.410d0, sigma5(0.935d0, 0.036d0), theta(:, 1), pi_theta(:, :, 1), dist_theta(:, 1))
    theta(:, 1) = exp(theta(:, 1))

    call discretize_AR(0.940d0**5d0, -0.345d0, sigma5(0.940d0, 0.036d0), theta(:, 2), pi_theta(:, :, 2), dist_theta(:, 2))
    theta(:, 2) = exp(theta(:, 2))

    call discretize_AR(rho_val(h3)**5d0, mu_val(m3), sigma5(rho_val(h3), sigma_val(s3)), theta(:, 3), pi_theta(:, :, 3), dist_theta(:, 3))
    theta(:, 3) = exp(theta(:, 3))

    ! initial guesses for macro variables
    taup   = 0.10d0
    tauc   = 0.19d0
    inc_bar = 0.61d0
    bqs = (/0.02d0, 0.10d0, 0.20d0/)
    BQ = 0.50d0
    BB = 0d0
    KC = 4.93d0
    LC = 5.25d0

    ! open files
    open(21, file='output.out')

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
    r1 = 0.286d0*inc_bar*2d0
    r2 = 0.456d0*inc_bar*2d0
    r3 = 1.786d0*inc_bar*2d0

    b1 = (t2-t1)/(r2-r1)
    b2 = (t3-t2)/(r3-r2)

  end subroutine


  !##############################################################################
  ! SUBROUTINE solve_household
  !
  ! Determines the solution to the household optimization problem
  !##############################################################################
  subroutine solve_household()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: ia, ip, iw, ie, is, ij
    real*8 :: xy(2), fret, limit

    do ij = JJ, 1, -1

      !call tick(calc)
      !write(*,*)'Optimize for age: ', ij

      ! set up communication variables
      ij_com = ij

      if (ij >= JE) then

        ! set up communication variables
        ie_com = 1
        iw_com = 1
        io_com = 0

        !$omp parallel do copyin(io_com, iw_com, ie_com, ij_com) collapse(3) &
        !$omp             schedule(guided) private(xy, fret) num_threads(numthreads)
        do is = 1, NS
          do ip = 0, NP
            do ia = 0, NA

              ! set up communication variables
              is_com = is
              ip_com = ip
              ia_com = ia

              ! get initial guess for the individual choices
              xy(1) = max(aplus(0, ia, ip, 1, 1, is, ij), 1d-4)

              call fminsearch(xy(1), fret, a_l, a_u, valuefunc_r)

              ! copy decisions
              aplus(:, ia, ip,  :,  :, is, ij) = xy(1)
              pplus(:, ia, ip,  :,  :, is, ij) = p(ip)
              c(:, ia, ip,  :,  :, is, ij) = max(c_com, 1d-10)
              l(:, ia, ip,  :,  :, is, ij) = 0d0
              k(:, ia, ip,  :,  :, is, ij) = 0d0
              oplus(:, ia, ip,  :,  :, is, ij) = 0d0
              pencon(:, ia, ip,  :,  :, is, ij) = 0d0
              inctax(:, ia, ip,  :,  :, is, ij) = inctax_com
              captax(:, ia, ip,  :,  :, is, ij) = captax_com
              VV(:, ia, ip,  :,  :, is, ij) = -fret

            enddo ! ia
          enddo ! ip
        enddo ! is
        !$omp end parallel do

      elseif (ij >= JR) then

        if (ent) then

          ! set up communication variables
          iw_com = 1
          io_com = 1

          !$omp parallel do copyin(io_com, iw_com, ij_com) collapse(4) &
          !$omp             schedule(guided) private(xy, fret) num_threads(numthreads)
          do is = 1, NS
            do ie = 1, NE
              do ip = 0, NP
                do ia = 0, NA

                  ! set up communication variables
                  is_com = is
                  ie_com = ie
                  ip_com = ip
                  ia_com = ia

                  ! get initial guess for the individual choices
                  xy(1) = max(aplus(1, ia, ip, 1, ie, is, ij), 1d-4)
                  xy(2) = max(k(1, ia, ip, 1, ie, is, ij), 1d-4)

                  limit = max(1.5d0*a(ia), 1d-4)

                  call fminsearch(xy, fret, (/a_l, 0d0/), (/a_u, limit/), valuefunc_e)

                  ! copy decisions
                  aplus(1, ia, ip, :, ie, is, ij) = xy(1)
                  k(1, ia, ip, :, ie, is, ij) = k_com
                  pplus(1, ia, ip, :, ie, is, ij) = pplus_com
                  c(1, ia, ip, :, ie, is, ij) = max(c_com, 1d-10)
                  l(1, ia, ip, :, ie, is, ij) = l_com
                  oplus(1, ia, ip, :, ie, is, ij) = oplus_com
                  pencon(1, ia, ip, :, ie, is, ij) = pencon_com
                  inctax(1, ia, ip, :, ie, is, ij) = inctax_com
                  captax(1, ia, ip, :, ie, is, ij) = captax_com
                  VV(1, ia, ip, :, ie, is, ij) = -fret

                enddo ! ia
              enddo ! ip
            enddo ! ie
          enddo ! is
          !$omp end parallel do

        endif

        ! set up communication variables
        ie_com = 1
        iw_com = 1
        io_com = 0

        !$omp parallel do copyin(io_com, iw_com, ie_com, ij_com) collapse(3) &
        !$omp             schedule(guided) private(xy, fret) num_threads(numthreads)
        do is = 1, NS
          do ip = 0, NP
            do ia = 0, NA

              ! set up communication variables
              is_com = is
              ip_com = ip
              ia_com = ia

              ! get initial guess for the individual choices
              xy(1) = max(aplus(0, ia, ip, 1, 1, is, ij), 1d-4)

              call fminsearch(xy(1), fret, a_l, a_u, valuefunc_r)

              ! copy decisions
              aplus(0, ia, ip,  :,  :, is, ij) = xy(1)
              pplus(0, ia, ip,  :,  :, is, ij) = p(ip)
              c(0, ia, ip,  :,  :, is, ij) = max(c_com, 1d-10)
              l(0, ia, ip,  :,  :, is, ij) = 0d0
              k(0, ia, ip,  :,  :, is, ij) = 0d0
              oplus(0, ia, ip,  :,  :, is, ij) = 0d0
              pencon(0, ia, ip,  :,  :, is, ij) = 0d0
              inctax(0, ia, ip,  :,  :, is, ij) = inctax_com
              captax(0, ia, ip,  :,  :, is, ij) = captax_com
              VV(0, ia, ip,  :,  :, is, ij) = -fret

            enddo ! ia
          enddo ! ip
        enddo ! is
        !$omp end parallel do

      elseif (ij >= 2) then

        !$omp parallel copyin(ij_com) private(xy, fret, limit) num_threads(numthreads)

        if (ent) then

          ! set up communication variables
          io_com = 1

          !$omp do collapse(4) schedule(guided)
          do is = 1, NS
            do ie = 1, NE
              do iw = 1, NW
                do ip = 0, NP
                  do ia = 0, NA

                    ! set up communication variables
                    is_com = is
                    ie_com = ie
                    iw_com = iw
                    ip_com = ip
                    ia_com = ia

                    ! get initial guess for the individual choices
                    xy(1) = max(aplus(1, ia, ip, iw, ie, is, ij), 1d-4)
                    xy(2) = max(k(1, ia, ip, iw, ie, is, ij), 1d-4)

                    limit = max(1.5d0*a(ia), 1d-4)

                    call fminsearch(xy, fret, (/a_l, 0d0/), (/a_u, limit/), valuefunc_e)

                    ! copy decisions
                    aplus(1, ia, ip, iw, ie, is, ij) = xy(1)
                    k(1, ia, ip, iw, ie, is, ij) = k_com
                    pplus(1, ia, ip, iw, ie, is, ij) = pplus_com
                    c(1, ia, ip, iw, ie, is, ij) = max(c_com, 1d-10)
                    l(1, ia, ip, iw, ie, is, ij) = l_com
                    oplus(1, ia, ip, iw, ie, is, ij) = oplus_com
                    pencon(1, ia, ip, iw, ie, is, ij) = pencon_com
                    inctax(1, ia, ip, iw, ie, is, ij) = inctax_com
                    captax(1, ia, ip, iw, ie, is, ij) = captax_com
                    VV(1, ia, ip, iw, ie, is, ij) = -fret

                  enddo ! ia
                enddo ! ip
              enddo ! iw
            enddo ! ie
          enddo ! is
          !$omp end do nowait

        endif

        ! set up communication variables
        io_com = 0

        !$omp do collapse(4) schedule(guided)
        do is = 1, NS
          do ie = 1, NE
            do iw = 1, NW
              do ip = 0, NP
                do ia = 0, NA

                  ! set up communication variables
                  is_com = is
                  ie_com = ie
                  iw_com = iw
                  ip_com = ip
                  ia_com = ia

                  ! get initial guess for the individual choices
                  xy(1) = max(aplus(0, ia, ip, iw, ie, is, ij), 1d-4)
                  xy(2) = max(l(0, ia, ip, iw, ie, is, ij), 1d-4)

                  call fminsearch(xy(:2), fret, (/a_l, 0d0/), (/a_u, 1d0/), valuefunc_w)

                  ! copy decisions
                  aplus(0, ia, ip, iw, ie, is, ij) = xy(1)
                  k(0, ia, ip, iw, ie, is, ij) = k_com
                  pplus(0, ia, ip, iw, ie, is, ij) = pplus_com
                  c(0, ia, ip, iw, ie, is, ij) = max(c_com, 1d-10)
                  l(0, ia, ip, iw, ie, is, ij) = l_com
                  oplus(0, ia, ip, iw, ie, is, ij) = oplus_com
                  pencon(0, ia, ip, iw, ie, is, ij) = pencon_com
                  inctax(0, ia, ip, iw, ie, is, ij) = inctax_com
                  captax(0, ia, ip, iw, ie, is, ij) = captax_com
                  VV(0, ia, ip, iw, ie, is, ij) = -fret

                enddo ! ia
              enddo ! ip
            enddo ! iw
          enddo ! ie
        enddo ! is
        !$omp end do
        !$omp end parallel

      elseif (ij == 1) then

        ! set up communication variables
        ip_com = 0
        ia_com = 0
        io_com = 0

	      !$omp parallel do copyin(io_com, ia_com, ip_com, ij_com) collapse(3) &
	      !$omp             schedule(guided) private(xy, fret) num_threads(numthreads)
        do is = 1, NS
          do ie = 1, NE
            do iw = 1, NW

              ! set up communication variables
              is_com = is
              ie_com = ie
              iw_com = iw

              ! get initial guess for the individual choices
              xy(1) = max(aplus(0, 0, 0, iw, ie, is, ij), 1d-4)
              xy(2) = max(l(0, 0, 0, iw, ie, is, ij), 1d-4)

              call fminsearch(xy(:2), fret, (/a_l, 0d0/), (/a_u, 1d0/), valuefunc_w)

              ! copy decisions
              aplus(:, :, :, iw, ie, is, ij) = xy(1)
              pplus(:, :, :, iw, ie, is, ij) = pplus_com
              c(:, :, :, iw, ie, is, ij) = max(c_com, 1d-10)
              l(:, :, :, iw, ie, is, ij) = l_com
              k(:, :, :, iw, ie, is, ij) = 0d0
              oplus(:, :, :, iw, ie, is, ij) = oplus_com
              pencon(:, :, :, iw, ie, is, ij) = pencon_com
              inctax(:, :, :, iw, ie, is, ij) = inctax_com
              captax(:, :, :, iw, ie, is, ij) = captax_com
              VV(:, :, :, iw, ie, is, ij) = -fret

            enddo ! iw
          enddo ! ie
        enddo ! is
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
    integer :: ia, ip, iw, ie, is, iw_p, ie_p

    !$omp parallel do collapse(4) schedule(dynamic,1) private(iw_p, ie_p) num_threads(numthreads)
    do is = 1, NS
      do ie = 1, NE
        do iw = 1, NW
          do ip = 0, NP
            do ia = 0, NA

              EV(:, ia, ip, iw, ie, is, ij) = 0d0
              do ie_p = 1, NE
                do iw_p = 1, NW
                  EV(0, ia, ip, iw, ie, is, ij) = EV(0, ia, ip, iw, ie, is, ij) &
                    +pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*VV(0, ia, ip, iw_p, ie_p, is, ij)
                  EV(1, ia, ip, iw, ie, is, ij) = EV(1, ia, ip, iw, ie, is, ij) &
                    +pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)*VV(1, ia, ip, iw_p, ie_p, is, ij)
                enddo ! iw_p
              enddo ! ie_p

              EV(:, ia, ip, iw, ie, is, ij) &
                = ((1d0-gamma)*EV(:, ia, ip, iw, ie, is, ij))**(1d0/(1d0-gamma))

            enddo ! ia
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
  subroutine get_distribution()

    implicit none

    !##### OTHER VARIABLES ####################################################
    integer :: io, ia, ip, iw, ie, is, ij, &
               io_p, iw_p, ie_p, ial, iar, ipl, ipr
    real*8 :: varpsi, varphi

    ! set distribution to zero
    m(:, :, :, :, :, :, :) = 0d0

    ! get initial distribution at age 1
    do is = 1, NS
      do ie = 1, NE
        do iw = 1, NW
          m(0, 0, 0, iw, ie, is, 1) = dist_eta(iw, is)*dist_theta(ie, is)*dist_skill(is)
        enddo ! iw
      enddo ! ie
    enddo ! is

    ! successively compute distribution over ages
    do ij = 2, JJ

      ! iterate over yesterdays gridpoints

      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ia = 0, NA
                do io = 0, NO

                  ! interpolate yesterday's savings decision
                  call linint_Grow(aplus(io, ia, ip, iw, ie, is, ij-1), &
                           a_l, a_u, a_grow, NA, ial, iar, varphi)

                  ! interpolate today's pension claims
                  call linint_Equi(pplus(io, ia, ip, iw, ie, is, ij-1), &
                                     p_l, p_u, NP, ipl, ipr, varpsi)

                  ! this year's occupation
                  io_p = int(oplus(io, ia, ip, iw, ie, is, ij-1))

                  ! redistribute households
                  do ie_p = 1, NE
                    do iw_p = 1, NW
                      m(io_p, ial, ipl, iw_p, ie_p, is, ij) = &
                         m(io_p, ial, ipl, iw_p, ie_p, is, ij) &
                         +varphi*varpsi*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is)&
                         *psi(is, ij)*m(io, ia, ip, iw, ie, is, ij-1)/(1d0+n_p)
                      m(io_p, ial, ipr, iw_p, ie_p, is, ij) = &
                         m(io_p, ial, ipr, iw_p, ie_p, is, ij) &
                         +varphi*(1d0-varpsi)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                         *psi(is, ij)*m(io, ia, ip, iw, ie, is, ij-1)/(1d0+n_p)
                      m(io_p, iar, ipl, iw_p, ie_p, is, ij) = &
                         m(io_p, iar, ipl, iw_p, ie_p, is, ij) &
                         +(1d0-varphi)*varpsi*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                         *psi(is, ij)*m(io, ia, ip, iw, ie, is, ij-1)/(1d0+n_p)
                      m(io_p, iar, ipr, iw_p, ie_p, is, ij) = &
                         m(io_p, iar, ipr, iw_p, ie_p, is, ij) &
                         +(1d0-varphi)*(1d0-varpsi)*pi_eta(iw, iw_p, is)*pi_theta(ie, ie_p, is) &
                         *psi(is, ij)*m(io, ia, ip, iw, ie, is, ij-1)/(1d0+n_p)
                    enddo ! iw_p
                  enddo ! ie_p

                enddo ! ia
              enddo ! ip
            enddo ! iw
          enddo ! ie
        enddo ! is
      enddo ! io

      !m(:, :, :, :, :, :, ij) = m(:, :, :, :, :, :, ij)/sum(m(:, :, :, :, :, :, ij))*pop(ij)

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
    integer :: io, ia, ip, iw, ie, is, ij, ial, iar
    real*8 :: LC_old, varphi

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
      do is = 1, NS
        do ie = 1, NE
          do iw = 1, NW
            do ip = 0, NP
              do ia = 0, NA
                do io = 0, NO

                  call linint_Grow(aplus(io, ia, ip, iw, ie, is, ij), a_l, a_u, a_grow, NA, ial, iar, varphi)

                  AA = AA + (varphi*a(ial) + (1d0-varphi)*a(iar)) &
                            *m(io, ia, ip, iw, ie, is, ij)/(1d0+n_p)
                  CC = CC + c(io, ia, ip, iw, ie, is, ij) &
                            *m(io, ia, ip, iw, ie, is, ij)
                  if (io == 0 .and. ij < JR) then
                    LC = LC + eff(ij, is)*eta(iw, is)*l(io, ia, ip, iw, ie, is, ij) &
                              *m(io, ia, ip, iw, ie, is, ij)
                    HH = HH + l(io, ia, ip, iw, ie, is, ij) &
                              *m(io, ia, ip, iw, ie, is, ij)
                    PC = PC + min(w*eff(ij, is)*eta(iw, is)*l(io, ia, ip, iw, ie, is, ij), sscc*inc_bar) &
                              *m(io, ia, ip, iw, ie, is, ij)
                  else
                    KE = KE + k(io, ia, ip, iw, ie, is, ij)*m(io, ia, ip, iw, ie, is, ij)
                    YE = YE + theta(ie, is)*(k(io, ia, ip, iw, ie, is, ij)**alpha*(eff(ij, is)*l_bar)**(1d0-alpha))**nu &
                              *m(io, ia, ip, iw, ie, is, ij)
                    if (ij < JR) PC = PC + phi*min(profent(k(io, ia, ip, iw, ie, is, ij), ij, ia, is, ie), sscc*inc_bar) &
                              *m(io, ia, ip, iw, ie, is, ij)
                    PE = PE + profent(k(io, ia, ip, iw, ie, is, ij), ij, ia, is, ie) &
                              *m(io, ia, ip, iw, ie, is, ij)
                    if (ij >= JR) then
                      PRE = PRE + profent(k(io, ia, ip, iw, ie, is, ij), ij, ia, is, ie) &
                                *m(io, ia, ip, iw, ie, is, ij)
                    endif
                  endif
                  PP = PP + pen(ip, ij)*m(io, ia, ip, iw, ie, is, ij)
                  bqs(is) = bqs(is) + (1d0+r)*(varphi*a(ial) + (1d0-varphi)*a(iar))*(1d0-psi(is, ij+1)) &
                                *m(io, ia, ip, iw, ie, is, ij)/(1d0+n_p)
                  TAc = TAc + tauc*c(io, ia, ip, iw, ie, is, ij) &
                            *m(io, ia, ip, iw, ie, is, ij)
                  TAr = TAr + captax(io, ia, ip, iw, ie, is, ij)*m(io, ia, ip, iw, ie, is, ij)
                  TAw = TAw + inctax(io, ia, ip, iw, ie, is, ij)*m(io, ia, ip, iw, ie, is, ij)
                  if (io == 1 .and. ij < JR) then
                    pop_e(is) = pop_e(is) + m(io, ia, ip, iw, ie, is, ij)
                  elseif (io == 1 .and. ij >= JR) then
                    pop_re(is) = pop_re(is) + m(io, ia, ip, iw, ie, is, ij)
                  elseif (io == 0 .and. ij < JR) then
                    pop_w(is) = pop_w(is) + m(io, ia, ip, iw, ie, is, ij)
                  else
                    pop_r(is) = pop_r(is) + m(io, ia, ip, iw, ie, is, ij)
                  endif
                  vv_coh(ij) = vv_coh(ij) + VV(io, ia, ip, iw, ie, is, ij) &
                                  *m(io, ia, ip, iw, ie, is, ij)

                enddo ! io
              enddo ! ia
            enddo ! ip
          enddo ! iw
        enddo ! ie
      enddo ! is
    enddo ! ij

    ! damping and other quantities
    LC = damp*LC + (1d0-damp)*LC_old

    if (smopec) then
      KC = damp*(LC*((r/(1d0-tauy)+delta)/alpha)**(1d0/(alpha-1d0)))+(1d0-damp)*KC
      BF = AA - KE - KC - BB - BP
      NEX = (n_p-r)*BF
    else
      KC = damp*(AA - KE - BB - BP)+(1d0-damp)*KC
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
    real*8 :: expend

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
    integer :: io, ia, ip, iw, ie, is, ij, io_p, iamax(JJ)
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

                  c_coh(io, ij) = c_coh(io, ij) + c(io, ia, ip, iw, ie, is, ij) &
                                      *m(io, ia, ip, iw, ie, is, ij)
                  a_coh(io, ij) = a_coh(io, ij) + a(ia)*m(io, ia, ip, iw, ie, is, ij)
                  io_p = int(oplus(io, ia, ip, iw, ie, is, ij))
                  o_coh(io, io_p, ij) = o_coh(io, io_p, ij) + m(io, ia, ip, iw, ie, is, ij)
                  os_coh(io, io_p, is, ij) = os_coh(io, io_p, is, ij) &
                                   + m(io, ia, ip, iw, ie, is, ij)
                  inc_coh(io, ij) = inc_coh(io, ij) + pen(ip, ij)*m(io, ia, ip, iw, ie, is, ij)
                  if (io == 1) then
                    k_coh(ij) = k_coh(ij) + k(io, ia, ip, iw, ie, is, ij) &
                                    *m(io, ia, ip, iw, ie, is, ij)
                    inc_coh(io, ij) = inc_coh(io, ij) + profent(k(io, ia, ip, iw, ie, is, ij), ij, ia, is, ie) &
                                          *m(io, ia, ip, iw, ie, is, ij)
                  else
                    inc_coh(io, ij) = inc_coh(io, ij) + w*eff(ij, is)*eta(iw, is)*l(io, ia, ip, iw, ie, is, ij) &
                                          *m(io, ia, ip, iw, ie, is, ij)
                  endif
                  if (aplus(io, ia, ip, iw, ie, is, ij) <= 1d-10) then
                    flc_coh(ij) = flc_coh(ij) + m(io, ia, ip, iw, ie, is, ij)
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
      do is = 1, NS
        os_coh(0, :, is, ij) = os_coh(0, :, is, ij)/sum(m(:, :, :, :, :, is, ij))
        os_coh(1, :, is, ij) = os_coh(1, :, is, ij)/sum(m(:, :, :, :, :, is, ij))
      end do ! is
    enddo ! ij

    ! Output
    write(21,'(a, i3/)')'EQUILIBRIUM YEAR '
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

    write(21, '(a,a)')' IJ   CONSw   CONSe   ASSw     ASSe    INCw    INCe    INVe     ENT    ENTs1  ENTs2   ENTs3', &
        '    ENTn    WORn     FLC      VALUE  IAMAX'
    write(21,'(a)')'------------------------------------------------------------------------------------------------&
                    -------------------------------------'
    do ij = 1, JJ
      write(21, '(i3, 14f8.3, f11.3, i7)')ij, c_coh(0, ij), c_coh(1, ij), a_coh(0, ij), a_coh(1, ij), &
          inc_coh(0, ij), inc_coh(1, ij), k_coh(ij), sum(o_coh(1, :, ij)), sum(os_coh(1, :, 1, ij)), sum(os_coh(1, :, 2, ij)), &
          sum(os_coh(1, :, 3, ij)), o_coh(0, 1, ij), o_coh(1, 0, ij), flc_coh(ij), vv_coh(ij), iamax(ij)
      if (ij == JR-1) write(21,'(a)')'------------------------------------------------------------------------------&
                                      -------------------------------------------------------'
    enddo

    if (gini_on) then

      write(21,'(a)')'----------------------------------------------------------------------------------------------&
                      ---------------------------------------'
      write(21,'(a/)')' '

      call gini_wealth
      call gini_income

      write(21, '(a)')'WEALTH    GINI     1%     5%    10%    20%    40%    60%    FLC'
      write(21,'(4x, f10.3, 8f7.2)')gini_w, percentiles_w(:)*100d0, sum(pop(:)*flc_coh(:))/sum(pop(:))*100d0

      write(21, '(a)')'INCOME    GINI     1%     5%    10%    20%    40%    60%'
      write(21,'(4x, f10.3, 6f7.2/)')gini_i, percentiles_i(:)*100d0

    endif

    write(21,'(a)')'------------------------------------------------------------------------------------------------&
                    -------------------------------------'
    write(21,'(a/)')'-----------------------------------------------------------------------------------------------&
                     --------------------------------------'

  end subroutine

end program
