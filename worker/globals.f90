module globals

  use toolbox
  use sorting

  implicit none

  ! number of years the household lives
  integer, parameter :: JJ = 16
  
  ! number of years the household can be a worker (+1)
  integer, parameter :: JR = 10
  
  ! number of transition periods
  integer, parameter :: TT = 48

  ! number of permanent skill classes
  integer, parameter :: NS = 3
  
  ! number of transitory shock process values (worker)
  integer, parameter :: NW = 5
  
  ! number of points on the asset grid (-1)
  integer, parameter :: NA = 11
  
  ! number of points on the pension claim grid (-1)
  integer, parameter :: NP = 4
  
  ! household parameters  
  real*8 :: gamma, sigma, beta

  ! production parameters
  real*8 :: alpha, delta

  ! numerical parameters
  real*8 :: a_l, a_u, a_grow
  real*8 :: p_l, p_u, p_grow
  real*8 :: damp, tol
  integer :: itermax
  
  ! counter variables
  integer :: iter

  ! macroeconomic variables
  real*8 :: gy, by
  real*8 :: r(0:TT), w(0:TT), inc_bar(0:TT), px(JJ, 0:TT)
  real*8 :: KK(0:TT), KC(0:TT), AA(0:TT), LC(0:TT), HH(0:TT), BB(0:TT), BF(0:TT), BQ(0:TT)
  real*8 :: YY(0:TT), YC(0:TT), CC(0:TT), II(0:TT), GG(0:TT), NX(0:TT)
  real*8 :: TAc(0:TT), TAr(0:TT), TAw(0:TT), TAy(0:TT)
  logical :: smopec = .false.
  
  ! government variables
  real*8 :: tauc(0:TT), taur(0:TT), tauy(0:TT), tauw(0:TT)
  real*8 :: kappa(0:TT), lambda(0:TT), mu(0:TT), pen(JJ, 0:NP, 0:TT), taup(0:TT), PP(0:TT), PC(0:TT), BP(0:TT) = 0d0
  logical :: pen_debt = .false.
  
  ! LSRA variables
  real*8 :: BA(0:TT) = 0d0, SV(0:TT) = 0d0, lsra_comp, lsra_all, Vstar
  logical :: lsra_on

  ! cohort aggregate variables
  real*8 :: pop_w(NS, 0:TT), pop_r(NS, 0:TT)
  real*8 :: c_coh(JJ, 0:TT), a_coh(JJ, 0:TT), x_coh(JJ, 0:TT), bx_coh(JJ, 0:TT), inc_coh(JJ, 0:TT)
  real*8 :: flc_coh(JJ, 0:TT)
  real*8 :: vv_coh(JJ, 0:TT) = 0d0

  ! the shock process
  real*8 :: dist_skill(NS)
  real*8 :: eta(NS, NW), dist_eta(NS, NW), pi_eta(NS, NW, NW)

  ! demographic and other model parameters
  real*8 :: eff(JJ, NS), rpop(JJ, NS, 0:TT), pop(JJ, 0:TT), psi(JJ+1, NS), psi_avg(JJ+1), beq(JJ, 0:TT), Gama(JJ), n_p
  
  ! individual variables
  real*8 :: a(0:NA), p(0:NP)
  real*8 :: aplus(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: xplus(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: pplus(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: c(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: l(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: mx(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: pencon(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: inctax(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: captax(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: m(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: VV(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  real*8 :: v(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT) = 0d0
  
  ! numerical variables
  real*8 :: EV(0:NA, 0:NA, 0:NP, NS, NW, JJ, 0:TT)
  integer :: ia_com, ix_com, ip_com, is_com, iw_com, ij_com, it_com
  real*8 :: c_com, l_com, mx_com, xplus_com, pplus_com, pencon_com, inctax_com, captax_com, DIFF(0:TT)

  ! statistical variables
  logical :: gini_on = .false.
  real*8 :: gini(0:TT, 0:1), percentiles(6, 0:TT, 0:1)
  
  !$omp threadprivate(ia_com, ix_com, ip_com, is_com, iw_com, ij_com, it_com)  
  !$omp threadprivate(c_com, l_com, mx_com, xplus_com, pplus_com, pencon_com, inctax_com, captax_com)
  
contains


  !##############################################################################
  ! FUNCTION valuefunc_w
  !
  ! Determines the value function of a worker
  !##############################################################################
  function valuefunc_w(x)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################
    real*8, intent(in) :: x(:)
    real*8 :: valuefunc_w
    
    !##### OTHER VARIABLES ####################################################
    real*8 :: a_plus, wage, v_ind, valuefunc_help, varphi, varchi, varpsi
    real*8 :: p_hat, temp1, temp2
    integer :: ij, iij, itj, itp, ial, iar, ixl, ixr, ipl, ipr

    ! tomorrow's assets
    a_plus = x(1)
    
    ! today's labor
    l_com  = x(2)

    ! investment in annuitized assets
    mx_com = x(3)
 
    ! calculate annuities
    p_hat = 0d0
    
    if (ij_com >= JR) then
      temp1 = 0d0
      do ij = ij_com, JJ
        temp2 = 1d0
        do iij = ij_com+1, ij
          itj = year(it_com, ij_com, iij)
          temp2 = temp2*(1d0+r(itj))/px(iij, itj)
        enddo
        temp1 = temp1 + 1d0/temp2
      enddo
      p_hat = a(ix_com)*(1d0+r(it_com))/px(ij_com, it_com)/temp1
    endif      

    ! calculate tomorrow's annuitized asset stock
    if (ij_com <= JR-1) then  
      xplus_com = a(ix_com)*(1d0+r(it_com))/px(ij_com, it_com) - p_hat + mx_com
    else
      xplus_com = a(ix_com)*(1d0+r(it_com))/px(ij_com, it_com) - p_hat
    endif

    ! get tomorrow's year
    itp = year(it_com, ij_com, ij_com+1)  
    
    ! get lsra transfer payment
    v_ind = v(ia_com, ix_com, ip_com, is_com, iw_com, ij_com, it_com)
     
    ! calculate the wage rate and next periods pension claims
    if (ij_com < JR) then
      wage = w(it_com)*eff(ij_com, is_com)*eta(is_com, iw_com)
      pplus_com = (p(ip_com)*dble(ij_com-1) + mu(it_com)*(lambda(it_com) &
             + (1d0-lambda(it_com))*min(wage*l_com/inc_bar(it_com), 2d0)))/dble(ij_com)
    else
      wage = 0d0
      pplus_com = p(ip_com)
    endif
       
    ! calculate contribution to pension system
    pencon_com = taup(it_com)*min(wage*l_com, 2d0*inc_bar(it_com))
    
    ! calculate income tax
    inctax_com = tauw(it_com)*wage*l_com
    
    ! calculate capital gains tax
    captax_com = taur(it_com)*r(it_com)*a(ia_com)
    
    ! calculate consumption
    c_com = ((1d0+r(it_com))*a(ia_com) + wage*l_com + beq(ij_com, it_com) + pen(ij_com, ip_com, it_com) + p_hat + v_ind &
         - pencon_com - inctax_com - captax_com - mx_com - a_plus)/(1d0+tauc(it_com))   
    
    ! calculate tomorrow's part of the value function and occupational decision
    valuefunc_w = 0d0
    valuefunc_help = 0d0

    if (ij_com < JJ) then
        
      ! interpolate next period's value function as a worker/retiree
      call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)
      call linint_Grow(xplus_com, a_l, a_u, a_grow, NA, ixl, ixr, varchi)
      if (p_grow > 0d0) call linint_Grow(pplus_com, p_l, p_u, p_grow, NP, ipl, ipr, varpsi)
      if (p_grow <= 0d0) call linint_Equi(pplus_com, p_l, p_u, NP, ipl, ipr, varpsi)

      valuefunc_w = (varphi*varchi*varpsi*EV(ial, ixl, ipl, is_com, iw_com, ij_com+1, itp) &
               + varphi*varchi*(1d0-varpsi)*EV(ial, ixl, ipr, is_com, iw_com, ij_com+1, itp) &
               + varphi*(1d0-varchi)*varpsi*EV(ial, ixr, ipl, is_com, iw_com, ij_com+1, itp) &
               + varphi*(1d0-varchi)*(1d0-varpsi)*EV(ial, ixr, ipr, is_com, iw_com, ij_com+1, itp) &
               + (1d0-varphi)*varchi*varpsi*EV(iar, ixl, ipl, is_com, iw_com, ij_com+1, itp) &
               + (1d0-varphi)*varchi*(1d0-varpsi)*EV(iar, ixl, ipr, is_com, iw_com, ij_com+1, itp) &
               + (1d0-varphi)*(1d0-varchi)*varpsi*EV(iar, ixr, ipl, is_com, iw_com, ij_com+1, itp) &
               + (1d0-varphi)*(1d0-varchi)*(1d0-varpsi)*EV(iar, ixr, ipr, is_com, iw_com, ij_com+1, itp)) &
              **(1d0-1d0/gamma)/(1d0-1d0/gamma)
      
    endif
    
    ! add today's part and discount
    valuefunc_w = -(util(c_com, l_com) + beta*psi(ij_com+1, is_com)*valuefunc_w)
    
  end function

  
  !##############################################################################
  ! FUNCTION util
  !
  ! Determines instantaneous utility
  !##############################################################################  
  function util(cons, lab)
  
    implicit none
    
    !##### INPUT/OUTPUT VARIABLES #############################################    
    real*8, intent(in) :: cons, lab
    real*8 :: util
    
    ! determine utility
    if (cons <= 0d0 .or. lab < 0d0 .or. lab >= 1d0) then
      if(cons <= 0d0)util=-1d18*(1d0+abs(cons))
      if(lab < 0d0)util=-1d18*(1d0+abs(lab))
      if(lab >= 1d0)util=-1d18*(1d0+abs(lab))
    else
      util = (cons**sigma*(1d0-lab)**(1d0-sigma))**(1d0-1d0/gamma)/&
          (1d0-1d0/gamma)
    endif

  end function
  
  
  !##############################################################################
  ! FUNCTION margu
  !
  ! Determines marginal utility
  !##############################################################################  
  function margu(cons, lab, it)

    !##### INPUT/OUTPUT VARIABLES #############################################   
    real*8, intent(in) :: cons, lab
    integer, intent(in) :: it
    real*8 :: margu

    ! determine marginal utility
    margu = sigma*(cons**sigma*(1d0-lab)**(1d0-sigma))**(1d0-1d0/gamma)/((1d0+tauc(it))*cons)

  end function
  
  
  !##############################################################################
  ! FUNCTION year
  !
  ! Calculates year at which age ij agent is ijj
  !##############################################################################
  function year(it, ij, ijj)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################    
    integer, intent(in) :: it, ij, ijj
    integer :: year
    
    ! calculate year
    year = it + ijj - ij

    if(it == 0 .or. year <= 0)year = 0
    if(it == TT .or. year >= TT)year = TT

  end function

  
  !##############################################################################
  ! FUNCTION check_grid
  !
  ! Checks for the maximum gridpoint used
  !############################################################################## 
  subroutine check_grid(iamax, ixmax, it)
  
    implicit none

    !##### INPUT/OUTPUT VARIABLES ############################################# 
    integer, intent(in) :: it
    integer :: iamax(JJ), ixmax(JJ)

    !##### OTHER VARIABLES ####################################################    
    integer :: ij, ia, ix, ip, is, iw, ie

    iamax = 0
    do ij = 1, JJ

      ! check for the maximum asset grid point used at a certain age
      do iw = 1, NW
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NA
              do ia = 0, NA
                if (m(ia, ix, ip, is, iw, ij, it) > 0d0) iamax(ij)=ia
                if (m(ia, ix, ip, is, iw, ij, it) > 0d0) ixmax(ij)=ix
              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
      enddo ! iw
      
    enddo

  end subroutine

  
  !##############################################################################
  ! FUNCTION sigma5
  !
  ! Converts sigma of an one-year process into sigma of an five-year process
  !##############################################################################   
  function sigma5(rho, sigma)
  
    implicit none
    
    !##### INPUT/OUTPUT VARIABLES ############################################# 
    real*8, intent(in) :: rho, sigma
    real*8 :: sigma5

    !##### OTHER VARIABLES ####################################################    
    real*8 :: Gamma5
    integer :: i, j
    
    Gamma5 = 0d0
    do i = 0, 4
      do j = 0, 4
        Gamma5 = Gamma5 + rho**abs(j-i)
      enddo
    enddo
    Gamma5 = Gamma5/25d0
    
    sigma5 = Gamma5*sigma*(1d0-rho**10d0)/(1d0-rho**2d0)
    
  end function
  
    
  !##############################################################################
  ! FUNCTION gini_income
  !
  ! Calculates GINI coefficient of income (ikind = 0 for wealth, = 1 for income)
  !##############################################################################
  subroutine compute_gini(it, ikind)
  
    implicit none
    
    !##### INPUT/OUTPUT VARIABLES #############################################    
    integer, intent(in) :: it, ikind
    
    !##### OTHER VARIABLES ####################################################    
    integer :: ij, ia, ix, ip, iw, is, ic, ICMAX
    real*8, allocatable :: xs(:), ys(:), xcum(:), ycum(:)
    
    if(allocated(xs))deallocate(xs)
    if(allocated(ys))deallocate(ys)
    if(allocated(xcum))deallocate(xcum)
    if(allocated(ycum))deallocate(ycum)
  
    allocate(xs((JR-1)*(NA+1)*(NA+1)*(NP+1)*NS*NW))
    allocate(ys((JR-1)*(NA+1)*(NA+1)*(NP+1)*NS*NW))
    allocate(xcum(0:(JR-1)*(NA+1)*(NA+1)*(NP+1)*NS*NW))
    allocate(ycum(0:(JR-1)*(NA+1)*(NA+1)*(NP+1)*NS*NW))
    
    ! check input
    if (ikind > 1 .or. ikind < 0) call error('gini', 'wrong input')
    
    ic = 1
    xs = 0d0
    ys = 0d0
    xcum = 0d0
    ycum = 0d0
    
    do ij = 1, JR-1
      do iw = 1, NW
        do is = 1, NS
          do ip = 0, NP
            do ix = 0, NA
              do ia = 0, NA
                if(m(ia, ix, ip, is, iw, ij, it) > 1d-14) then
                  ys(ic) = m(ia, ix, ip, is, iw, ij, it)*pop(ij, it)
                  if (ikind == 0) then
                    xs(ic) = a(ia)
                  elseif (ikind == 1) then
                    xs(ic) = a(ia)*r(it) + pen(ij, ip, it) + w(it)*eff(ij, is)*eta(is, iw)*l(ia, ix, ip, is, iw, ij, it)
                  endif
                  ic = ic + 1
                endif
              enddo ! ia
            enddo ! ix
          enddo ! ip
        enddo ! is
      enddo ! iw
    enddo ! ij
    
    ! get array size and normalize ys
    ICMAX = ic - 1
    ys(1:ICMAX) = ys(1:ICMAX)/sum(ys(1:ICMAX))
     
    ! sort array
    call bubble_sort(xs(1:ICMAX), ys(1:ICMAX))
 
    ! calculate cumulative distributions
    xcum(0) = 0d0
    ycum(0) = 0d0
    do ic = 1, ICMAX
      xcum(ic) = xcum(ic-1) + xs(ic)*ys(ic)
      ycum(ic) = ycum(ic-1) + ys(ic)
    enddo
    
    gini(it, ikind) = 0d0
    do ic = 1, ICMAX
      gini(it, ikind) = gini(it, ikind) + ys(ic)*(xcum(ic-1)+xcum(ic))
    enddo
    gini(it, ikind) = 1d0-gini(it, ikind)/xcum(ICMAX)
      
    percentiles(:, it, ikind) = 0d0
    do ic = ICMAX, 1, -1
      if (1d0-ycum(ic) > 0.01 .and. percentiles(1, it, ikind) <= 0d0) then
        percentiles(1, it, ikind) = (xcum(ICMAX)-xcum(ic))/xcum(ICMAX)
      endif
      if (1d0-ycum(ic) > 0.05 .and. percentiles(2, it, ikind) <= 0d0) then
        percentiles(2, it, ikind) = (xcum(ICMAX)-xcum(ic))/xcum(ICMAX)
      endif
      if (1d0-ycum(ic) > 0.10 .and. percentiles(3, it, ikind) <= 0d0) then
        percentiles(3, it, ikind) = (xcum(ICMAX)-xcum(ic))/xcum(ICMAX)
      endif
      if (1d0-ycum(ic) > 0.20 .and. percentiles(4, it, ikind) <= 0d0) then
        percentiles(4, it, ikind) = (xcum(ICMAX)-xcum(ic))/xcum(ICMAX)
      endif
      if (1d0-ycum(ic) > 0.40 .and. percentiles(5, it, ikind) <= 0d0) then
        percentiles(5, it, ikind) = (xcum(ICMAX)-xcum(ic))/xcum(ICMAX)
      endif
      if (1d0-ycum(ic) > 0.60 .and. percentiles(6, it, ikind) <= 0d0) then
        percentiles(6, it, ikind) = (xcum(ICMAX)-xcum(ic))/xcum(ICMAX)
      endif
    enddo

  end subroutine

  
end module
