module globals

    use toolbox

    implicit none

    ! number of years the household lives
    integer, parameter :: JJ = 16

    ! number of years at retirement
    integer, parameter :: JR = 9

    ! number of skills
    integer, parameter :: NS = 3

    ! number of points on the asset grid
    integer, parameter :: NA = 20

    ! number of points on earning points grid
    integer, parameter :: NE = 4

    ! number of points on the housing grid
    integer, parameter :: NH = 20

    ! number of transitory productivity shocks
    integer, parameter :: NP = 3

    ! household preference parameters
    real*8, parameter :: gamma = 0.2d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: nu = 5d0                  ! risk aversion
    real*8, parameter :: beta = 0.96d0**5
    real*8, parameter :: theta = 0.6d0             ! consumption share
    real*8, parameter :: chi = 0.8d0    !0.65      ! ownership preference
    real*8, parameter :: q_1 = 0.0d0                 ! bequest motive
    real*8, parameter :: q_2 = 10.0d0              ! bequest motive


    ! shock process parameters
    real*8, parameter :: sigma_skill = 0.23d0
    real*8, parameter :: rho = 0.98d0
    real*8, parameter :: sigma_eps = 0.05d0

    ! demographic parameters
    real*8, parameter :: n_p   = (1d0+0.005d0)**5-1d0

    ! housing parameters
    real*8, parameter :: delta_h = 1d0-(1d0-0.02d0)**5
    real*8, parameter :: h_min = 2d0                ! minimum house size
    real*8, parameter :: xi = 0.7d0                 ! maximum mortgage share
    real*8, parameter :: phi_1 = 0.03d0
    real*8, parameter :: phi_2 = 0.06d0


    ! production parameters
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta_k = 1d0-(1d0-0.02d0)**5
    real*8, parameter :: Omega = 1.48d0

    ! government policy parameters
    real*8, parameter :: kappa = 0.0d0       ! pension replacement rate
    real*8, parameter :: mu = 0d0            ! paygo pension
    real*8, parameter :: lambda = 0.0d0      ! pension progressivity

    ! size of the total asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 100d0
    real*8, parameter :: a_grow = 0.075d0

    ! size of the liquid asset grid
    real*8, parameter :: al_l    = a_l
    real*8, parameter :: al_u    = a_u*0.5d0
    real*8, parameter :: al_grow = a_grow

    ! size of the earnings point grid
    real*8, parameter :: ep_l    = 0.0d0
    real*8, parameter :: ep_u    = 2d0
    real*8, parameter :: ep_grow = 0.01d0

    ! size of the real estate grid
    real*8, parameter :: h_l = h_min
    real*8, parameter :: h_u = a_u/(1d0-xi)
    real*8, parameter :: h_grow = a_grow

    ! numerical parameters
    integer, parameter :: itermax = 200
    real*8, parameter :: sig = 1d-6
    real*8, parameter :: damp = 0.3d0


    ! the shock process
    real*8 :: dist_skill(NS)
    real*8 :: dist_eta(NS, NP), eta(NS, NP), pi(NS, NP, NP)

    ! demographic and other model parameters
    real*8 :: rpop(0:JJ+1), eff(JJ, NS), psi(JJ+1, NS), workpop

    ! numerical variables
    integer :: iamax(JJ+1), ihmax(JJ+1)
    integer :: ij_com, is_com, ia_com, ie_com, ih_com, ip_com, ia_p_com, iter
    real*8 :: cons_com, ch_com
    real*8 :: DIFF

    ! temporary variables
    real*8 :: aplus_t(JJ+1, NS, 0:NA, 0:NE, 0:NH, NP, 0:1), hplus_t(JJ+1, NS, 0:NA, 0:NE, 0:NH, NP, 0:1)
    real*8 :: ch_t(JJ+1, NS, 0:NA, 0:NE, 0:NH, NP, 0:1), c_t(JJ, NS, 0:NA, 0:NE, 0:NH, NP, 0:1)
    real*8 :: V_t(JJ, NS, 0:NA, 0:NE, 0:NH, NP, 0:1)


    ! individual variables
    real*8 :: a(0:NA), h(0:NH), ep(0:NE), al(0:NA), b(JJ)
    real*8 :: aplus(JJ, NS, 0:NA, 0:NE, 0:NH, NP), hplus(JJ, NS, 0:NA, 0:NE, 0:NH, NP)
    real*8 :: c(JJ, NS, 0:NA, 0:NE, 0:NH, NP), ch(JJ, NS, 0:NA, 0:NE, 0:NH, NP)
    real*8 :: omegaplus(JJ, NS, 0:NA, 0:NE, 0:NH, NP, 0:NA)
    real*8 :: penc(JJ, NS, 0:NA, 0:NE, 0:NH, NP), penp(JJ, NS, 0:NA, 0:NE, 0:NH, NP)
    real*8 :: y(JJ, NS, 0:NA, 0:NE, 0:NH, NP), yg(JJ, NS, 0:NA, 0:NE, 0:NH, NP)
    real*8 :: phi(JJ+1, NS, 0:NA, 0:NE, 0:NH, NP)
    real*8 :: V(JJ, NS, 0:NA, 0:NE, 0:NH, NP)
    real*8 :: SO(JJ, NS, 0:NA, 0:NE, 0:NH, NP, 0:NA), SR(JJ, NS, 0:NA, 0:NE, 0:NH, NP, 0:NA)

    ! cohort average variables
    real*8 ::  a_coh(JJ+1), h_coh(JJ+1), c_coh(JJ), al_coh(JJ+1), ch_coh(JJ), tr_coh(JJ+1)
    real*8 ::  l_coh(JJ), penc_coh(JJ), penp_coh(JJ), v_coh(JJ),  shr(JJ+1, 0:1)


    ! macroeconomic variables
    real*8 :: r, rn, w, wn, p, ph
    real*8 :: ybar
    real*8 :: AA, AAL, BB, BQ, HH, HR, PBEN, PCON
    real*8 :: YY, CC, GG, II, KK, LL, IHO, IHR, TRG


    ! government variables
    integer :: tax
    real*8 :: tauc, taur, tauw, taup
    real*8 :: gy, by, taxrev(4)

contains



    ! solve the household's decision of how much wealth to store in real estate
    subroutine solve_owner(ij, is, ia, ie, ih, ip, ia_p)

        implicit none

        integer, intent(in) :: ij, is, ia, ie, ih, ip, ia_p
        real*8 :: x_in, fret, omega_low

        ! set up communication variables
        ij_com = ij; is_com = is; ia_com = ia; ie_com = ie; ih_com = ih; ip_com = ip; ia_p_com = ia_p

        if (a(ia_p) > ((1d0-xi)*h_min + tr(h(ih), h_min))) then

           ! get lower limit
           omega_low = h_min*(1d0-xi)/a(ia_p)
           x_in = max(omegaplus(ij, is, ia, ie, ih, ip, ia_p-1), omega_low)


           ! solve the household problem using fminsearch
           call fminsearch(x_in, fret, omega_low, 1d0, real_o)

           ! portfolio share for housing
           omegaplus(ij, is, ia, ie, ih, ip, ia_p) = x_in
           SO(ij, is, ia, ie, ih, ip, ia_p) = -fret
       else

           ! portfolio share for housing
           omegaplus(ij, is, ia, ie, ih, ip, ia_p) = 0d0
           SO(ij, is, ia, ie, ih, ip, ia_p) = 1d-10**egam/egam

       endif

    end subroutine



    ! Compute value of household from beeing a future renter
    subroutine solve_renter(ij, is, ia, ie, ih, ip, ia_p)

        implicit none

        integer, intent(in) :: ij, is, ia, ie, ih, ip, ia_p
        integer :: ial, iar, iel, ier, ipp
        real*8 :: al_p, ep_p, varphi_a, varphi_e, dist, EV


        al_p  = a(ia_p) - tr(h(ih),0d0)

        ! calculate linear interpolation for future assets
        call linint_Grow(al_p, al_l, al_u, al_grow, NA, ial, iar, varphi_a)

        ! restrict values to grid just in case
        ial = min(ial, NA)
        iar = min(iar, NA)
        varphi_a = max(min(varphi_a, 1d0),0d0)

        if (al_p >= 0d0) then


            if(ij+1 >= JR)then

                ! get expected future utility from value function (consumption)
                EV = (varphi_a      *(egam*V(ij+1, is, ial, ie, 0, ip))**(1d0/egam) + &
                      (1d0-varphi_a)*(egam*V(ij+1, is, iar, ie, 0, ip))**(1d0/egam))**egam/egam

            else

                do ipp = 1, NP

                    ! get tomorrow's earning points
                    ep_p = ep(ie) + mu*(lambda+(1d0-lambda)*  &
                           min(y(ij, is, ia, ie, ih, ip)/ybar, 2d0))/dble(JR-1d0)

                    ! derive interpolation weights
                    call linint_Grow(ep_p, ep_l, ep_u, ep_grow, NE, iel, ier, varphi_e)

                    iel = min(iel, NE)
                    ier = min(ier, NE)
                    varphi_e = max(min(varphi_e, 1d0), 0d0)

                    ! get distributional weight
                    dist = pi(is, ip, ipp)

                    ! get expected future utility from value function (consumption)
                    EV =  varphi_a*varphi_e*(egam*V(ij+1, is, ial, iel, 0, ipp))**(1d0/egam) + &
                          varphi_a*(1d0-varphi_e)*(egam*V(ij+1, is, ial, ier, 0, ipp))**(1d0/egam) + &
                          (1d0-varphi_a)*varphi_e*(egam*V(ij+1, is, iar, iel, 0, ipp))**(1d0/egam) + &
                          (1d0-varphi_a)*(1d0-varphi_e)*(egam*V(ij+1, is, iar, ier, 0, ipp))**(1d0/egam)

                    EV = EV + dist*EV**egam/egam
                enddo
            endif

            SR(ij, is, ia, ie, ih, ip, ia_p) = psi(ij+1, is)*EV + (1d0-psi(ij+1, is))*q_1*(1d0+max(al_p, 1d-10)/q_2)**egam/egam

        else

            SR(ij, is, ia, ie, ih, ip, ia_p)  = 1d-10**egam/egam

        endif


    end subroutine


    ! solve the household's consumption-savings decision
    subroutine solve_consumption_o(ij, is, ia, ie, ih, ip)

      implicit none

      integer, intent(in) :: ij, is, ia, ie, ih, ip
      real*8 :: x_in, fret, amin, varphi_a, omegah
      integer :: ial, iar

      ! set up communication variables
      ij_com = ij; is_com = is; ia_com = ia; ie_com = ie; ih_com = ih; ip_com = ip

      amin = (1d0-xi)*h_min + tr(h(ih), h_min)

      if (yg(ij, is, ia, ie, ih, ip) < amin)then

          hplus_t(ij, is, ia, ie, ih, ip, 1) = 0d0
          x_in = 0d0
          cons_com = 0d0
          ch_com = 0d0
          fret = -1d-10**egam/egam
      else

          ! get best initial guess from future period
          x_in = max(aplus_t(ij+1, is, ia, ie, ih, ip, 1), amin)

          ! solve the household problem using rootfinding
          call fminsearch(x_in, fret, amin, a_u, cons_o)

          call linint_Grow(x_in, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

          omegah = varphi_a*omegaplus(ij, is, ia, ie, ih, ip, ial) + (1d0-varphi_a)*omegaplus(ij, is, ia, ie, ih, ip, iar)
          hplus_t(ij, is, ia, ie, ih, ip, 1) = omegah*x_in/(1d0-xi)

      endif

      ! copy decisions
      aplus_t(ij, is, ia, ie, ih, ip, 1) = x_in
      c_t(ij, is, ia, ie, ih, ip, 1) = cons_com
      ch_t(ij, is, ia, ie, ih, ip, 1) = ch_com
      V_t(ij, is, ia, ie, ih, ip, 1) = -fret

  end subroutine






  ! solve the household's consumption-savings decision
  subroutine solve_consumption_r(ij, is, ia, ie, ih, ip)

      implicit none

      integer, intent(in) :: ij, is, ia, ie, ih, ip
      real*8 :: x_in, fret, amin

      ! set up communication variables
      ij_com = ij; is_com = is; ia_com = ia; ie_com = ie; ih_com = ih; ip_com = ip

      amin = tr(h(ih), 0d0)

      if(yg(ij, is, ia, ie, ih, ip) < amin)then
          x_in = 0d0
          cons_com = 0d0
          ch_com = 0d0
          fret = -1d-10**egam/egam
      else

          ! get best initial guess from future period
          x_in = max(aplus_t(ij+1, is, ia, ie, ih, ip, 0), amin)

          ! solve the household problem using rootfinding
          call fminsearch(x_in, fret, amin, a_u, cons_r)

      endif

      ! copy decisions
      aplus_t(ij, is, ia, ie, ih, ip, 0) = x_in
      hplus_t(ij, is, ia, ie, ih, ip, 0) = 0d0
      c_t(ij, is, ia, ie, ih, ip, 0) = cons_com
      ch_t(ij, is, ia, ie, ih, ip, 0) = ch_com
      V_t(ij, is, ia, ie, ih, ip, 0) = -fret

  end subroutine


!   Functions for renters and owners


    ! function that computes adjustment cost
    function tr(h, h_p)

       implicit none

       ! input variable
       real*8, intent(in) :: h, h_p
       real*8 :: tr


       tr = 0d0
       if ((1-delta_h)*h > h_p) then

          tr = phi_1*((1-delta_h)*h - h_p)

       elseif (h_p > h) then

          tr = phi_2*(h_p - h)

       endif


    end function



  ! the value function with respect to next period real estate
  function real_o(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: real_o, omega_p, ep_p, h_p, al_p, EV, varphi_a, varphi_h, varphi_e, dist
      integer :: ial, iar, iel, ier, ihl, ihr, ip_p

      ! store real estate share
      omega_p  = x_in

      ! get tomorrow's housing and liquid assets
      h_p = omega_p*a(ia_p_com)/(1d0-xi)
      al_p = (1d0-omega_p)*a(ia_p_com) - tr(h(ih_com), h_p)

      ! derive interpolation weights
      call linint_Grow(h_p, h_l, h_u, h_grow, NH-1, ihl, ihr, varphi_h)
      call linint_Grow(al_p, al_l, al_u, al_grow, NA, ial, iar, varphi_a)

      ! restrict values to grid just in case (note: grid for owners starts at 1!)
      ihl = min(ihl+1, NH)
      ihr = min(ihr+1, NH)
      varphi_h = max(min(varphi_h, 1d0), 0d0)

      ial = min(ial, NA)
      iar = min(iar, NA)
      varphi_a = max(min(varphi_a, 1d0), 0d0)

      if (al_p > 0d0) then

         if(ij_com+1 >= JR)then

             ! get expected future utility from value function (consumption)
             if(varphi_a <= varphi_h)then
                 EV = (varphi_a           *(egam*V(ij_com+1, is_com, ial, ie_com, ihl, ip_com))**(1d0/egam) + &
                       (varphi_h-varphi_a)*(egam*V(ij_com+1, is_com, iar, ie_com, ihl, ip_com))**(1d0/egam) + &
                       (1d0-varphi_h)     *(egam*V(ij_com+1, is_com, iar, ie_com, ihr, ip_com))**(1d0/egam))**egam/egam
             else
                 EV = (varphi_h           *(egam*V(ij_com+1, is_com, ial, ie_com, ihl, ip_com))**(1d0/egam) + &
                       (varphi_a-varphi_h)*(egam*V(ij_com+1, is_com, ial, ie_com, ihr, ip_com))**(1d0/egam) + &
                       (1d0-varphi_a)     *(egam*V(ij_com+1, is_com, iar, ie_com, ihr, ip_com))**(1d0/egam))**egam/egam
             endif

         else

             ! get tomorrow's earning points
             ep_p = ep(ie_com) + mu*(lambda+(1d0-lambda)*  &
                    min(y(ij_com, is_com, ia_com, ie_com, ih_com, ip_com)/ybar, 2d0))/dble(JR-1d0)

             ! derive interpolation weights
             call linint_Grow(ep_p, ep_l, ep_u, ep_grow, NE, iel, ier, varphi_e)

             iel = min(iel, NE)
             ier = min(ier, NE)
             varphi_e = max(min(varphi_e, 1d0), 0d0)

             do ip_p = 1, NP

                 ! get distributional weight
                 dist = pi(is_com, ip_com, ip_p)

                 ! get expected future utility from value function (consumption)
                 EV =  varphi_a*varphi_e*varphi_h*(egam*V(ij_com+1, is_com, ial, iel, ihl, ip_p))**(1d0/egam) + &
                       varphi_a*varphi_e*(1d0-varphi_h)*(egam*V(ij_com+1, is_com, ial, iel, ihr, ip_p))**(1d0/egam) + &
                       varphi_a*(1d0-varphi_e)*varphi_h*(egam*V(ij_com+1, is_com, ial, ier, ihl, ip_p))**(1d0/egam) + &
                       varphi_a*(1d0-varphi_e)*(1d0-varphi_h)*(egam*V(ij_com+1, is_com, ial, ier, ihr, ip_p))**(1d0/egam) + &
                       (1d0-varphi_a)*varphi_e*varphi_h*(egam*V(ij_com+1, is_com, iar, iel, ihl, ip_p))**(1d0/egam) + &
                       (1d0-varphi_a)*varphi_e*(1d0-varphi_h)*(egam*V(ij_com+1, is_com, iar, iel, ihr, ip_p))**(1d0/egam) + &
                       (1d0-varphi_a)*(1d0-varphi_e)*varphi_h*(egam*V(ij_com+1, is_com, iar, ier, ihl, ip_p))**(1d0/egam) + &
                       (1d0-varphi_a)*(1d0-varphi_e)*(1d0-varphi_h)*(egam*V(ij_com+1, is_com, iar, ier, ihr, ip_p))**(1d0/egam)

                 EV = EV + dist*EV**egam/egam
             enddo
         endif

         real_o = -(psi(ij_com+1, is_com)*EV + (1d0-psi(ij_com+1, is_com))*q_1*(1d0+max(al_p+h_p, 1d-10)/q_2)**egam/egam)

      else

         real_o = 1000d0*abs(al_p)

      endif

  end function


  ! the optimizing function regarding saving/consumption of next period owner
  function cons_o(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: cons_o, a_plus, tomorrow, varphi_a
      integer :: ial, iar

      ! calculate tomorrow's assets
      a_plus  = x_in

      ! calculate linear interpolation for future assets
      call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

      ! restrict values to grid just in case
      ial = min(ial, NA)
      iar = min(iar, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)
write(*,*)'1'
      ! get next period value function
      tomorrow = max(varphi_a*(egam*SO(ij_com, is_com, ia_com, ie_com, ih_com, ip_com, ial))**(1d0/egam) +  &
               (1d0-varphi_a)*(egam*SO(ij_com, is_com, ia_com, ie_com, ih_com, ip_com, iar))**(1d0/egam), 1d-10)**egam/egam
write(*,*)'2'
      ! maximize value function for current renter
      if (ih_com == 0) then

           ch_com = (1d0-theta)*(yg(ij_com, is_com, ia_com, ie_com, ih_com, ip_com)- a_plus)/ph
           cons_com = (yg(ij_com, is_com, ia_com, ie_com, ih_com, ip_com) - a_plus - ph*ch_com)/p
           if(ch_com <= 0d0 .or. cons_com <= 0d0)then
               cons_o = -1d-10**egam/egam*(1d0+abs(cons_com))
           else
               cons_o = -((cons_com**theta*(chi*ch_com)**(1d0-theta))**egam/egam + beta*tomorrow)
           endif

      ! maximize value function for current owner
      else

           cons_com = (yg(ij_com, is_com, ia_com, ie_com, ih_com, ip_com) - a_plus)/p
           ch_com = h(ih_com)

           if(cons_com <= 0d0)then
              cons_o = -1d-10**egam/egam*(1d0+abs(cons_com))
           else
              cons_o = -((cons_com**theta*h(ih_com)**(1d0-theta))**egam/egam + beta*tomorrow)
           endif

      endif

  end function


  ! the optimizing function regarding saving/consumption of next period renter
  function cons_r(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: cons_r, a_plus, tomorrow, varphi_a
      integer :: ial, iar

      ! calculate tomorrow's assets
      a_plus  = x_in

      ! calculate linear interpolation for future part of first order condition
      call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

      ! restrict values to grid just in case
      ial = min(ial, NA)
      iar = min(iar, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)

      ! get next period value function
      tomorrow = max(varphi_a*(egam*SR(ij_com, is_com, ia_com, ie_com, ih_com, ip_com, ial))**(1d0/egam)   +  &
               (1d0-varphi_a)*(egam*SR(ij_com, is_com, ia_com, ie_com, ih_com, ip_com, iar))**(1d0/egam), 1d-10)**egam/egam

      ! maximize value function for current renter
      if (ih_com == 0) then

          ch_com = (1d0-theta)*(yg(ij_com, is_com, ia_com, ie_com, ih_com, ip_com) - a_plus)/ph
          cons_com = (yg(ij_com, is_com, ia_com, ie_com, ih_com, ip_com) - a_plus - ph*ch_com)/p

          if((ch_com <= 0d0) .or. (cons_com <= 0d0))then
             cons_r = -1d-10**egam/egam*(1d0+abs(cons_com))
          else
             cons_r = -((cons_com**theta*(chi*ch_com)**(1d0-theta))**egam/egam + beta*tomorrow)
          endif

      ! maximize value function for current owner
      else

          cons_com = (yg(ij_com, is_com, ia_com, ie_com, ih_com, ip_com) - a_plus)/p
          ch_com = h(ih_com)

          if(cons_com <= 0d0)then
              cons_r = -1d-10**egam/egam*(1d0+abs(cons_com))
          else
              cons_r = -((cons_com**theta*h(ih_com)**(1d0-theta))**egam/egam + beta*tomorrow)
          endif
      endif

  end function



    ! Converts sigma of an one-year process into sigma of an five-year process
    function sigma5(rho, sigma)

       implicit none

       real*8, intent(in) :: rho, sigma
       real*8 :: sigma5

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

end module
