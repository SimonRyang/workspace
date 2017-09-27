!##############################################################################
! MODULE renter_foc
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################
module renter_foc

    use toolbox
    use globals

    implicit none

contains

  ! the first order condition regarding portfolio choice
  function liq_port_r(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: liq_port_r, earnings, omega_p, R_port, X_p, varphi_x, dist, Q_temp, EV
      integer :: ic_p, ixl, ixr, iw, isr

      ! store portfolio share
      omega_p  = x_in

      ! initialize Q_temp
      Q_temp = (1d0 - psi(ij_com+1, ic_com))*nu**(1d0/gamma)*max(l(il_p_com), 1d-10)**egam/egam

      if(ij_com+1 >= JR)then
          do isr = 1, NSR
              do ic_p = 0, NC

                  ! get return on the portfolio
                  R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))

                  ! get tomorrow's cash-on-hand (epsilon^+ = 0)
                  X_p = R_port*l(il_p_com) + pen(ij_com+1) - dble(ic_p)*ltc(ij_com+1)

                  ! derive interpolation weights
                  call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_x)

                  ! restrict values to grid just in case
                  ixl = min(ixl, NX)
                  ixr = min(ixr, NX)
                  varphi_x = max(min(varphi_x, 1d0), 0d0)

                  ! get distributional weight
                  dist = dist_epsvtheta(isr)*dist_ltc(ij_com+1, ic_com, ic_p)

                  EV = varphi_X*(egam*V(ij_com+1, ixl, ic_p, 0))**(1d0/egam) + &
                           (1d0-varphi_X)*(egam*V(ij_com+1, ixr, ic_p, 0))**(1d0/egam)

                  Q_temp = Q_temp + dist*psi(ij_com+1, ic_com)*EV**egam/egam

              enddo
          enddo
      else
          do iw = 1, NW
              do isr = 1, NSR

                  ! get return on the portfolio
                  R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))

                  ! derive labor earnings
                  earnings  = w*eff(ij_com+1)*zeta(iw)

                  ! get tomorrow's cash on hand
                  X_p = R_port*l(il_p_com)/eps(isr) + earnings

                  ! derive interpolation weights
                  call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_x)

                  ! restrict values to grid just in case
                  ixl = min(ixl, NX)
                  ixr = min(ixr, NX)
                  varphi_x = max(min(varphi_x, 1d0), 0d0)

                  ! get distributional weight
                  dist = dist_zeta(iw)*dist_epsvtheta(isr)

                  ! get expected future utility from value function (consumption)
                  EV = varphi_X*(egam*V(ij_com+1, ixl, 0, 0))**(1d0/egam) + &
                           (1d0-varphi_X)*(egam*V(ij_com+1, ixr, 0, 0))**(1d0/egam)

                  Q_temp = Q_temp + dist*psi(ij_com+1, ic_com)*(eps(isr)*EV)**egam/egam

              enddo
          enddo
      endif

      liq_port_r = - Q_temp

  end function


  ! the first order condition regarding consumption
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
      tomorrow = max(varphi_a*(egam*S(ij_com, ic_com, ial, ih_com, 0))**(1d0/egam)   +  &
                     (1d0-varphi_a)*(egam*S(ij_com, ic_com, iar, ih_com, 0))**(1d0/egam), 1d-10)**egam/egam
      
      ! maximize value function for current renter (next period renter)
      if (ih_com == 0) then
 
          ch_com = (1d0-theta)*(X(ix_com) - a_plus)/ph
          cons_com = X(ix_com) - a_plus - ph*ch_com
 
          if(cons_com <= 0d0)then
             cons_r = -1d-10**egam/egam*(1d0+abs(cons_com))
          else
             cons_r = -((cons_com**theta*(chi*ch_com)**(1d0-theta))**egam/egam + beta*tomorrow)
          endif
 
      ! maximize value function for current owner (next period renter)
      else
      
          cons_com = X(ix_com) - a_plus 
          ch_com = h(ih_com)
          
          if(cons_com <= 0d0)then
              cons_r = -1d-10**egam/egam*(1d0+abs(cons_com))
           else
              if (ic_com == 0) then 
                 cons_r = -((cons_com**theta*h(ih_com)**(1d0-theta))**egam/egam + beta*tomorrow)
              else
                 cons_r = -((cons_com**theta*(chi*h(ih_com))**(1d0-theta))**egam/egam + beta*tomorrow)
              endif   
           endif
      endif 

  end function

end module
