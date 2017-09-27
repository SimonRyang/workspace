!##############################################################################
! MODULE owner_foc
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################

module owner_foc

    use toolbox
    use globals

    implicit none

contains

  ! the first order condition regarding portfolio choice
  function liq_port_o(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: liq_port_o, earnings, h_p, omega_p, R_port, X_p, varphi_h, varphi_x, dist, Q_temp, EV
      integer :: ic_p, ihl, ihr, ixl, ixr, isr, iw

      ! store portfolio share
      omega_p  = x_in

      ! initialize Q_temp
      Q_temp = (1d0 - psi(ij_com+1, ic_com))*nu**(1d0/gamma)*max(l(il_p_com)+d(id_p_com),1d-10)**egam/egam

      if(ij_com+1 >= JR)then
          do isr = 1, NSR
              do ic_p = 0, NC

                  ! get return on the portfolio
                  if (omega_p >= 0d0) then 
                     R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
                  else
                     R_port = 1d0 + r_f + rp
                  endif      

                  ! get tomorrow's cash-on-hand (epsilon^+ = 0)
                  h_p = d(id_p_com)/(1d0-xi)
                  X_p = R_port*(l(il_p_com) - xi*h_p) + pen(ij_com+1)+(1d0-delta_h)*h_p - dble(ic_p)*ltc(ij_com+1)

                  ! derive interpolation weights
                  call linint_Grow(h_p, h_l, h_u, h_grow, NH-1, ihl, ihr, varphi_h)
                  call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_x)

                  ! restrict values to grid just in case (achtung: grid für homeowner beginnt bei 1!!!)
                  ihl = min(ihl+1, NH)
                  ihr = min(ihr+1, NH)
                  varphi_h = max(min(varphi_h, 1d0), 0d0)

                  ! restrict values to grid just in case
                  ixl = min(ixl, NX)
                  ixr = min(ixr, NX)
                  varphi_x = max(min(varphi_x, 1d0), 0d0)

                  ! get distributional weight
                  dist = dist_epsvtheta(isr)*dist_ltc(ij_com+1, ic_com, ic_p)

                  ! get expected future utility from value function (consumption)
                  if(varphi_X <= varphi_h)then
                      EV = varphi_X           *(egam*V(ij_com+1, ixl, ic_p, ihl))**(1d0/egam) + &
                           (varphi_h-varphi_X)*(egam*V(ij_com+1, ixr, ic_p, ihl))**(1d0/egam) + &
                           (1d0-varphi_h)     *(egam*V(ij_com+1, ixr, ic_p, ihr))**(1d0/egam)
                  else
                      EV = varphi_h           *(egam*V(ij_com+1, ixl, ic_p, ihl))**(1d0/egam) + &
                           (varphi_X-varphi_h)*(egam*V(ij_com+1, ixl, ic_p, ihr))**(1d0/egam) + &
                           (1d0-varphi_X)     *(egam*V(ij_com+1, ixr, ic_p, ihr))**(1d0/egam)
                  endif

                  Q_temp = Q_temp + dist*psi(ij_com+1, ic_com)*EV**egam/egam

              enddo
          enddo
      else
          do iw = 1, NW
              do isr = 1, NSR

                  ! get return on the portfolio
                  if (omega_p >= 0d0) then 
                     R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
                  else
                     R_port = 1d0 + r_f + rp
                  endif      

                  ! derive labor earnings
                  earnings  = w*eff(ij_com+1)*zeta(iw)

                  ! get tomorrow's housing
                  h_p = d(id_p_com)/(eps(isr)*(1d0-xi))

                  ! get tomorrow's cash on hand
                  X_p = R_port*(l(il_p_com) - xi*d(id_p_com)/(1d0-xi))/eps(isr) + earnings + (1d0 - delta_h)*h_p

                  ! derive interpolation weights
                  call linint_Grow(h_p, h_l, h_u, h_grow, NH-1, ihl, ihr, varphi_h)
                  call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_x)

                  ! restrict values to grid just in case
                  ihl = min(ihl+1, NH)
                  ihr = min(ihr+1, NH)
                  varphi_h = max(min(varphi_h, 1d0), 0d0)

                  ! restrict values to grid just in case
                  ixl = min(ixl, NX)
                  ixr = min(ixr, NX)
                  varphi_x = max(min(varphi_x, 1d0), 0d0)

                  ! get distributional weight
                  dist = dist_zeta(iw)*dist_epsvtheta(isr)

                  ! get expected future utility from value function (consumption)
                  if(varphi_X <= varphi_h)then
                      EV = varphi_X           *(egam*V(ij_com+1, ixl, 0, ihl))**(1d0/egam) + &
                           (varphi_h-varphi_X)*(egam*V(ij_com+1, ixr, 0, ihl))**(1d0/egam) + &
                           (1d0-varphi_h)     *(egam*V(ij_com+1, ixr, 0, ihr))**(1d0/egam)
                  else
                      EV = varphi_h           *(egam*V(ij_com+1, ixl, 0, ihl))**(1d0/egam) + &
                           (varphi_X-varphi_h)*(egam*V(ij_com+1, ixl, 0, ihr))**(1d0/egam) + &
                           (1d0-varphi_X)     *(egam*V(ij_com+1, ixr, 0, ihr))**(1d0/egam)
                  endif

                  Q_temp = Q_temp+ dist*psi(ij_com+1, ic_com)*(eps(isr)*EV)**egam/egam

              enddo
          enddo
      endif

      liq_port_o = -Q_temp

  end function


  ! the first order condition with respect to next period real estate
  function real_o(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: real_o, ad_p, aliq_p, h_p, S_temp, omega_h, varphi_d, varphi_l, aliq_temp
      integer :: idl, idr, ill, ilr

      ! store real estate share
      omega_h  = x_in

      ! determine future liquid wealth and future downpayment
      ad_p = (1d0-xi)*h_min + omega_h*(a(ia_p_com)-(1d0-xi)*h_min)
      h_p =  ad_p/(1d0-xi)
      aliq_temp = a(ia_p_com)- ad_p - tr(h(ih_com), h_p)
      aliq_p = max(aliq_temp, 0d0)

      ! derive interpolation weights
      call linint_Grow(aliq_p, l_l, l_u, l_grow, NL, ill, ilr, varphi_l)
      call linint_Grow(ad_p, d_l, d_u, d_grow, ND-1, idl, idr, varphi_d)

      ! restrict values to grid just in case
      ill = min(ill, NL)
      ilr = min(ilr, NL)
      varphi_l = max(min(varphi_l, 1d0),0d0)

      ! restrict values to grid just in case
      idl = min(idl+1, ND)
      idr = min(idr+1, ND)
      varphi_d = max(min(varphi_d, 1d0), 0d0)

      ! get optimal investment strategy
      if(varphi_l <= varphi_d)then
          S_temp = varphi_l           *(egam*Q(ij_com, ic_com, ill, idl))**(1d0/egam) + &
                   (varphi_d-varphi_l)*(egam*Q(ij_com, ic_com, ilr, idl))**(1d0/egam) + &
                   (1d0-varphi_d)     *(egam*Q(ij_com, ic_com, ilr, idr))**(1d0/egam)
      else
          S_temp = varphi_d           *(egam*Q(ij_com, ic_com, ill, idl))**(1d0/egam) + &
                   (varphi_l-varphi_d)*(egam*Q(ij_com, ic_com, ill, idr))**(1d0/egam) + &
                   (1d0-varphi_l)     *(egam*Q(ij_com, ic_com, ilr, idr))**(1d0/egam)
      endif

      real_o = - S_temp**egam/egam + 100d0*abs(aliq_p-aliq_temp)

  end function


  ! the first order condition regarding consumption
  function cons_o(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: cons_o, a_plus, tomorrow, varphi_a
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
      tomorrow = max(varphi_a*(egam*S(ij_com, ic_com, ial, ih_com, 1))**(1d0/egam) +  &
                     (1d0-varphi_a)*(egam*S(ij_com, ic_com, iar, ih_com, 1))**(1d0/egam), 1d-10)**egam/egam


      ! maximize value function for current renter (next period owner)
      if (ih_com == 0) then

           ch_com = (1d0-theta)*(X(ix_com)- a_plus)/ph
           cons_com = X(ix_com) - a_plus - ph*ch_com 
           if(cons_com <= 0d0)then
               cons_o = -1d-10**egam/egam*(1d0+abs(cons_com))
           else
               cons_o = -((cons_com**theta*(chi*ch_com)**(1d0-theta))**egam/egam + beta*tomorrow)
           endif                            

      ! maximize value function for current owner (next period owner)
      else

           cons_com = X(ix_com) - a_plus
           ch_com = h(ih_com)

           if(cons_com <= 0d0)then
              cons_o = -1d-10**egam/egam*(1d0+abs(cons_com))
           else
              if (ic_com == 0) then 
                 cons_o = -((cons_com**theta*h(ih_com)**(1d0-theta))**egam/egam + beta*tomorrow)
              else
                 cons_o = -((cons_com**theta*(chi*h(ih_com))**(1d0-theta))**egam/egam + beta*tomorrow)
              endif   
           endif

      endif

  end function

end module
