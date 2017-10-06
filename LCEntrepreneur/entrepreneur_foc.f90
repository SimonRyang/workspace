!##############################################################################
! MODULE owner_foc
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################

module entrepreneur_foc

    use toolbox
    use globals

    implicit none

contains

  ! the first order condition with respect to next period real estate
  function real_o(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: real_o, ad_p, a_p, k_p, EV_temp, S_temp, omega_k, varphi_k, varphi_a, a_temp
      integer :: ikl, ikr, ial, iar

      ! store real estate share
      omega_k  = x_in

      ! determine future liquid wealth and future downpayment
      ad_p = (1d0-xi)*k_min + omega_k*(X(ix_p_com)-(1d0-xi)*k_min)
      k_p =  ad_p/(1d0-xi)
      a_temp = X(ix_p_com) - ad_p - tr(k(ik_com), k_p)
      a_p = max(a_temp, 0d0)

      ! derive interpolation weights
      call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
      call linint_Grow(k_p, k_l, k_u, k_grow, NK, ikl, ikr, varphi_k)

      ! restrict values to grid just in case
      ial = min(ial, NA)
      iar = min(iar, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)

      ! restrict values to grid just in case
      ikl = min(ikl, NK)
      ikr = min(ikr, NK)
      varphi_k = max(min(varphi_k, 1d0), 0d0)

      S_temp = (1d0-psi(ij_com+1))*mu_b*max(X(ix_p_com), 1d-10)**egam/egam

      ! get optimal investment strategy
      if(varphi_a <= varphi_k)then
          EV_temp = varphi_a           *(egam*EV(ij_com+1, ial, ikl, iw_com, ie_com))**(1d0/egam) + &
                    (varphi_k-varphi_a)*(egam*EV(ij_com+1, iar, ikl, iw_com, ie_com))**(1d0/egam) + &
                    (1d0-varphi_k)     *(egam*EV(ij_com+1, iar, ikr, iw_com, ie_com))**(1d0/egam)
      else
          EV_temp = varphi_k           *(egam*EV(ij_com+1, ial, ikl, iw_com, ie_com))**(1d0/egam) + &
                    (varphi_a-varphi_k)*(egam*EV(ij_com+1, ial, ikr, iw_com, ie_com))**(1d0/egam) + &
                    (1d0-varphi_a)     *(egam*EV(ij_com+1, iar, ikr, iw_com, ie_com))**(1d0/egam)
      endif

      S_temp = S_temp + psi(ij_com+1)*EV_temp**egam/egam

      real_o = - S_temp + 100d0*abs(a_p-a_temp)

  end function


  ! the first order condition regarding consumption
  function cons_e(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: cons_e, X_plus, tomorrow, varphi_x
      integer :: ixl, ixr

      ! calculate tomorrow's assets
      X_plus  = x_in

      ! calculate linear interpolation for future part of first order condition
      call linint_Grow(X_plus, X_l, X_u, X_grow, NX, ixl, ixr, varphi_x)

      ! restrict values to grid just in case
      ixl = min(ixl, NX)
      ixr = min(ixr, NX)
      varphi_x = max(min(varphi_x, 1d0),0d0)

      ! get next period value function
      tomorrow = max(varphi_x      *(egam*S(ij_com, ixl, ik_com, iw_com, ie_com, 1))**(1d0/egam) +  &
                     (1d0-varphi_x)*(egam*S(ij_com, ixr, ik_com, iw_com, ie_com, 1))**(1d0/egam), 1d-10)**egam/egam


      ! maximize value function for current worker (next period entrepreneur)
      if (ik_com == 0) then

        cons_com = (1d0+r)*a(ia_com) + w*eff(ij_com)*eta(iw_com) + pen(ij_com) - X_plus

         if(cons_com <= 0d0)then
             cons_e = -1d-10**egam/egam*(1d0+abs(cons_com))
         else
             cons_e = -(cons_com**egam/egam + beta*tomorrow)
         endif

      ! maximize value function for current owner (next period owner)
      else

        cons_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + theta(ie_com)*(k(ik_com)**alpha*eff(ij_com)**(1d0-alpha))**nu + (1d0-delta_k)*k(ik_com) + pen(ij_com) - X_plus

         if(cons_com <= 0d0)then
            cons_e = -1d-10**egam/egam*(1d0+abs(cons_com))
         else
            cons_e = -(cons_com**egam/egam + beta*tomorrow)
         endif

      endif

  end function

end module
