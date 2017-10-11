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
  function inv_o(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in

      ! variable declarations
      real*8 :: inv_o, ad_p, a_p, k_p, EV_temp, S_temp, omega_k, varphi_k, varphi_a, a_temp
      integer :: ikl_p, ikr_p, ial_p, iar_p

      ! store real estate share
      omega_k  = x_in

      ! determine future liquid wealth and future downpayment
      ad_p = (1d0-xi)*k_min + omega_k*(X(ix_p_com)-(1d0-xi)*k_min)
      k_p =  ad_p/(1d0-xi)
      a_temp = X(ix_p_com) - ad_p - tr(k(ik_com), k_p)
      a_p = max(a_temp, 0d0)

      ! derive interpolation weights
      call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial_p, iar_p, varphi_a)
      call linint_Grow(k_p, k_l, k_u, k_grow, NK-1, ikl_p, ikr_p, varphi_k)

      ! restrict values to grid just in case
      ial_p = min(ial_p, NA)
      iar_p = min(iar_p, NA)
      varphi_a = max(min(varphi_a, 1d0),0d0)

      ! restrict values to grid just in case
      ikl_p = min(ikl_p+1, NK)
      ikr_p = min(ikr_p+1, NK)
      varphi_k = max(min(varphi_k, 1d0), 0d0)

      S_temp = (1d0-psi(ij_com+1))*mu_b*max(X(ix_p_com), 1d-10)**egam/egam

      ! get optimal investment strategy
      if(varphi_a <= varphi_k)then
          EV_temp = varphi_a           *(egam*EV(ij_com+1, ial_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                    (varphi_k-varphi_a)*(egam*EV(ij_com+1, iar_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                    (1d0-varphi_k)     *(egam*EV(ij_com+1, iar_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam)
      else
          EV_temp = varphi_k           *(egam*EV(ij_com+1, ial_p, ip_p_com, ikl_p, iw_com, ie_com))**(1d0/egam) + &
                    (varphi_a-varphi_k)*(egam*EV(ij_com+1, ial_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam) + &
                    (1d0-varphi_a)     *(egam*EV(ij_com+1, iar_p, ip_p_com, ikr_p, iw_com, ie_com))**(1d0/egam)
      endif

      S_temp = S_temp + psi(ij_com+1)*EV_temp**egam/egam

      inv_o = - S_temp + 100d0*abs(a_p-a_temp)

  end function


  ! the first order condition regarding consumption
  function cons_e(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in(:)

      ! variable declarations
      real*8 :: cons_e, X_plus, ind_o, income, tomorrow, varphi_x, varphi_p
      integer :: ixl_p, ixr_p, ipl_p, ipr_p

      ! calculate tomorrow's assets
      X_plus  = x_in(1)
      lab_com = max(x_in(2), 0d0)

      ! current occupation
      ind_o = abs(dble(ik_com > 0))

      income = (1d0-ind_o)*w*eff(ij_com)*eta(iw_com)*lab_com + &
               ind_o*theta(ie_com)*(k(ik_com)**alpha*(eff(ij_com)*lab_com)**(1d0-alpha))**nu + (1d0-delta_k)*k(ik_com)

      p_plus_com = (p(ip_com)*dble(ij_com-1) + (1d0-(1d0-phi)*ind_o)*mu*(lambda + (1d0-lambda)*min(w*eff(ij_com)*eta(iw_com)*lab_com, p_u)))/dble(ij_com)


      cons_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + income + pen(ij_com, ip_com) - X_plus


      ! 
      ! ! maximize value function for current worker (next period entrepreneur)
      ! if (ik_com == 0) then
      !   cons_com = (1d0+r)*a(ia_com) + w*eff(ij_com)*eta(iw_com)*lab_com + pen(ij_com, ip_com) - X_plus
      !   p_plus_com = (p(ip_com)*dble(ij_com-1) + mu*(lambda + (1d0-lambda)*min(w*eff(ij_com)*eta(iw_com)*lab_com, p_u)))/dble(ij_com)
      !
      ! ! maximize value function for current owner (next period owner)
      ! else
      !   cons_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + theta(ie_com)*(k(ik_com)**alpha*(eff(ij_com)*lab_com)**(1d0-alpha))**nu + (1d0-delta_k)*k(ik_com) + pen(ij_com, ip_com) - X_plus
      !   p_plus_com = (p(ip_com)*dble(ij_com-1) + phi*mu*(lambda + (1d0-lambda)*min(theta(ie_com)*(k(ik_com)**alpha*(eff(ij_com)*lab_com)**(1d0-alpha))**nu, p_u)))/dble(ij_com)
      !
      ! endif

      if (ij_com >= JR) then
        p_plus_com = p(ip_com)
      endif

      ! calculate linear interpolation for future part of first order condition
      call linint_Grow(X_plus, X_l, X_u, X_grow, NX, ixl_p, ixr_p, varphi_x)
      call linint_Equi(p_plus_com, p_l, p_u, NP, ipl_p, ipr_p, varphi_p)

      ! restrict values to grid just in case
      ixl_p = min(ixl_p, NX)
      ixr_p = min(ixr_p, NX)
      varphi_x = max(min(varphi_x, 1d0),0d0)

      ! restrict values to grid just in case
      ipl_p = min(ipl_p, NP)
      ipr_p = min(ipr_p, NP)
      varphi_p = max(min(varphi_p, 1d0),0d0)

      ! get next period value function
      tomorrow = max(varphi_x*varphi_p              *(egam*S(ij_com, ixl_p, ipl_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                     varphi_x*(1d0-varphi_p)        *(egam*S(ij_com, ixl_p, ipr_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                     (1d0-varphi_x)*varphi_p        *(egam*S(ij_com, ixr_p, ipl_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam) +  &
                     (1d0-varphi_x)*(1d0-varphi_p)  *(egam*S(ij_com, ixr_p, ipr_p, ik_com, iw_com, ie_com, io_p_com))**(1d0/egam), 1d-10)**egam/egam

      if(cons_com <= 0d0)then
         cons_e = -1d-18**egam/egam*(1d0+abs(cons_com))
      elseif(lab_com < 0d0) then
        cons_e = -1d-18**egam/egam*(1d0+abs(lab_com))
      elseif(lab_com >= 1d0) then
        cons_e = -1d-18**egam/egam*lab_com
      else
         cons_e = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + beta*tomorrow)
      endif

  end function

end module
