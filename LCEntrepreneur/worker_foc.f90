!##############################################################################
! MODULE renter_foc
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################
module worker_foc

    use toolbox
    use globals

    implicit none

contains

  ! the first order condition regarding consumption
  function cons_w(x_in)

      implicit none

      ! input variable
      real*8, intent(in) :: x_in(:)

      ! variable declarations
      real*8 :: cons_w, X_plus, income, tomorrow, varphi_x, varphi_p
      integer :: io, ixl_p, ixr_p, ipl_p, ipr_p

      ! calculate tomorrow's assets
      X_plus  = x_in(1)
      lab_com = x_in(2)

      ! current occupation
      io = abs(ik_com > 0)

      income = (1d0-dble(io))*w*eff(ij_com)*eta(iw_com)*lab_com + &
               dble(io)*theta(ie_com)*(k(ik_com)**alpha*(eff(ij_com)*lab_com)**(1d0-alpha))**nu + (1d0-delta_k)*k(ik_com)

      p_plus_com = (p(ip_com)*dble(ij_com-1) + (1d0-(1d0-phi)*dble(io))*mu*(lambda + (1d0-lambda)*min(w*eff(ij_com)*eta(iw_com)*lab_com, p_u)))/dble(ij_com)


      cons_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + income + pen(ij_com, ip_com) - X_plus

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
         cons_w = -1d-10**egam/egam*(1d0+abs(cons_com))
      else
         cons_w = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + beta*tomorrow)
      endif

  end function

end module
