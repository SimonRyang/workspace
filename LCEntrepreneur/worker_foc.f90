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
      real*8 :: cons_w, X_plus, tomorrow, varphi_x
      integer :: ixl, ixr

      ! calculate tomorrow's assets
      X_plus  = x_in(1)
      lab_com = x_in(2)

      ! calculate linear interpolation for future part of first order condition
      call linint_Grow(X_plus, X_l, X_u, X_grow, NX, ixl, ixr, varphi_x)

      ! restrict values to grid just in case
      ixl = min(ixl, NX)
      ixr = min(ixr, NX)
      varphi_x = max(min(varphi_x, 1d0),0d0)

      ! get next period value function
      tomorrow = max(varphi_x      *(egam*S(ij_com, ixl, ik_com, iw_com, ie_com, 0))**(1d0/egam) +  &
                     (1d0-varphi_x)*(egam*S(ij_com, ixr, ik_com, iw_com, ie_com, 0))**(1d0/egam), 1d-10)**egam/egam

      ! maximize value function for current worker (next period worker)
      if (ik_com == 0) then

          cons_com = (1d0+r)*a(ia_com) + w*eff(ij_com)*eta(iw_com)*lab_com + pen(ij_com) - X_plus

          if(cons_com <= 0d0)then
             cons_w = -1d-10**egam/egam*(1d0+abs(cons_com))
          else
             cons_w = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + beta*tomorrow)
          endif

      ! maximize value function for current entrepreneur (next period worker)
      else

          cons_com = (1d0+r)*(a(ia_com)-xi*k(ik_com)) + theta(ie_com)*(k(ik_com)**alpha*(eff(ij_com)*lab_com)**(1d0-alpha))**nu + (1d0-delta_k)*k(ik_com) + pen(ij_com) - X_plus

          if(cons_com <= 0d0)then
              cons_w = -1d-10**egam/egam*(1d0+abs(cons_com))
           else
              cons_w = -((cons_com**sigma*(1d0-lab_com)**(1d0-sigma))**egam/egam + beta*tomorrow)
           endif
      endif

  end function

end module
