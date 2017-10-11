!##############################################################################
! MODULE renter_solve
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################

module worker_solve

    use toolbox
    use globals
    use worker_foc

    implicit none

contains

    ! solve the household's decision of how much wealth to invest into capital
    subroutine solve_worker(ij, ix_p, ip_p, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, ix_p, ip_p, ik, iw, ie
        integer :: ial, iar
        real*8 :: a_plus, EV_temp, S_temp, varphi_a

        a_plus  = X(ix_p)

       ! calculate linear interpolation for future assets
       call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

       ! restrict values to grid just in case
       ial = min(ial, NA)
       iar = min(iar, NA)
       varphi_a = max(min(varphi_a, 1d0),0d0)

       S_temp = (1d0-psi(ij+1))*mu_b*max(a_plus, 1d-10)**egam/egam

       EV_temp = varphi_a      *(egam*EV(ij+1, ial, ip_p, 0, iw, ie))**(1d0/egam) + &
                 (1d0-varphi_a)*(egam*EV(ij+1, iar, ip_p, 0, iw, ie))**(1d0/egam)

       S_temp = S_temp + psi(ij+1)*EV_temp**egam/egam

       S(ij, ix_p, ip_p, ik, iw, ie, 0) = S_temp

    end subroutine


  ! solve the household's consumption-savings decision
  subroutine solve_consumption_w(ij, ia, ip, ik, iw, ie)

      implicit none

      integer, intent(in) :: ij, ia, ip, ik, iw, ie
      real*8 :: x_in(2), fret

      ! set up communication variables
      ij_com = ij; ia_com = ia; ip_com = ip; ik_com = ik; iw_com = iw; ie_com = ie

      ! get best initial guess from future period
      x_in(1) = max(X_plus_t(ij+1, ia, ip, ik, iw, ie, 0), 1d-4)
      x_in(2) = max(l_t(ij+1, ia, ip, ik, iw, ie, 0), 0.33d0)

      ! solve the household problem using rootfinding
      call fminsearch(x_in, fret, (/X_l, 0d0/), (/X_u, 0.8d0/), cons_w)

      ! copy decisions
      X_plus_t(ij, ia, ip, ik, iw, ie, 0) = x_in(1)
      a_plus_t(ij, ia, ip, ik, iw, ie, 0) = x_in(1)
      p_plus_t(ij, ia, ip, ik, iw, ie, 0) = p_plus_com
      k_plus_t(ij, ia, ip, ik, iw, ie, 0) = 0d0
      c_t(ij, ia, ip, ik, iw, ie, 0) = cons_com
      l_t(ij, ia, ip, ik, iw, ie, 0) = x_in(2)
      V_t(ij, ia, ip, ik, iw, ie, 0) = -fret

  end subroutine


end module
