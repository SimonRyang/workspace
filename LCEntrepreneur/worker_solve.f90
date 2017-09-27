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
    subroutine solve_worker(ij, ix_p, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, ix_p, ik, iw, ie
        integer :: ial, iar, iw_p, ie_p
        real*8 :: a_p, EV, S_temp, varphi_a

        a_p  = X(ix_p)

       ! calculate linear interpolation for future assets
       call linint_Grow(a_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

       ! restrict values to grid just in case
       ial = min(ial, NA)
       iar = min(iar, NA)
       varphi_a = max(min(varphi_a, 1d0),0d0)

       S_temp = (1d0 - psi(ij_com+1))*mu_b*max(a_p, 1d-10)**egam/egam

       ! get value function
       do iw_p = 1, NW
         do ie_p = 1, NE
           EV = varphi_a*(egam*V(ij_com+1, ial, 0, iw_p, ie_p))**(1d0/egam) + &
                    (1d0-varphi_a)*(egam*V(ij_com+1, iar, 0, iw_p, ie_p))**(1d0/egam)
           S_temp = S_temp + pi_eta(iw, iw_p)*pi_theta(ie, ie_p)*psi(ij_com+1)*EV**egam/egam
         enddo
       enddo

       S(ij, ix_p, ik, iw, ie, 0) = S_temp

    end subroutine


  ! solve the household's consumption-savings decision
  subroutine solve_consumption_w(ij, ia, ik, iw, ie)

      implicit none

      integer, intent(in) :: ij, ia, ik, iw, ie
      real*8 :: x_in, fret

      ! set up communication variables
      ij_com = ij; ia_com = ia; ik_com = ik; iw_com = iw; ie_com = ie

      ! get best initial guess from future period
      x_in = max(X_plus_t(ij+1, ia, ik, iw, ie, 0), 1d-4)

      ! solve the household problem using rootfinding
      call fminsearch(x_in, fret, 0d0, a_u, cons_w)

      ! copy decisions
      X_plus_t(ij, ia, ik, iw, ie, 0) = x_in
      a_plus_t(ij, ia, ik, iw, ie, 0) = x_in
      k_plus_t(ij, ia, ik, iw, ie, 0) = 0d0
      c_t(ij, ia, ik, iw, ie, 0) = cons_com
      V_t(ij, ia, ik, iw, ie, 0) = -fret

  end subroutine


end module
