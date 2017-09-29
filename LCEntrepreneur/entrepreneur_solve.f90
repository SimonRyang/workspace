!##############################################################################
! MODULE owner_solve
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################

module entrepreneur_solve

    use toolbox
    use globals
    use entrepreneur_foc

    implicit none

contains

    ! solve the household's decision of how much wealth to invest into firm capital
    subroutine solve_entrepreneur(ij, ix_p, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, ix_p, ik, iw, ie
        real*8 :: x_in, fret

        ! set up communication variables
        ij_com = ij; ix_p_com = ix_p; ik_com = ik; iw_com = iw; ie_com = ie

        if (X(ix_p) > (1d0-xi)*k_min) then

           ! get best guess for the root of foc_real
           if(ix_p > 0)then
              x_in = omega_k(ij, ix_p-1, ik, iw, ie)
           else
              x_in = 0d0
           endif

           ! solve the household problem using fminsearch
           call fminsearch(x_in, fret, 0d0, 1d0, real_o)

           ! portfolio share for capital
           omega_k(ij, ix_p, ik, iw, ie) = x_in
           S(ij, ix_p, ik, iw, ie, 1) = -fret

        else

          omega_k(ij, ix_p, ik, iw, ie) = 0d0
          S(ij, ix_p, ik, iw, ie, 1) = 1d-10**egam/egam

        endif

    end subroutine



    ! solve the household's consumption-savings decision
    subroutine solve_consumption_e(ij, ia, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, ia, ik, iw, ie
        real*8 :: x_in, fret, varphi_x, k_p
        integer :: ixl, ixr

        ! set up communication variables
        ij_com = ij; ia_com = ia; ik_com = ik; iw_com = iw; ie_com = ie

        ! get best initial guess from future period
        x_in = max(X_plus_t(ij+1, ia, ik, iw, ie, 1), 1d-4)

        ! solve the household problem using rootfinding
        call fminsearch(x_in, fret, 0d0, X_u, cons_e)

        call linint_Grow(x_in, x_l, x_u, x_grow, NX, ixl, ixr, varphi_x)

        ! restrict values to grid just in case
        ixl = min(ixl, NX)
        ixr = min(ixr, NX)
        varphi_x = max(min(varphi_x, 1d0),0d0)

        ! get next period's capital size
        k_p = ((1d0-xi)*k_min + (varphi_x*omega_k(ij, ixl, ik, iw, ie) +  &
               (1d0-varphi_x)*omega_k(ij, ixr, ik, iw, ie))*(x_in-(1d0-xi)*k_min))/(1d0-xi)

        X_plus_t(ij, ia, ik, iw, ie, 1) = x_in
        a_plus_t(ij, ia, ik, iw, ie, 1) = x_in - (1d0-xi)*k_p
        k_plus_t(ij, ia, ik, iw, ie, 1) = k_p
        c_t(ij, ia, ik, iw, ie, 1) = cons_com
        V_t(ij, ia, ik, iw, ie, 1) = -fret

    end subroutine

end module
