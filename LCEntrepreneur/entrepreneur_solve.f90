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
    subroutine solve_entrepreneur(ij, ix_p, ip_p, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, ix_p, ip_p, ik, iw, ie
        real*8 :: x_in, fret

        ! set up communication variables
        ij_com = ij; ix_p_com = ix_p; ip_p_com = ip_p; ik_com = ik; iw_com = iw; ie_com = ie

        if (X(ix_p) > (1d0-xi)*k_min + tr(k(ik), k_min)) then

           ! get best guess for the root of foc_real
           if(ix_p > 0)then
              x_in = omega_k(ij, ix_p-1, ip_p, ik, iw, ie)
           else
              x_in = 1d-4
           endif

           ! solve the household problem using fminsearch
           call fminsearch(x_in, fret, 0d0, 1d0, real_o)

           ! portfolio share for capital
           omega_k(ij, ix_p, ip_p, ik, iw, ie) = x_in
           S(ij, ix_p, ip_p, ik, iw, ie, 1) = -fret

        else

          omega_k(ij, ix_p, ip_p, ik, iw, ie) = 0d0
          S(ij, ix_p, ip_p, ik, iw, ie, 1) = 1d-10**egam/egam

        endif

    end subroutine


    ! solve the household's consumption-savings decision
    subroutine solve_consumption_e(ij, ia, ip, ik, iw, ie)

        implicit none

        integer, intent(in) :: ij, ia, ip, ik, iw, ie
        real*8 :: x_in(2), fret, varphi_x, varphi_p, k_p
        integer :: ixl_p, ixr_p, ipl_p, ipr_p

        ! set up communication variables
        ij_com = ij; ia_com = ia; ip_com = ip; ik_com = ik; iw_com = iw; ie_com = ie; io_p_com = 1

        ! get best initial guess from future period
        x_in(1) = max(X_plus_t(ij+1, ia, ip, ik, iw, ie, 1), 1d-4)
        x_in(2) = max(l_t(ij+1, ia, ip, ik, iw, ie, 1), 0.33d0)

        ! solve the household problem using rootfinding
        call fminsearch(x_in, fret, (/X_l, 0d0/), (/X_u, 0.8d0/), cons_w)

        call linint_Grow(x_in(1), x_l, x_u, x_grow, NX, ixl_p, ixr_p, varphi_x)
        call linint_Equi(p_plus_com, p_l, p_u, NP, ipl_p, ipr_p, varphi_p)

        ! restrict values to grid just in case
        ixl_p = min(ixl_p, NX)
        ixr_p = min(ixr_p, NX)
        varphi_x = max(min(varphi_x, 1d0),0d0)

        ipl_p = min(ipl_p, NP)
        ipr_p = min(ipr_p, NP)
        varphi_p = max(min(varphi_p, 1d0),0d0)

        ! get next period's capital size
        k_p = ((1d0-xi)*k_min + (varphi_x*varphi_p            *omega_k(ij, ixl_p, ipl_p, ik, iw, ie) +  &
                                 varphi_x*(1d0-varphi_p)      *omega_k(ij, ixl_p, ipr_p, ik, iw, ie) +  &
                                 (1d0-varphi_x)*varphi_p      *omega_k(ij, ixr_p, ipl_p, ik, iw, ie) +  &
                                 (1d0-varphi_x)*(1d0-varphi_p)*omega_k(ij, ixr_p, ipr_p, ik, iw, ie))*(x_in(1)-(1d0-xi)*k_min))/(1d0-xi)

        X_plus_t(ij, ia, ip, ik, iw, ie, 1) = x_in(1)
        a_plus_t(ij, ia, ip, ik, iw, ie, 1) = x_in(1) - (1d0-xi)*k_p
        k_plus_t(ij, ia, ip, ik, iw, ie, 1) = k_p
        c_t(ij, ia, ip, ik, iw, ie, 1) = cons_com
        l_t(ij, ia, ip, ik, iw, ie, 1) = x_in(2)
        V_t(ij, ia, ip, ik, iw, ie, 1) = -fret

    end subroutine

end module
