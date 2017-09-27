!##############################################################################
! MODULE owner_solve
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################

module owner_solve

    use toolbox
    use globals
    use owner_foc

    implicit none

contains

    ! solve the household's portfolio decision
    subroutine solve_portfolio_o(ij, ic, il_p, id_p)

        implicit none

        integer, intent(in) :: ij, ic, il_p, id_p
        real*8 :: x_in, fret

        ! set up communication variables
        ij_com = ij; ic_com = ic; il_p_com = il_p; id_p_com = id_p

        if (il_p == 0) then 

          omega_plus(ij, ic, 0, id_p) = omega_plus(ij, ic, 1, id_p)
          Q(ij, ic, il_p, id_p) = -liq_port_o(omega_plus(ij, ic, 0, id_p))


        elseif(l(il_p) > xi*d(id_p)/(1d0-xi)) then

           ! get best guess for the root of foc_port
           x_in = omega_plus(ij+1, ic, il_p, id_p)

           ! solve the household problem using rootfinding
           call fminsearch(x_in, fret, 0d0, 1d0, liq_port_o)

           omega_plus(ij, ic, il_p, id_p) = x_in
           Q(ij, ic, il_p, id_p) = -fret

        else 

          omega_plus(ij, ic, il_p, id_p) = 0d0
          Q(ij, ic, il_p, id_p) = -liq_port_o(-1d0)

        endif  

    end subroutine


    ! solve the household's decision of how much wealth to store in real estate
    subroutine solve_realestate(ij, ic, ia_p, ih)

        implicit none

        integer, intent(in) :: ij, ic, ia_p, ih
        real*8 :: x_in, fret

        ! set up communication variables
        ij_com = ij; ic_com = ic; ia_p_com = ia_p; ih_com = ih

        if (a(ia_p) > ((1d0-xi)*h_min + tr(h(ih), h_min))) then 

           ! get best guess for the root of foc_real
           if(ia_p > 0)then
              x_in = omega_h(ij, ic, ia_p-1, ih)
           else
              x_in = 0d0
           endif

           ! solve the household problem using fminsearch
           call fminsearch(x_in, fret, 0d0, 1d0, real_o)

           ! portfolio share for housing
           omega_h(ij, ic, ia_p, ih) = x_in
           S(ij, ic, ia_p, ih, 1) = -fret
       else 

           ! portfolio share for housing
           omega_h(ij, ic, ia_p, ih) = 0d0
           S(ij, ic, ia_p, ih, 1) = 1d-10**egam/egam

       endif  

    end subroutine



    ! solve the household's consumption-savings decision
    subroutine solve_consumption_o(ij, ix, ic, ih)

        implicit none

        integer, intent(in) :: ij, ix, ic, ih
        real*8 :: x_in, fret, amin, varphi_a, h_p
        integer :: ial, iar

        ! set up communication variables
        ij_com = ij;ix_com = ix; ic_com = ic; ih_com = ih

        amin = (1d0-xi)*h_min - tr(h(ih_com), h_min)

        if(X(ix) < amin)then
            x_in = 0d0
            h_p = 0d0
            cons_com = 0d0
            ch_com = 0d0
            fret = -1d-10**egam/egam
        else

            ! get best initial guess from future period
            x_in = max(a_plus_t(ij+1, ix, ic, ih, 1), amin)

            ! solve the household problem using rootfinding
            call fminsearch(x_in, fret, amin, a_u, cons_o)

            call linint_Grow(x_in, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

            ! restrict values to grid just in case
            ial = min(ial, NA)
            iar = min(iar, NA)
            varphi_a = max(min(varphi_a, 1d0),0d0)

            ! get next period's house size
            h_p = max((varphi_a*omega_h(ij, ic, ial, ih) +  &
                 (1d0-varphi_a)*omega_h(ij, ic, iar, ih))*x_in, 1d-10)/(1d0-xi)

        endif


        a_plus_t(ij, ix, ic, ih, 1) = x_in
        h_plus_t(ij, ix, ic, ih, 1) = h_p
        c_t(ij, ix, ic, ih, 1) = cons_com
        ch_t(ij, ix, ic, ih, 1) = ch_com
        V_t(ij, ix, ic, ih, 1) = -fret

    end subroutine

end module
