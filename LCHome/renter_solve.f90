!##############################################################################
! MODULE renter_solve
!
! copyright: Hans Fehr and Maurice Hofmann
!            University of Wuerzburg
!##############################################################################

module renter_solve

    use toolbox
    use globals
    use renter_foc

    implicit none

contains

  ! solve the household's portfolio decision
  subroutine solve_portfolio_r(ij, ic, il_p, id_p)

      implicit none

      integer, intent(in) :: ij, ic, il_p, id_p
      real*8 :: x_in, fret


      if (id_p /= 0) write(*,*)'Fehler!!'  

      ! set up communication variables
      ij_com = ij; ic_com = ic; il_p_com = il_p

      if (il_p == 0) then 
        
          omega_plus(ij, ic, 0, id_p) = omega_plus(ij, ic, 1, id_p)
          Q(ij, ic, il_p, id_p) = -liq_port_r(omega_plus(ij, ic, 0, id_p))
      
      else

         ! get best guess for the root of foc_port
         x_in = omega_plus(ij+1, ic, il_p, id_p)

         ! solve the household problem using rootfinding
         call fminsearch(x_in, fret, 0d0, 1d0, liq_port_r)

         omega_plus(ij, ic, il_p, id_p) = x_in
         Q(ij, ic, il_p, 0) = -fret

     endif    

  end subroutine


    ! solve the household's decision of how much wealth to store in real estate
    subroutine solve_renter(ij, ic, ia_p, ih)

        implicit none

        integer, intent(in) :: ij, ic, ia_p, ih
        integer :: ill, ilr   
        real*8 :: al_p, varphi_a


        al_p  = a(ia_p) - tr(h(ih),0d0)
 
        if (al_p >= 0d0) then
           
           ! calculate linear interpolation for future assets
           call linint_Grow(al_p, l_l, l_u, l_grow, NL, ill, ilr, varphi_a)

           ! restrict values to grid just in case
           ill = min(ill, NL)
           ilr = min(ilr, NL)
           varphi_a = max(min(varphi_a, 1d0),0d0)
        
           ! get value function
           S(ij, ic, ia_p, ih, 0) = max(varphi_a*(egam*Q(ij, ic, ill, 0))**(1d0/egam)  + &
                                     (1d0-varphi_a)*(egam*Q(ij, ic, ilr, 0))**(1d0/egam), 1d-10)**egam/egam

        else 

           S(ij, ic, ia_p, ih, 0) =  1d-10**egam/egam
          
        endif


    end subroutine


  ! solve the household's consumption-savings decision
  subroutine solve_consumption_r(ij, ix, ic, ih)

      implicit none

      integer, intent(in) :: ij, ix, ic, ih
      real*8 :: x_in, fret, amin

      ! set up communication variables
      ij_com = ij; ix_com = ix; ic_com = ic; ih_com = ih

      amin = tr(h(ih), 0d0)

      if(X(ix) < amin)then
          x_in = 0d0
          cons_com = 0d0
          ch_com = 0d0
          fret = -1d-10**egam/egam
      else

          ! get best initial guess from future period
          x_in = max(a_plus_t(ij+1, ix, ic, ih, 0), amin)

          ! solve the household problem using rootfinding
          call fminsearch(x_in, fret, amin, a_u, cons_r)
      endif

      ! copy decisions
      a_plus_t(ij, ix, ic, ih, 0) = x_in
      h_plus_t(ij, ix, ic, ih, 0) = 0d0
      c_t(ij, ix, ic, ih, 0) = cons_com
      ch_t(ij, ix, ic, ih, 0) = ch_com
      V_t(ij, ix, ic, ih, 0) = -fret

  end subroutine
  

end module
