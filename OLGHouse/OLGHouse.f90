!##############################################################################
! PROGRAM OLG-House
!
! copyright: Hans Fehr & Maurice Hofmann
!            University of Wuerzburg
!            maurice.hofmann@uni-wuerzburg.de
!##############################################################################
include'toolbox.f90'
include'globals.f90'

program OLG_House

    ! modules
    use globals

    implicit none

    ! initialize remaining variables
    call initialize()

    ! calculate initial equilibrium
    call get_SteadyState()

    call toc
!    call plotoutput()

    call output()
    close(21)

contains


! computes the initial steady state of the economy
subroutine get_SteadyState()

    implicit none

    ! start timer
    call tic()

    ! iterate until value function converges
    do iter = 1, itermax

        ! get new prices
        call prices()

        ! solve the household problem
        call solve_household()

        ! calculate the distribution of households over state space
        call get_distribution()

        ! aggregate individual decisions
        call aggregation()

        ! determine the government parameters
        call government()

        write(*,'(i4,2i7,6f8.2,f12.5)')iter, maxval(iamax), maxval(ihmax),&
                                        (/5d0*KK, CC, II, IHO/)/YY*100d0, &
                                        ((1d0+r)**0.2d0-1d0)*100d0, w, DIFF/YY*100d0

        if(abs(DIFF/YY)*100d0 < sig) return
    enddo

    write(*,*)'No Convergence'


end subroutine


! initializes the remaining model parameters and variables
subroutine initialize

    implicit none

    integer :: ij

    write(*,'(/a/)')'INITIAL EQUILIBRIUM'
    write(*,'(a)')'ITER  IAMAX  IHMAX     K/Y     C/Y    II/Y       r       w        DIFF'

    ! initialize survival probabilities
    psi(:, 1)=(/   1.00000000d0,&
                   0.99834294d0,&
                   0.99809443d0,&
                   0.99745495d0,&
                   0.99635413d0,&
                   0.99410298d0,&
                   0.98972953d0,&
                   0.98185396d0,&
                   0.97070373d0,&
                   0.95530594d0,&
                   0.93417914d0,&
                   0.90238714d0,&
                   0.83653436d0,&
                   0.71048182d0,&
                   0.52669353d0,&
                   0.31179803d0,&
                   0.00000001d0/)

    psi(:, 2) = psi(:, 1)
    psi(:, 3) = psi(:, 1)

    ! initialize age earnings process
    eff(1,:) = (/ 1.2987778,   1.4327164,   1.3882564  /)
    eff(2,:) = (/ 1.5794954,   1.8210024,   2.1841104  /)
    eff(3,:) = (/ 1.6404434,   1.9747812,   2.9655702  /)
    eff(4,:) = (/ 1.6908550,   2.0647004,   3.3290738  /)
    eff(5,:) = (/ 1.7507724,   2.1559744,   3.4171474  /)
    eff(6,:) = (/ 1.7586790,   2.2020510,   3.4497238  /)
    eff(7,:) = (/ 1.7611338,   2.2484878,   3.4046532  /)
    eff(8,:) = (/ 1.8054554,   2.2359332,   3.3062074  /)
    eff(9,:) = (/ 1.7423268,   2.1737906,   3.1235630  /)
    eff(JR:JJ,:) = 0d0

    ! distribution of skill classes
    dist_skill = (/0.263d0, 0.545d0, 0.192d0/)

    call discretize_AR(0.95666d0**5d0, 0.0d0, sigma5(0.95666d0, 0.02321d0), eta(1, :), pi(1, :, :), dist_eta(1, :))
    eta(1, :) = exp(eta(1, :))/sum(dist_eta(1,:)*exp(eta(1, :)))

    call discretize_AR(0.95687d0**5d0, 0.0d0, sigma5(0.95687d0, 0.02812d0), eta(2, :), pi(2, :, :), dist_eta(2, :))
    eta(2, :) = exp(eta(2, :))/sum(dist_eta(2,:)*exp(eta(2, :)))

    call discretize_AR(0.95828d0**5d0, 0.0d0, sigma5(0.95828d0, 0.03538d0), eta(3, :), pi(3, :, :), dist_eta(3, :))
    eta(3, :) = exp(eta(3, :))/sum(dist_eta(3,:)*exp(eta(3, :)))


    ! set up population structure
    rpop(0) = 1d0+n_p
    do ij = 1, JJ+1
        rpop(ij) = rpop(ij-1)*psi(ij, 1)/(1d0+n_p)
    enddo

    workpop = 0d0
    ! initialize workforce
    do ij = 1, JJ
       if(ij < JR)workpop = workpop + rpop(ij)
    enddo

    ! initialize asset grid
    call grid_Cons_Grow(a, a_l, a_u, a_grow)

    ! initialize liquid asset grid
    call grid_Cons_Grow(al, al_l, al_u, al_grow)

    ! initialize earning points grid
    call grid_Cons_Grow(ep, ep_l, ep_u, ep_grow)

    ! initialize housing grid
    call grid_Cons_Grow(h(1:NH), h_l, h_u, h_grow)
    h(0) = 0d0


    ! specify which should be the endogenous budget-balancing tax rate
    tax   = 1

    ! initialize tax rates
    tauc  = 0.0d0
    taup  = 0.0d0
    taur  = 0.0d0
    tauw  = 0.0d0

    ! specify the extent of the public good and government debt
    gy    = 0.0d0
    by    = 0.0d0
    by    = by/5d0

    ! set starting values
    KK = 10d0
    LL = 10d0
    BQ = 3d0

    ! initial guess bequests
    b = 0d0
    do ij = 1, JR-1
       b(ij) = BQ/workpop
    enddo

    ! initial guess average income
    ybar = 2d0

    ! initialize value functions
    SO = 1d-10**egam/egam; SR = 1d-10**egam/egam;  V = 1d-10**egam/egam

    ! initialize policy functions
    aplus = 0d0; hplus = 0; c = 0d0; ch = 0d0
    omegaplus = 0d0

    ! initialize temporary policy and value functions
    aplus_t = 0d0; hplus_t = 0d0; ch_t = 0d0; c_t = 0d0; V_t = 0d0


    ! open files
    open(21, file='output.out')
!    open(22, file='output_owner.out')
!    open(23, file='output_renter.out')

end subroutine


! compute prices and distribution of bequests for next iteration step
subroutine prices()

    implicit none

    integer :: ij

    ! calculate new prices
    r = Omega*alpha*(KK/LL)**(alpha-1d0)-delta_k
    w = Omega*(1d0-alpha)*(KK/LL)**alpha
    rn = r*(1d0 - taur)
    wn = w*(1d0 - tauw - taup)
    p  = 1d0 + tauc
    ph = r + delta_h


    ! compute bequest per capita within workforce for next iteration step
    do ij = 1, JR-1
       b(ij) = BQ/workpop
    enddo

    !write(*,*)r,w,rn,wn,p,ph

end subroutine


! determines the solution to the household optimization problem
subroutine solve_household()

    implicit none

    integer :: ij, is, ia, ie, ih, ip, ia_p, is_max, ip_max
    real*8 :: al_p

    call tic

    ! solve household problem recursively
    do ij = JJ, 1, -1

        ! check about how many is to iterate
        if(ij >= JR)then
            is_max = 1
            ip_max = 1
        else
            is_max = NS
            ip_max = NP
        endif

        if(ij == JJ)then

            omegaplus(JJ, :, :, :, :, :, :) = 0d0
            hplus(JJ, :, :, :, :, :) = 0d0
            aplus_t(JJ+1, :, : , :, :, :, :) = 0d0
            hplus_t(JJ+1, :, : , :, :, :, :) = 0d0

            do ih = 0, NH
                do ia_p = 0, NA
                    al_p = a(ia_p)-tr(h(ih),0d0)
                    if (al_p > 0d0) then
                         SO(JJ, :, :, :, ih, :, ia_p) = q_1*(1d0+al_p/q_2)**egam
                         SR(JJ, :, :, :, ih, :, ia_p) = q_1*(1d0+al_p/q_2)**egam
                    else
                        SO(JJ, :, :, :, ih, :, ia_p) = 1d-10**egam/egam
                        SR(JJ, :, :, :, ih, :, ia_p) = 1d-10**egam/egam
                    endif
                enddo
            enddo

            do ia = 0, NA
                do ie = 0, NE
                    do ih = 0, NH

                        if (q_1 == 0d0) then

                            aplus(JJ, :, ia, ie, ih, :) = 0d0
                            y(JJ, :, ia, ie, ih, : ) = 0d0
                            penp(JJ, :, ia, ie, ih, : ) = kappa*ybar*ep(ie)
                            penc(JJ, :, ia, ie, ih, : ) = 0d0
                            yg(JJ, :, ia, ie, ih, : ) = (1d0+ rn)*(al(ia)-xi*h(ih)) + (1d0-delta_h)*h(ih) + b(ij) + &
                                                        penp(ij, 1, ia, ie, ih, 1)

                            if (ih == 0) then
                                ch_com = max((1d0-theta)*yg(JJ, 1, ia, ie, ih, 1)/ph, 1d-10)
                                cons_com = max((yg(JJ, 1, ia, ie, ih, 1)-ph*ch_com)/p, 1d-10)
                                V(JJ, :, ia, ie, ih, :) = (cons_com**theta*(chi*ch_com)**(1d0-theta))**egam/egam
                                c(JJ, :, ia, ie, ih, :) = cons_com
                                ch(JJ, :, ia, ie, ih, :) = ch_com
                            else
                                ch_com = h(ih)
                                cons_com = max(yg(JJ, 1, ia, ie, ih, 1)/p, 1d-10)
                                V(JJ, :, ia, ie, ih, :) = (cons_com**theta*ch_com**(1d0-theta))**egam/egam
                                c(JJ, :, ia, ie, ih, :) = cons_com
                                ch(JJ, :, ia, ie, ih, :) = ch_com
                            endif

                        else

                            ! with bequest motive we assume future renter
                            call solve_consumption_r(JJ, 1, ia, ie, ih, 1)

                            V(JJ, :, ia, ie, ih, :) = V_t(JJ, 1, ia, ie, ih, 1, 0)
                            aplus(JJ, :, ia, ie, ih, :) = aplus_t(JJ, 1, ia, ie, ih, 1, 0)
                            c(JJ, :, ia, ie, ih, :) = c_t(JJ, 1, ia, ie, ih, 1, 0)
                            ch(JJ, :, ia, ie, ih, :) = ch_t(JJ, 1, ia, ie, ih, 1, 0)
                        endif
                    enddo
                enddo
            enddo
        else

            ! get optimal share of wealth stored in real estate
            do is = 1, is_max
                do ia = 0, NA
                    do ie = 0, NE
                        do ih = 0, NH
                            do ip = 1, ip_max

                	            ! income and pension system
                                if (ij >= JR)then
                                     y(ij, :, ia, ie, ih, :) = 0d0
                                     penc(ij, :, ia, ie, ih, :) = 0d0
                                     penp(ij, :, ia, ie, ih, :) = kappa*ybar*ep(ie)
                                     yg(JJ, :, ia, ie, ih, : ) = (1d0 + rn)*(al(ia)-xi*h(ih)) + (1d0-delta_h)*h(ih) + b(ij) + &
                                                                  penp(ij, 1, ia, ie, ih, 1)
                                else
                                     y(ij, is, ia, ie, ih, ip) = w*eff(ij, is)*eta(is, ip)
                                     penc(ij, is, ia, ie, ih, ip) = taup*min(y(ij, is, ia, ie, ih, ip), 2d0*ybar)
                                     penp(ij, is, ia, ie, ih, ip) = 0d0
                                     yg(JJ, :, ia, ie, ih, : ) = (1d0 + rn)*(al(ia)-xi*h(ih)) + (1d0-delta_h)*h(ih) + b(ij) + &
                                                                  w*eff(ij, is)*eta(is, ip) - penc(ij, is, ia, ie, ih, ip)
                                endif

                                do ia_p = 0, NA

                                    ! next period homeowner
                                    call solve_owner(ij, is, ia, ie, ih, ip, ia_p)

                                    ! next period renter
                                    call solve_renter(ij, is, ia, ie, ih, ip, ia_p)
                                enddo

                            enddo
                        enddo
                    enddo
                enddo
            enddo

            ! solve the consumption savings problem
            do is = 1, is_max
                do ia = 0, NA
                    do ie = 0, NE
                        do ih = 0, NH
                            do ip = 1, ip_max


                                ! next period homeowner
                                call solve_consumption_o(ij, is, ia, ie, ih, ip)

                                ! next period renter
                                call solve_consumption_r(ij, is, ia, ie, ih, ip)
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            ! decision whether to be owner or renter next period
            do is = 1, is_max
                do ia = 0, NA
                    do ie = 0, NE
                        do ih = 0, NH
                            do ip = 1, ip_max

                                ! decision on whether to be homeowner or renter next period
                                if (V_t(ij, is, ia, ie, ih, ip, 1) >= V_t(ij, is, ia, ie, ih, ip, 0) ) then
                                    hplus(ij, is, ia, ie, ih, ip) = hplus_t(ij, is, ia, ie, ih, ip, 1)
                                    aplus(ij, is, ia, ie, ih, ip) = aplus_t(ij, is, ia, ie, ih, ip, 1)
                                    c(ij, is, ia, ie, ih, ip) = c_t(ij, is, ia, ie, ih, ip, 1)
                                    ch(ij, is, ia, ie, ih, ip) = ch_t(ij, is, ia, ie, ih, ip, 1)
                                    V(ij, is, ia, ie, ih, ip) = V_t(ij, is, ia, ie, ih, ip, 1)
                                else
                                    hplus(ij, is, ia, ie, ih, ip) = hplus_t(ij, is, ia, ie, ih, ip, 0)
                                    aplus(ij, is, ia, ie, ih, ip) = aplus_t(ij, is, ia, ie, ih, ip, 0)
                                    c(ij, is, ia, ie, ih, ip) = c_t(ij, is, ia, ie, ih, ip, 0)
                                    ch(ij, is, ia, ie, ih, ip) = ch_t(ij, is, ia, ie, ih, ip, 0)
                                    V(ij, is, ia, ie, ih, ip) = V_t(ij, is, ia, ie, ih, ip, 0)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo

        endif
        write(*,'(a,i3,a)')'Age: ',ij,' DONE!'

    enddo

    call toc


end subroutine



! determines the invariant distribution of households
subroutine get_distribution()

   implicit none

   integer :: ij, is, ia, ie, ih, ip, ip_p, ial, iar, iel, ier, ihl, ihr
   real*8 :: h_p, ep_p, varphi_a, varphi_h, varphi_e

   ! set distribution to zero
   phi(:, :, :, :, :, :) = 0d0

   ! get initial distribution at age 1
   do is = 1, NS
      do ip = 1, NP
         phi(1, is, 0, 0, 0, ip) = dist_skill(is)*dist_eta(is, ip)
      enddo
   enddo

   ! successively compute distribution over ages
   do ij = 1, JJ

       ! iterate over yesterdays gridpoints
       do is = 1, NS
           do ia = 0, NA
               do ie = 0, NE
                   do ih = 0, NH
                       do ip = 1, NP

                           h_p = hplus(ij, is, ia, ie, ih, ip)

                           ! interpolate yesterday's savings decision
                           call linint_Grow(aplus(ij, is, ia, ie, ih, ip), a_l, a_u, a_grow, NA, ial, iar, varphi_a)

                           ial = min(ial, NA)
                           iar = min(iar, NA)
                           varphi_a = max(min(varphi_a, 1d0), 0d0)

                           if (h_p>= h_min) then
                              call linint_Grow(h_p, h_l, h_u, h_grow, NH-1, ihl, ihr, varphi_h)

                              ihl = min(ihl+1, NH)
                              ihr = min(ihr+1, NH)
                              varphi_h = max(min(varphi_h, 1d0), 0d0)

                              ! redistribute households in period ij+1 based on their decision in ij and shock ip_p
                              if (ij+1>=JR) then

                                  phi(ij+1, is, ial, ie, ihl, ip) = phi(ij+1, is, ial, ie, ihl, ip) + &
                                                                      varphi_a*varphi_h*phi(ij, is, ia, ie, ih, ip)
                                  phi(ij+1, is, ial, ie, ihr, ip) = phi(ij+1, is, ial, ie, ihr, ip) + &
                                                                      varphi_a*(1d0-varphi_h)*phi(ij, is, ia, ie, ih, ip)
                                  phi(ij+1, is, iar, ie, ihl, ip) = phi(ij+1, is, iar, ie, ihl, ip) + &
                                                                      (1d0-varphi_a)*varphi_h*phi(ij, is, ia, ie, ih, ip)
                                  phi(ij+1, is, iar, ie, ihr, ip) = phi(ij+1, is, iar, ie, ihr, ip) + &
                                                                      (1d0-varphi_a)*(1d0-varphi_h)*phi(ij, is, ia, ie, ih, ip)

                              else

                                  ! get tomorrow's earning points
                                  ep_p = ep(ie) + mu*(lambda+(1d0-lambda)*min(y(ij, is, ia, ie, ih, ip)/ybar, 2d0))/dble(JR-1d0)

                                  ! derive interpolation weights
                                  call linint_Grow(ep_p, ep_l, ep_u, ep_grow, NE, iel, ier, varphi_e)

                                  iel = min(iel, NE)
                                  ier = min(ier, NE)
                                  varphi_e = max(min(varphi_e, 1d0), 0d0)

                                  do ip_p = 1, NP

                                      phi(ij+1, is, ial, iel, ihl, ip_p) = phi(ij+1, is, ial, iel, ihl, ip_p) + pi(is, ip, ip_p)&
                                                             *varphi_a*varphi_e*varphi_h*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, ial, iel, ihr, ip_p) = phi(ij+1, is, ial, iel, ihr, ip_p) + pi(is, ip, ip_p)&
                                                             *varphi_a*varphi_e*(1d0-varphi_h)*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, ial, ier, ihl, ip_p) = phi(ij+1, is, ial, ier, ihl, ip_p) + pi(is, ip, ip_p)&
                                                             *varphi_a*(1d0-varphi_e)*varphi_h*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, ial, ier, ihr, ip_p) = phi(ij+1, is, ial, ier, ihr, ip_p) + pi(is, ip, ip_p)&
                                                             *varphi_a*(1d0-varphi_e)*(1d0-varphi_h)*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, iar, iel, ihl, ip_p) = phi(ij+1, is, iar, iel, ihl, ip_p) + pi(is, ip, ip_p)&
                                                             *(1d0-varphi_a)*varphi_e*varphi_h*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, iar, iel, ihr, ip_p) = phi(ij+1, is, iar, iel, ihr, ip_p) + pi(is, ip, ip_p)&
                                                             *(1d0-varphi_a)*varphi_e*(1d0-varphi_h)*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, iar, ier, ihl, ip_p) = phi(ij+1, is, iar, ier, ihl, ip_p) + pi(is, ip, ip_p)&
                                                             *(1d0-varphi_a)*(1d0-varphi_e)*varphi_h*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, iar, ier, ihr, ip_p) = phi(ij+1, is, iar, ier, ihr, ip_p) + pi(is, ip, ip_p)&
                                                        *(1d0-varphi_a)*(1d0-varphi_e)*(1d0-varphi_h)*phi(ij, is, ia, ie, ih, ip)

                                  enddo
                              endif

                           else

                              if (ij+1>=JR) then

                                  phi(ij+1, is, ial, ie, 0, ip) = phi(ij+1, is, ial, ie, 0, ip) + &
                                                                      varphi_a*phi(ij, is, ia, ie, ih, ip)
                                  phi(ij+1, is, iar, ie, 0, ip) = phi(ij+1, is, iar, ie, 0, ip) + &
                                                                      (1d0-varphi_a)*phi(ij, is, ia, ie, ih, ip)

                              else

                                  ! get tomorrow's earning points
                                  ep_p = ep(ie) + mu*(lambda+(1d0-lambda)*min(y(ij, is, ia, ie, ih, ip)/ybar, 2d0))/dble(JR-1d0)

                                  ! derive interpolation weights
                                  call linint_Grow(ep_p, ep_l, ep_u, ep_grow, NE, iel, ier, varphi_e)

                                  iel = min(iel, NE)
                                  ier = min(ier, NE)
                                  varphi_e = max(min(varphi_e, 1d0), 0d0)

                                  do ip_p = 1, NP

                                      phi(ij+1, is, ial, iel, 0, ip_p) = phi(ij+1, is, ial, iel, 0, ip_p) + pi(is, ip, ip_p)&
                                                             *varphi_a*varphi_e*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, ial, ier, 0, ip_p) = phi(ij+1, is, ial, iel, 0, ip_p) + pi(is, ip, ip_p)&
                                                             *varphi_a*(1d0-varphi_e)*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, iar, iel, 0, ip_p) = phi(ij+1, is, iar, iel, 0, ip_p) + pi(is, ip, ip_p)&
                                                             *(1d0-varphi_a)*varphi_e*phi(ij, is, ia, ie, ih, ip)
                                      phi(ij+1, is, iar, ier, 0, ip_p) = phi(ij+1, is, iar, ier, 0, ip_p) + pi(is, ip, ip_p)&
                                                             *(1d0-varphi_a)*(1d0-varphi_e)*phi(ij, is, ia, ie, ih, ip)
                                  enddo

                              endif

                           endif

                       enddo
                   enddo
               enddo
           enddo
       enddo

       if (sum(phi(ij, :, :, :, :, :)) < 1d0-1d-10)write(*,'(a, i5, f8.2)')'X: ', ij, sum(phi(ij, :, :, :, :, :))

   enddo

   ! check maximum grid points used
   call check_grid(iamax, ihmax)

end subroutine


! subroutine for calculating macroecomic aggregates
subroutine aggregation()

    implicit none

    integer :: ij, is, ia, ie, ih, ip
    real*8 ::  b_old(JJ)

    ! copy bequests from previous iteration for damping
    b_old = b

    ! reset macroeconomic aggregates in each iteration step
    AA = 0d0; BQ = 0d0; CC = 0d0; AAL = 0d0
    HH = 0d0; HR = 0d0; LL = 0d0

    a_coh = 0d0; c_coh = 0d0; h_coh = 0d0; al_coh = 0d0; ch_coh = 0d0
    penp_coh = 0d0; penc_coh = 0d0; v_coh = 0d0; l_coh = 0d0; tr_coh = 0d0
    shr = 0d0

    ! compute macroeconomic aggregates (normalized to the youngest cohort alive in period t)
    do ij = 1, JJ+1
        do is = 1, NS
            do ia = 0, NA
                do ie = 0, NE
                    do ih = 0, NH
                        do ip = 1, NP

                            a_coh(ij) = a_coh(ij) + a(ia)*phi(ij, is, ia, ie, ih, ip)/psi(ij, is)
                            c_coh(ij) = c_coh(ij) + c(ij, is, ia, ie, ih, ip)*phi(ij, is, ia, ie, ih, ip)
                            h_coh(ij) = h_coh(ij) + h(ih)*phi(ij, is, ia, ie, ih, ip)/psi(ij, is)
                            l_coh(ij) = l_coh(ij) + eff(ij, is)*eta(is, ip)*phi(ij, is, ia, ie, ih, ip)
                            al_coh(ij) = al_coh(ij) + al(ia)*phi(ij, is, ia, ie, ih, ip)/psi(ij, is)
                            ch_coh(ij) = ch_coh(ij) + ch(ij, is, ia, ie, ih, ip)*phi(ij, is, ia, ie, ih, ip)
                            tr_coh(ij+1) = tr_coh(ij+1) + tr(h(ih), hplus(ij, is, ia, ie, ih, ip))/psi(ij, is)
                            penp_coh(ij) = penp_coh(ij) + penp(ij, is, ia, ie, ih, ip)*phi(ij, is, ia, ie, ih, ip)
                            penc_coh(ij) = penc_coh(ij) + penc(ij, is, ia, ie, ih, ip)*phi(ij, is, ia, ie, ih, ip)
                            v_coh(ij) = v_coh(ij) + V(ij, is, ia, ie, ih, ip)*phi(ij, is, ia, ie, ih, ip)
                            if (ih == 0) then
                               shr(ij, 0) = shr(ij, 0) + phi(ij, is, ia, ie, ih, ip)
                            else
                               shr(ij, 1) = shr(ij, 1) + phi(ij, is, ia, ie, ih, ip)
                            endif

                        enddo
                    enddo
                enddo
            enddo
       enddo
        AA = AA + a_coh(ij)*rpop(ij)
        CC = CC + c_coh(ij)*rpop(ij)
        HH = HH + h_coh(ij)*rpop(ij)
        LL = LL + l_coh(ij)*rpop(ij)
        AAL = AAL + al_coh(ij)*rpop(ij)
        TRG = TRG + tr_coh(ij)*rpop(ij)
        HR = HR + ch_coh(ij)*rpop(ij)
        PBEN = PBEN + penp_coh(ij)*rpop(ij)
        PCON = PCON + penc_coh(ij)*rpop(ij)
        BQ = BQ + (1d0-psi(ij, 1))*((1d0+rn)*a_coh(ij)+(1d0-delta_h)*h_coh(ij))*rpop(ij)/psi(ij, 1)

    enddo

    ! get average income
    ybar = w*LL/workpop

    ! compute stock of capital
    KK = damp*(AAL-xi*HH-BB)+(1d0-damp)*KK

    !write(*,*)KK, LL, BQ, ybar

    ! labor supply, investions and total output
    II = (n_p+delta_k)*KK
    IHO = (n_p+delta_h)*HH
    IHR = (n_p + delta_h)*HR
    YY = Omega*KK**alpha*LL**(1d0-alpha)

    ! compute gap on goods market
    DIFF = YY-CC-GG-II-IHO-IHR-TRG

end subroutine


! subroutine for calculating government parameters
subroutine government()

    implicit none

    real*8 :: expend, tauc_old, tauh_old, taup_old, taur_old, tauw_old

    ! copy tax rate from previous iteration for damping
    tauc_old = tauc; taup_old = taup
    taur_old = taur; tauw_old = tauw

    ! set government spending and debt
    BB = by*YY
    GG = gy*YY

    ! calculate total government expenditure
    expend = GG + (1d0+r)*BB - (1d0+n_p)*BB

    ! get budget balancing tax rate
    if(tax == 1)then
        tauc = (expend - (tauw*w*LL + taur*r*(AAL-xi*HH)))/CC
        tauc = damp*tauc + (1d0-damp)*tauc_old
    elseif(tax == 2)then
        taur = (expend - (tauc*CC + tauw*w*LL))/(r*(AAL-xi*HH))
        taur = damp*taur + (1d0-damp)*taur_old
    else
        tauw = (expend - (tauc*CC + taur*r*(AAL-xi*HH)))/(w*LL)
        tauw = damp*tauw + (1d0-damp)*tauw_old
    endif

    ! compute total tax revenue
    taxrev(1) = tauc*CC; taxrev(2) = taur*r*(AAL-xi*HH)
    taxrev(3) = tauw*w*LL; taxrev(4) = sum(taxrev(1:3))

    ! obtain aggregated contribution basis
    PCON = max(PCON/taup,1d-10)

    ! get budget balancing pension contribution rate
    taup = PBEN/PCON

    ! damping pension contribution rate
    taup = damp*taup + (1d0-damp)*taup_old

    !write(*,*)BB, GG, tauc, taup

end subroutine

! subroutine that checks for the maximum asset gridpoint used
subroutine check_grid(iamax_local,ihmax_local)

    implicit none

    integer :: iamax_local(JJ+1), ihmax_local(JJ+1)
    integer :: ij, is, ia, ie, ih, ip

    iamax_local = 0; ihmax_local = 0

    do ij = 1, JJ+1

        ! check for the maximum grid points used at a certain age
        do is = 1, NS
            do ia = 0, NA
                do ie = 0, NE
                    do ih = 0, NH
                        do ip = 1, NP
                            if(phi(ij, is, ia, ie, ih, ip) > 1d-8)iamax_local(ij) = ia
                            if(phi(ij, is, ia, ie, ih, ip) > 1d-8)ihmax_local(ij) = ih
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

end subroutine


! subroutine for writing output
subroutine output()

    implicit none

    integer :: ij, is, ia, ie, ih, ip


    ! Output
    write(21,'(a/)')'STEADY STATE EQUILIBRIUM'
    write(21,'(a)')'CAPITAL        A      HH      HR      K       B      BQ       r    p.a.'
    write(21,'(8x,8f8.2)')AA, HH, HR, KK, BB,  BQ, r, ((1d0+r)**(1d0/5d0)-1d0)*100d0
    write(21,'(a,8f6.2/)')'(in %)  ',(/AA, HH, HR, KK, BB, BQ/)/YY*500d0

    write(21,'(a)')'LABOR          L    ybar       w      ph'
    write(21,'(8x, 4f8.2/)')LL, ybar, w, ph

    write(21,'(a)')'GOODS          Y       C       I     IH     HR      TR       G      DIFF'
    write(21,'(8x,7f8.2,f8.3)')YY,CC,II,IHO,IHR,TRG,GG,DIFF
    write(21,'(a,8f8.2/)')'(in %)  ',(/YY, CC, II, IHO, IHR, TRG, GG, DIFF/)/YY*100d0

!    write(21,'(a)')'GOV         TAUC    ITAX    TAUx    TAUH   TOTAL      HA      HS      WF   (r-n)B'
!    write(21,'(8x,9f8.2)')taxrev(1:5), GHA, GHS, WP,(r-n_p)*BB
!    write(21,'(a,4f8.2/)')'(rate)  ',(/tauc,  0.0, taur, tauh/)*100d0

!    write(21,'(a)')'PENS        TAUP     PC      PP'
!    write(21,'(a,3f8.2/)')'(in %)  ',(/taup, PC/YY, PP/YY/)*100d0

    write(21, '(a,a)')' IJ      CONS     ASSETS    OWNING   RENTING    TRCOST     HDEBT    LABOR  EARNINGS    PENB     PENC/P  ',&
       '   RSHR      OSHR    IAMAX     IHMAX'
    do ij = 1, JJ
        write(21,'(i3,12f10.3,2i10)')ij, c_coh(ij), a_coh(ij), h_coh(ij), ch_coh(ij), tr_coh(ij), al_coh(ij)-xi*h_coh(ij), &
              w*l_coh(ij),w*l_coh(ij)+r*(al_coh(ij)-xi*h_coh(ij)), penp_coh(ij), penc_coh(ij), shr(ij,0), shr(ij,1), &
              iamax(ij), ihmax(ij)
    enddo
    write(21,'(13x,2f10.3,70x,2f10.3/)')a_coh(JJ+1),h_coh(JJ+1),shr(JJ+1,0),shr(JJ+1,1)


    ! plot output
!    call plotoutput()

end subroutine


end program
