module clock
    
    implicit none
    
    integer :: calc
    
contains

    subroutine tick(t)
    
        integer, intent(OUT) :: t

        call system_clock(t)
        
    end subroutine tick

    ! returns time in seconds from now to time described by t
    subroutine tock(t)
    
        integer, intent(in) :: t
        real*8 :: calctime
        integer :: now, clock_rate
        
        call system_clock(now,clock_rate)

        calctime = real(now - t)/real(clock_rate)

        if (calctime < 60) then
            write(*,'(/a,f7.3,a)')'Time elapsed: ', calctime, ' s'
        elseif (calctime < 3600) then
            write(*,'(/a,i3,a,f7.3,a)')'Time elapsed: ', floor(calctime/60d0), ' min ', mod(calctime, 60d0), ' s'
        else
            write(*,'(/a,i3,a,i3,a,f7.3,a)')'Time elapsed: ', floor(calctime/3600d0), ' h ', floor(mod(calctime, 3600d0)/60d0), ' min ', mod(calctime, 60d0), ' s'
        endif
        
    end subroutine tock
    
end module
