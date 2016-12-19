!##############################################################################
! MODULE sorting
! 
! Module that contains different sorting routines.
!
! copyright: Fabian Kindermann
!            University of Wuerzburg
!            kindermann.fabian@uni-wuerzburg.de
!
!##############################################################################
module sorting


!##############################################################################
! Declaration of modules used
!##############################################################################

! for assertion of equality in dimensions
use toolbox

implicit none

save


!##############################################################################
! Interface declarations
!##############################################################################


!##############################################################################
! INTERFACE bubble_sort
! 
! Sorts one or two arrays by means of bubble sort. Bubble sort is pretty slow,
!     but nearly needs no additional space. If bubble sort is called with
!     two arrays, the first array is sorted and the second is changed in the
!     same way.
!##############################################################################
interface bubble_sort

    ! define methods used
    module procedure bubble_sort1_r, bubble_sort2_r, &
        bubble_sort1_i, bubble_sort2_i
        
end interface


!##############################################################################
! INTERFACE quick_sort
! 
! Sorts one or two arrays by means of quick sort. Quick sort is on average
!     one of the fastest sorting routines. However, on worst case it might
!     be as slow as bubble sort. In opposite to bubble sort, quick sort is
!     recursive
!##############################################################################
interface quick_sort

    ! define methods used
    module procedure quick_sort1_r, quick_sort2_r, &
        quick_sort1_i , quick_sort2_i
        
end interface




!##############################################################################
! Subroutines and functions                                               
!##############################################################################

contains


!##############################################################################
! SUBROUTINE bubble_sort1_r
! 
! Sorts an array by means of bubble_sort.
!##############################################################################
subroutine bubble_sort1_r(a)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    real*8, intent(inout) :: a(:)
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: j1, j2, n
    real*8 :: dummy
    
        
    !##### ROUTINE CODE #######################################################
	
	! get size of array
	n = size(a, 1)
	
	! sucessively sort the array
	do j1 = n, 2, -1
        do j2 = 1, j1-1
        
            ! if left entry greater than right entry -> change
            if(a(j2) > a(j2+1))then
                dummy = a(j2)
                a(j2) = a(j2+1)
                a(j2+1) = dummy
            endif
        enddo
    enddo
    
    ! check whether array is sorted
    call is_sorted_r(a, 'bubble_sort')

end subroutine bubble_sort1_r


!##############################################################################
! SUBROUTINE bubble_sort2_r
! 
! Sorts two arrays by means of bubble_sort.
!##############################################################################
subroutine bubble_sort2_r(a, b)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    real*8, intent(inout) :: a(:)
    
    ! additional array that will be changed in the same way
    real*8, intent(inout) :: b(:)
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: j1, j2, n
    real*8 :: dummy
    
        
    !##### ROUTINE CODE #######################################################
	
	! get size of array
	n = assert_eq(size(a, 1), size(b, 1), 'bubble_sort')
	
	! sucessively sort the array
	do j1 = n, 2, -1
        do j2 = 1, j1-1
        
            ! if left entry greater than right entry -> change
            if(a(j2) > a(j2+1))then
                dummy = a(j2)
                a(j2) = a(j2+1)
                a(j2+1) = dummy
                
                dummy = b(j2)
                b(j2) = b(j2+1)
                b(j2+1) = dummy
            endif
        enddo
    enddo
    
    ! check whether array is sorted
    call is_sorted_r(a, 'bubble_sort')

end subroutine bubble_sort2_r


!##############################################################################
! SUBROUTINE bubble_sort1_i
! 
! Sorts an array by means of bubble_sort.
!##############################################################################
subroutine bubble_sort1_i(a)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    integer, intent(inout) :: a(:)
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: j1, j2, n
    integer :: dummy
    
        
    !##### ROUTINE CODE #######################################################
	
	! get size of array
	n = size(a, 1)
	
	! sucessively sort the array
	do j1 = n, 2, -1
        do j2 = 1, j1-1
        
            ! if left entry greater than right entry -> change
            if(a(j2) > a(j2+1))then
                dummy = a(j2)
                a(j2) = a(j2+1)
                a(j2+1) = dummy             
            endif
        enddo
    enddo
    
    ! check whether array is sorted
    call is_sorted_i(a, 'bubble_sort')

end subroutine bubble_sort1_i


!##############################################################################
! SUBROUTINE bubble_sort2_i
! 
! Sorts two arrays by means of bubble_sort.
!##############################################################################
subroutine bubble_sort2_i(a, b)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    integer, intent(inout) :: a(:)
    
    ! additional array that will be changed in the same way
    integer, intent(inout) :: b(:)
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: j1, j2, n
    integer :: dummy
    
        
    !##### ROUTINE CODE #######################################################
	
	! get size of array
	n = assert_eq(size(a, 1), size(b, 1), 'bubble_sort')
	
	! sucessively sort the array
	do j1 = n, 2, -1
        do j2 = 1, j1-1
        
            ! if left entry greater than right entry -> change
            if(a(j2) > a(j2+1))then                
                dummy = a(j2)
                a(j2) = a(j2+1)
                a(j2+1) = dummy
                
                dummy = b(j2)
                b(j2) = b(j2+1)
                b(j2+1) = dummy
            endif
        enddo
    enddo

    ! check whether array is sorted
    call is_sorted_i(a, 'bubble_sort')

end subroutine bubble_sort2_i


!##############################################################################
! SUBROUTINE quick_sort1_r
! 
! Initializing subroutine for quicksort.
!##############################################################################
subroutine quick_sort1_r(a)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    real*8, intent(inout) :: a(:)    
    
        
    !##### ROUTINE CODE #######################################################
    
    ! call quick sort
    call quick_sort1_help_r(a, 1, size(a, 1))
    
    ! check whether array is sorted
    call is_sorted_r(a, 'quick_sort')
    
end subroutine quick_sort1_r


!##############################################################################
! SUBROUTINE quick_sort1_help_r
! 
! Sorts an array by means of quick_sort.
!##############################################################################
recursive subroutine quick_sort1_help_r(a, left, right)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    real*8, intent(inout) :: a(:)
    
    ! left element to start
    integer, intent(in) :: left
    
    ! right element to start
    integer, intent(in) :: right
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: l, r
    real*8 :: pivot, dummy
    
        
    !##### ROUTINE CODE #######################################################
	
	! check whether there is still something to sort
	if(left < right)then
	
	    ! get pivot element
	    pivot = a(right)
	    
	    ! set left and right counter
	    l = left
	    r = right - 1
	    
	    do
	        
	        ! increment left counter as long as element is lower
	        ! than pivot element
	        do while(a(l) <= pivot .and. l < right)
	            l = l + 1
	        enddo
	        
	        ! decrement right counter as long as element is greater
	        ! than pivot element
	        do while(a(r) >= pivot .and. r > left)
	            r = r - 1
	        enddo
	        
	        ! change l and r element if l < r, else there is nothing more
	        ! to bring into right order
	        if(l < r)then
	            dummy = a(l)
                a(l) = a(r)
                a(r) = dummy
	        endif
	        
	        if(l >= r)exit
	    
	    enddo
	    
	    ! put pivot element to the right position
	    dummy = a(l)
        a(l) = a(right)
        a(right) = dummy
	    
	    ! sort the left and right part of the array
	    call quick_sort1_help_r(a, left, l-1)
	    call quick_sort1_help_r(a, l+1, right)
	    
	endif

end subroutine quick_sort1_help_r



!##############################################################################
! SUBROUTINE quick_sort2_r
! 
! Initializing subroutine for quicksort.
!##############################################################################
subroutine quick_sort2_r(a, b)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    real*8, intent(inout) :: a(:)
    
    ! additional array that will be changed in the same way
    real*8, intent(inout) :: b(:)    
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: n
    
        
    !##### ROUTINE CODE #######################################################
    
    ! check for equal sizes
    n = assert_eq(size(a, 1), size(b, 1), 'quick_sort')
    
    ! call quick sort
    call quick_sort2_help_r(a, b, 1, size(a, 1))
    
    ! check whether array is sorted
    call is_sorted_r(a, 'quick_sort')
    
end subroutine quick_sort2_r


!##############################################################################
! SUBROUTINE quick_sort2_help_r
! 
! Sorts an array by means of quick_sort.
!##############################################################################
recursive subroutine quick_sort2_help_r(a, b, left, right)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    real*8, intent(inout) :: a(:)
    
    ! additional array that will be changed in the same way
    real*8, intent(inout) :: b(:)
    
    ! left element to start
    integer, intent(in) :: left
    
    ! right element to start
    integer, intent(in) :: right
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: l, r, n
    real*8 :: pivot, dummy
    
        
    !##### ROUTINE CODE #######################################################        
	
	! check whether there is still something to sort
	if(left < right)then
	
	    ! get pivot element
	    pivot = a(right)
	    
	    ! set left and right counter
	    l = left
	    r = right - 1
	    
	    do
	        
	        ! increment left counter as long as element is lower
	        ! than pivot element
	        do while(a(l) <= pivot .and. l < right)
	            l = l + 1
	        enddo
	        
	        ! decrement right counter as long as element is greater
	        ! than pivot element
	        do while(a(r) >= pivot .and. r > left)
	            r = r - 1
	        enddo
	        
	        ! change l and r element if l < r, else there is nothing more
	        ! to bring into right order
	        if(l < r)then
	            dummy = a(l)
                a(l) = a(r)
                a(r) = dummy
                
                dummy = b(l)
                b(l) = b(r)
                b(r) = dummy
	        endif
	        
	        if(l >= r)exit
	    
	    enddo
	    
	    ! put pivot element to the right position
	    dummy = a(l)
        a(l) = a(right)
        a(right) = dummy
	    
	    dummy = b(l)
        b(l) = b(right)
        b(right) = dummy
	    
	    ! sort the left and right part of the array
	    call quick_sort2_help_r(a, b, left, l-1)
	    call quick_sort2_help_r(a, b, l+1, right)
	    
	endif

end subroutine quick_sort2_help_r


!##############################################################################
! SUBROUTINE quick_sort1_i
! 
! Initializing subroutine for quicksort.
!##############################################################################
subroutine quick_sort1_i(a)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    integer, intent(inout) :: a(:)    
    
        
    !##### ROUTINE CODE #######################################################
    
    ! call quick sort
    call quick_sort1_help_i(a, 1, size(a, 1))
    
    ! check whether array is sorted
    call is_sorted_i(a, 'quick_sort')
    
end subroutine quick_sort1_i


!##############################################################################
! SUBROUTINE quick_sort1_help_i
! 
! Sorts an array by means of quick_sort.
!##############################################################################
recursive subroutine quick_sort1_help_i(a, left, right)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    integer, intent(inout) :: a(:)
    
    ! left element to start
    integer, intent(in) :: left
    
    ! right element to start
    integer, intent(in) :: right
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: l, r, pivot, dummy
    
        
    !##### ROUTINE CODE #######################################################
	
	! check whether there is still something to sort
	if(left < right)then
	
	    ! get pivot element
	    pivot = a(right)
	    
	    ! set left and right counter
	    l = left
	    r = right - 1
	    
	    do
	        
	        ! increment left counter as long as element is lower
	        ! than pivot element
	        do while(a(l) <= pivot .and. l < right)
	            l = l + 1
	        enddo
	        
	        ! decrement right counter as long as element is greater
	        ! than pivot element
	        do while(a(r) >= pivot .and. r > left)
	            r = r - 1
	        enddo
	        
	        ! change l and r element if l < r, else there is nothing more
	        ! to bring into right order
	        if(l < r)then
	            dummy = a(l)
                a(l) = a(r)
                a(r) = dummy
            endif
	        
	        if(l >= r)exit
	    
	    enddo
	    
	    ! put pivot element to the right position
	    dummy = a(l)
        a(l) = a(right)
        a(right) = dummy
	    
	    ! sort the left and right part of the array
	    call quick_sort1_help_i(a, left, l-1)
	    call quick_sort1_help_i(a, l+1, right)
	    
	endif

end subroutine quick_sort1_help_i



!##############################################################################
! SUBROUTINE quick_sort2_i
! 
! Initializing subroutine for quicksort.
!##############################################################################
subroutine quick_sort2_i(a, b)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    integer, intent(inout) :: a(:)
    
    ! additional array that will be changed in the same way
    integer, intent(inout) :: b(:)    
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: n
    
        
    !##### ROUTINE CODE #######################################################
    
    ! check for equal sizes
    n = assert_eq(size(a, 1), size(b, 1), 'quick_sort')
    
    ! call quick sort
    call quick_sort2_help_i(a, b, 1, size(a, 1))
    
    ! check whether array is sorted
    call is_sorted_i(a, 'quick_sort')
    
end subroutine quick_sort2_i


!##############################################################################
! SUBROUTINE quick_sort2_help_i
! 
! Sorts an array by means of quick_sort.
!##############################################################################
recursive subroutine quick_sort2_help_i(a, b, left, right)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! array that should be sorted
    integer, intent(inout) :: a(:)
    
    ! additional array that will be changed in the same way
    integer, intent(inout) :: b(:)
    
    ! left element to start
    integer, intent(in) :: left
    
    ! right element to start
    integer, intent(in) :: right
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: l, r, n, pivot, dummy
    
        
    !##### ROUTINE CODE #######################################################        
	
	! check whether there is still something to sort
	if(left < right)then
	
	    ! get pivot element
	    pivot = a(right)
	    
	    ! set left and right counter
	    l = left
	    r = right - 1
	    
	    do
	        
	        ! increment left counter as long as element is lower
	        ! than pivot element
	        do while(a(l) <= pivot .and. l < right)
	            l = l + 1
	        enddo
	        
	        ! decrement right counter as long as element is greater
	        ! than pivot element
	        do while(a(r) >= pivot .and. r > left)
	            r = r - 1
	        enddo
	        
	        ! change l and r element if l < r, else there is nothing more
	        ! to bring into right order
	        if(l < r)then
	            dummy = a(l)
                a(l) = a(r)
                a(r) = dummy
                
	            dummy = b(l)
                b(l) = b(r)
                b(r) = dummy
	        endif
	        
	        if(l >= r)exit
	    
	    enddo
	    
	    ! put pivot element to the right position
	    dummy = a(l)
        a(l) = a(right)
        a(right) = dummy
        
	    dummy = b(l)
        b(l) = b(right)
        b(right) = dummy
	    
	    ! sort the left and right part of the array
	    call quick_sort2_help_i(a, b, left, l-1)
	    call quick_sort2_help_i(a, b, l+1, right)
	    
	endif

end subroutine quick_sort2_help_i



!##############################################################################
! SUBROUTINE change_r
! 
! Changes two elements of an integer array.
!##############################################################################
subroutine change_r(a, e1, e2)


    !##### INPUT/OUTPUT VARIABLES #############################################

    ! respective array
    real*8, intent(inout) :: a(:)
    
    ! first element to change
    integer, intent(in) :: e1
    
    ! second element to change
    integer, intent(in) :: e2
    
    
    !##### OTHER VARIABLES ####################################################
    
    real*8 :: dummy
    
        
    !##### ROUTINE CODE #######################################################
    
    dummy   = a(e1)
    a(e1)   = a(e2)
    a(e2)   = dummy

end subroutine change_r


!##############################################################################
! SUBROUTINE change_i
! 
! Changes two elements of an integer array.
!##############################################################################
subroutine change_i(a, e1, e2)


    !##### INPUT/OUTPUT VARIABLES #############################################

    ! respective array
    integer, intent(inout) :: a(:)
    
    ! first element to change
    integer, intent(in) :: e1
    
    ! second element to change
    integer, intent(in) :: e2
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: dummy
    
        
    !##### ROUTINE CODE #######################################################
    
    dummy   = a(e1)
    a(e1)   = a(e2)
    a(e2)   = dummy

end subroutine change_i


!##############################################################################
! SUBROUTINE is_sorted_r
! 
! Checks whether an array is sorted.
!##############################################################################
subroutine is_sorted_r(a, routine)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! respective array
    real*8, intent(in) :: a(:)
    
    ! the routine name
    character(len = *) :: routine
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: j
    
        
    !##### ROUTINE CODE #######################################################
    
    do j = 1, size(a, 1)-1
    
        ! check whether a(j) <= a(j+1)
        if(a(j) > a(j+1))then
            call error(routine, 'sorting was not succesful')
        endif
    
    enddo

end subroutine is_sorted_r


!##############################################################################
! SUBROUTINE is_sorted_i
! 
! Checks whether an array is sorted.
!##############################################################################
subroutine is_sorted_i(a, routine)

    implicit none

    !##### INPUT/OUTPUT VARIABLES #############################################

    ! respective array
    integer, intent(in) :: a(:)
    
    ! the routine name
    character(len = *) :: routine
    
    
    !##### OTHER VARIABLES ####################################################
    
    integer :: j
    
        
    !##### ROUTINE CODE #######################################################
    
    do j = 1, size(a, 1)-1
    
        ! check whether a(j) <= a(j+1)
        if(a(j) > a(j+1))then
            call error(routine, 'sorting was not succesful')
        endif
    
    enddo

end subroutine is_sorted_i

end module sorting
