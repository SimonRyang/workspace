!##############################################################################
!##############################################################################
! MODULE toolbox
!
! Contains all the toolbox subroutines etc.
!
! copyright: Hans Fehr and Fabian Kindermann
!            University of Wuerzburg
!            contact@ce-fortran.com
!##############################################################################
!##############################################################################
module toolbox


!##############################################################################
!##############################################################################
! Declaration of Variables
!##############################################################################
!##############################################################################

implicit none

! starting time for cpu timer
real*8, private :: starttime_cpu

! coefficients for normal distribution
real*8, private :: tbox_ncoeffs(503)

! coefficients for inverse normal distribution
real*8, private :: tbox_icoeffs(503)

! derivative of normal cdf at 5
real*8, private :: tbox_dernormal

! derivative of inverse normal cdf at 0.999999
real*8, private :: tbox_dernormalInv

! is normal distribution loaded
logical, private :: tbox_normal_loaded = .false.

! should the random tbox_seed be set
logical, private :: tbox_seed = .true.

! Level of tolerance for all routines
real*8,  private  :: tbox_gftol = 1e-10

! Maximum number of iterations
integer, private  :: tbox_itermax_min = 1000

! Maximum number of iterations for brent_pow
integer, parameter, private  :: tbox_tbox_itermax_pow_b = 1000

! Level of tolerance for all routines
real*8,  private  :: tbox_gftol_root = 1e-8

! Maximum number of iterations for broydn
integer, private  :: itermax_root = 200

! control variables for gnuplot
logical, private :: gnu_addtoplot = .false.
logical, private :: gnu_dolegend = .false.
logical, private :: gnu_histogram = .false.
integer, private :: gnu_nmax, gnu_mmax

! plot data for gnuplot
real*8, allocatable, private :: gnu_x(:, :), gnu_y(:, :)
real*8, allocatable, private :: gnu_x_temp(:, :), gnu_y_temp(:, :)

! style data for gnuplot
character(LEN = 2000), dimension(1000) :: gnu_definitions




!##############################################################################
!##############################################################################
! Interface declarations
!##############################################################################
!##############################################################################


!##############################################################################
! INTERFACE assert_eq
!
! Interface for equality assertions by assert_eqx functions.
!##############################################################################
interface assert_eq

    module procedure assert_eq2, assert_eq3, assert_eq4, assert_eq5, &
        assert_eqn

end interface


!##############################################################################
! INTERFACE mat_mult
!
! Multiplies up to 6 matrices.
!##############################################################################
interface mat_mult

    module procedure mat_mult2, mat_mult3, mat_mult4, mat_mult5, mat_mult6

end interface


!##############################################################################
! INTERFACE normal_discrete
!
! Discretizes normal distribution.
!##############################################################################
interface normal_discrete

    module procedure normal_discrete_1, normal_discrete_2

end interface


!##############################################################################
! INTERFACE log_normal_discrete
!
! Discretizes log-normal distribution.
!##############################################################################
interface log_normal_discrete

    module procedure log_normal_discrete_1, log_normal_discrete_2

end interface


!##############################################################################
! INTERFACE simulate_uniform
!
! Simulates uniformly distributed random variables.
!##############################################################################
interface simulate_uniform

    module procedure simulate_uniform_1, simulate_uniform_n

end interface


!##############################################################################
! INTERFACE simulate_normal
!
! Simulates normallly distributed random variables.
!##############################################################################
interface simulate_normal

    module procedure simulate_normal_1, simulate_normal_n

end interface


!##############################################################################
! INTERFACE simulate_log-normal
!
! Simulates log-normallly distributed random variables.
!##############################################################################
interface simulate_log_normal

    module procedure simulate_log_normal_1, simulate_log_normal_n

end interface


!##############################################################################
! INTERFACE simulate_Gamma
!
! Simulates Gamma distributed random variables.
!##############################################################################
interface simulate_Gamma

    module procedure simulate_Gamma_1, simulate_Gamma_n

end interface


!##############################################################################
! INTERFACE simulate_beta
!
! Simulates beta distributed random variables.
!##############################################################################
interface simulate_beta

    module procedure simulate_beta_1, simulate_beta_n

end interface


!##############################################################################
! INTERFACE poly_interpol
!
! Interface for onedimensional polynomial interpolation.
!##############################################################################
interface poly_interpol

    module procedure poly_interpol_1, poly_interpol_m

end interface


!##############################################################################
! INTERFACE fminsearch
!
! Finds minimum of a function.
!##############################################################################
interface fminsearch

    module procedure brent, powell

end interface


!##############################################################################
! INTERFACE fzero
!
! Finds root of a function.
!##############################################################################
interface fzero

    module procedure newton_interpol, broydn

end interface


!##############################################################################
! INTERFACE grid_Inv_Equi
!
! Interface for inverting gridpoints.
!##############################################################################
interface grid_Inv_Equi

    module procedure grid_Inv_Equi_1, grid_Inv_Equi_m

end interface


!##############################################################################
! INTERFACE grid_Inv_Cheb
!
! Interface for inverting gridpoints.
!##############################################################################
interface grid_Inv_Cheb

    module procedure grid_Inv_Cheb_1, grid_Inv_Cheb_m

end interface


!##############################################################################
! INTERFACE grid_Inv_Grow
!
! Interface for inverting gridpoints.
!##############################################################################
interface grid_Inv_Grow

    module procedure grid_Inv_Grow_1, grid_Inv_Grow_m

end interface


!##############################################################################
! INTERFACE linint_Equi
!
! Interface for inverting gridpoints.
!##############################################################################
interface linint_Equi

    module procedure linint_Equi_1, linint_Equi_m

end interface


!##############################################################################
! INTERFACE linint_Cheb
!
! Interface for inverting gridpoints.
!##############################################################################
interface linint_Cheb

    module procedure linint_Cheb_1, linint_Cheb_m

end interface


!##############################################################################
! INTERFACE linint_Grow
!
! Interface for inverting gridpoints.
!##############################################################################
interface linint_Grow

    module procedure linint_Grow_1, linint_Grow_m

end interface


!##############################################################################
! INTERFACE spline_interp
!
! Interface for one- or multidimensional spline interpolation.
!##############################################################################
interface spline_interp

    module procedure spline_interp1, spline_interp2, spline_interp3, &
        spline_interp4, spline_interp5, spline_interp6, spline_interp7

end interface


!##############################################################################
! INTERFACE spline_eval
!
! Interface for evaluation of one- or multidimensional spline.
!##############################################################################
interface spline_eval

    module procedure &
        spline1, spline1_grid, &
        spline2, spline2_grid, &
        spline3, spline3_grid, &
        spline4, spline4_grid, &
        spline5, spline5_grid, &
        spline6, spline6_grid, &
        spline7, spline7_grid

end interface


!##############################################################################
! INTERFACE spline
!
! Interface for complete interpolation and evaluation of one- or
!     multidimensional spline.
!##############################################################################
interface spline

    module procedure &
        spline1_complete, spline1_complete_m , &
        spline2_complete, spline2_complete_m, &
        spline3_complete, spline3_complete_m , &
        spline4_complete, spline4_complete_m, &
        spline5_complete, spline5_complete_m, &
        spline6_complete, spline6_complete_m, &
        spline7_complete, spline7_complete_m

end interface


contains















!##############################################################################
!##############################################################################
! MODULE errwarn
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE error
    !
    ! Throws error message and stops program.
    !##############################################################################
    subroutine error(routine, message)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! routine in which error occured
        character(len=*), intent(in) :: routine

        ! error message
        character(len=*), intent(in) :: message


        !##### ROUTINE CODE #######################################################

        ! write error message
        write(*,'(/a,a,a,a/)')'ERROR ',routine,': ',message

        ! stop program
        stop

    end subroutine error


    !##############################################################################
    ! SUBROUTINE warning
    !
    ! Throws warning message
    !##############################################################################
    subroutine warning(routine, message)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! routine in which warning occured
        character(len=*), intent(in) :: routine

        ! warning message
        character(len=*), intent(in) :: message


        !##### ROUTINE CODE #######################################################

        ! write warning message
        write(*,'(/a,a,a,a/)')'WARNING ',routine,': ',message

    end subroutine warning















!##############################################################################
!##############################################################################
! MODULE assertions
!##############################################################################
!##############################################################################


    !##############################################################################
    ! FUNCTION assert_eq2
    !
    ! Checks equality for two integers.
    !##############################################################################
    function assert_eq2(n1, n2, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer, intent(in) :: n1, n2

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer :: assert_eq2


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (n1 == n2)then
            assert_eq2 = n1

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq2')
        end if

    end function assert_eq2


    !##############################################################################
    ! FUNCTION assert_eq3
    !
    ! Checks equality for three integers.
    !##############################################################################
    function assert_eq3(n1, n2, n3, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer, intent(in) :: n1, n2, n3

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer :: assert_eq3


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3)then
            assert_eq3 = n1

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq3')
        end if

    end function assert_eq3


    !##############################################################################
    ! FUNCTION assert_eq4
    !
    ! Checks equality for four integers.
    !##############################################################################
    function assert_eq4(n1, n2, n3, n4, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer, intent(in) :: n1, n2, n3, n4

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer :: assert_eq4


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4)then
            assert_eq4 = n1

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq4')
        end if

    end function assert_eq4


    !##############################################################################
    ! FUNCTION assert_eq5
    !
    ! Checks equality for five integers.
    !##############################################################################
    function assert_eq5(n1, n2, n3, n4, n5, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer, intent(in) :: n1, n2, n3, n4, n5

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer :: assert_eq5


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4 .and. n4 == n5)then
            assert_eq5 = n1

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eq5')
        end if

    end function assert_eq5


    !##############################################################################
    ! FUNCTION assert_eqn
    !
    ! Checks equality for n integers.
    !##############################################################################
    function assert_eqn(nn, string)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! integers to compare for equality
        integer, intent(in) :: nn(:)

        ! routine from which error should be thrown
        character(len=*), intent(in) :: string

        ! return value
        integer :: assert_eqn


        !##### ROUTINE CODE #######################################################

        ! if equality, set return value to n1
        if (all(nn(2:) == nn(1)))then
            assert_eqn = nn(1)

        ! else throw error message
        else
            call error(string, 'an assertion failed in assert_eqn')
        end if

    end function assert_eqn















!##############################################################################
!##############################################################################
! MODULE clock
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE tic
    !
    ! Starts cpu timer.
    !##############################################################################
    subroutine tic()


        !##### ROUTINE CODE #######################################################

        ! get cpu time
        call cpu_time(starttime_cpu)

    end subroutine tic


    !##############################################################################
    ! SUBROUTINE toc
    !
    ! Stops cpu timer.
    !##############################################################################
    subroutine toc(file)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! optional file identifier
        integer, intent(in), optional :: file


        !##### OTHER VARIABLES ####################################################

        real*8 :: time
        integer :: outfile
        real*8 :: times(4)


        !##### ROUTINE CODE #######################################################

        ! get output file identifier
        if(present(file))then
            outfile = file
        else
            outfile = 0
        endif

        ! get cpu time
        call cpu_time(time)

        ! calculate time difference
        time = time - starttime_cpu

        ! get number of days
        times(1) = floor(time/(24d0*60d0*60d0))
        time = time - times(1)*24d0*60d0*60d0

        ! get number of hours
        times(2) = floor(time/(60d0*60d0))
        time = time - times(2)*60d0*60d0

        ! get number of minutes
        times(3) = floor(time/60d0)
        time = time - times(3)*60d0

        ! get number of seconds
        times(4) = time

        call outTime(times, outfile)

    end subroutine toc


    !##############################################################################
    ! SUBROUTINE outTime
    !
    ! Writes time to file.
    !##############################################################################
    subroutine outTime(times, file)

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! time as integer array
        real*8, intent(in) :: times(4)

        ! the output file identifier
        integer, intent(in) :: file

        !##### OTHER VARIABLES ####################################################

        character(len=200) :: output1, output2


        !##### ROUTINE CODE #######################################################

        ! set up output
        write(output1, '(a)')'Time elapsed: '
        output2 = output1

        ! write time values
        if(times(1) > 0d0)then
            write(output1, '(a,1x,i3,a)') trim(output2), int(times(1)), ' d  '
            output2 = output1
        endif
        if(times(2) > 0d0)then
            write(output1, '(a,1x,i3,a)') trim(output2), int(times(2)), ' h  '
            output2 = output1
        endif
        if(times(3) > 0d0) then
            write(output1, '(a,1x,i3,a)')trim(output2), int(times(3)), ' min  '
            output2 = output1
        endif

        write(output1, '(a,1x,f7.3,a)')trim(output2), times(4), ' s  '

        if(file > 0) then
            write(file, '(/a/)')trim(output1)
        else
            write(*, '(/a/)')trim(output1)
        endif

    end subroutine outTime














!##############################################################################
!##############################################################################
! MODULE clock
!##############################################################################
!##############################################################################


    !##############################################################################
    ! FUNCTION grid_Cons_Equi
    !
    ! Constructs a whole equidistant grid on [left,right].
    !##############################################################################
    function grid_Cons_Equi(left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value at the grid point x \in [0,n]
        real*8 :: grid_Cons_Equi(0:n)


        !##### OTHER VARIABLES ####################################################

        real*8 :: h
        integer :: j


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Cons_Equi', &
            'left interval point greater than right point')

        ! calculate distance between grid points
        h = (right-left)/n

        ! calculate grid value
        grid_Cons_Equi = h*(/(dble(j), j=0,n)/)+left

    end function grid_Cons_Equi


    !##############################################################################
    ! FUNCTION grid_Cons_Cheb
    !
    ! Constructs a whole grid on [left,right] using chebychev nodes.
    !##############################################################################
    function grid_Cons_Cheb(left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value at the grid point x \in [0,n]
        real*8 :: grid_Cons_Cheb(0:n)


        !##### OTHER VARIABLES ####################################################

        integer :: j
        real*8, parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Cons_Cheb', &
            'left interval point greater than right point')

        ! calculate grid value
        grid_Cons_Cheb = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-(/(dble(j), j=0,n)/)+0.5d0)/dble(n+1)*pi)

    end function grid_Cons_Cheb


    !##############################################################################
    ! FUNCTION grid_Cons_Grow
    !
    ! Constructs a growing grid on [left, right].
    !##############################################################################
    function grid_Cons_Grow(left, right, growth, n)

        implicit none

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! growth rate
        real*8, intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value at the grid point x \in [0,n]
        real*8 :: grid_Cons_Grow(0:n)


        !##### OTHER VARIABLES ####################################################

        real*8 :: h
        integer :: j


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Cons_Grow', &
            'left interval point greater than right point')

        ! check for growth
        if(growth <= 0d0)call error('grid_Cons_Grow', &
            'growth rate must be greater than zero')

        ! calculate factor
        h = (right-left)/((1+growth)**n-1)

        ! calculate grid value
        grid_Cons_Grow = h*((1+growth)**(/(dble(j), j=0,n)/)-1)+left


    end function grid_Cons_Grow


    !##############################################################################
    ! FUNCTION grid_Val_Equi
    !
    ! Calculates single gridpoint of an equidistant grid.
    !##############################################################################
    function grid_Val_Equi(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value at the grid point x \in [0,n]
        real*8 :: grid_Val_Equi


        !##### OTHER VARIABLES ####################################################

        real*8 :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Val_Equi', &
            'left interval point greater than right point')

        ! calculate distance between grid points
        h = (right-left)/n

        ! calculate grid value
        grid_Val_Equi = h*x+left

    end function grid_Val_Equi


    !##############################################################################
    ! FUNCTION grid_Val_Cheb
    !
    ! Calculates single gridpoint of a Chebychev grid.
    !##############################################################################
    function grid_Val_Cheb(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value at the grid point x \in [0,n]
        real*8 :: grid_Val_Cheb


        !##### OTHER VARIABLES ####################################################

        real*8, parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Val_Cheb', &
            'left interval point greater than right point')

        ! calculate grid value
        grid_Val_Cheb = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-x+0.5d0)/dble(n+1)*pi)

    end function grid_Val_Cheb


    !##############################################################################
    ! FUNCTION grid_Val_Grow
    !
    ! Calculates single gridpoint of a growing grid.
    !##############################################################################
    function grid_Val_Grow(x, left, right, growth, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! growth rate
        real*8, intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value at the grid point x \in [0,n]
        real*8 :: grid_Val_Grow


        !##### OTHER VARIABLES ####################################################

        real*8 :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Val', &
            'left interval point greater than right point')

        ! check for growth
        if(growth <= 0d0)call error('grid_Val_Grow', &
            'growth rate must be greater than zero')

        ! calculate factor
        h = (right-left)/((1+growth)**n-1)

        ! calculate grid value
        grid_Val_Grow = h*((1+growth)**x-1)+left

    end function grid_Val_Grow


    !##############################################################################
    ! FUNCTION grid_Inv_Equi_1
    !
    ! Calculates inverse of gridpoints of an equidistant grid.
    !##############################################################################
    function grid_Inv_Equi_1(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real*8 :: grid_Inv_Equi_1


        !##### OTHER VARIABLES ####################################################

        real*8 :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Equi', &
            'left interval point greater than right point')

        ! calculate distance between grid points
        h = (right-left)/n

        ! calculate grid value
        grid_Inv_Equi_1 = (x-left)/h

    end function grid_Inv_Equi_1


    !##############################################################################
    ! FUNCTION grid_Inv_Cheb_1
    !
    ! Calculates inverse of gridpoints of a Chebychev grid.
    !##############################################################################
    function grid_Inv_Cheb_1(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real*8 :: grid_Inv_Cheb_1


        !##### OTHER VARIABLES ####################################################

        real*8, parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Cheb', &
            'left interval point greater than right point')

        ! calculate grid value
        grid_Inv_Cheb_1 = dble(n) + 0.5d0 - acos((2d0*x-(left+right)) &
            /(right-left))*dble(n+1)/pi

    end function grid_Inv_Cheb_1


    !##############################################################################
    ! FUNCTION grid_Inv_Grow_1
    !
    ! Calculates inverse of gridpoints of a growing grid.
    !##############################################################################
    function grid_Inv_Grow_1(x, left, right, growth, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! growth rate
        real*8, intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real*8 :: grid_Inv_Grow_1


        !##### OTHER VARIABLES ####################################################

        real*8 :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Grow', &
            'left interval point greater than right point')

        ! check for growth
        if(growth <= 0d0)call error('grid_Inv_Grow', &
            'growth rate must be greater than zero')

        ! calculate factor
        h = (right-left)/((1+growth)**n-1d0)

        ! calculate grid value
        grid_Inv_Grow_1 = log((x-left)/h+1d0)/log(1d0+growth)

    end function grid_Inv_Grow_1


    !##############################################################################
    ! FUNCTION grid_Inv_Equi_m
    !
    ! Calculates inverse of gridpoints of an equidistant grid.
    !##############################################################################
    function grid_Inv_Equi_m(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x(:)

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real*8 :: grid_Inv_Equi_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Equi', &
            'left interval point greater than right point')

        ! calculate distance between grid points
        h = (right-left)/n

        ! calculate grid value
        grid_Inv_Equi_m = (x-left)/h

    end function grid_Inv_Equi_m


    !##############################################################################
    ! FUNCTION grid_Inv_Cheb_m
    !
    ! Calculates inverse of gridpoints of a Chebychev grid.
    !##############################################################################
    function grid_Inv_Cheb_m(x, left, right, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x(:)

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real*8 :: grid_Inv_Cheb_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8, parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Cheb', &
            'left interval point greater than right point')

        ! calculate grid value
        grid_Inv_Cheb_m = dble(n) + 0.5d0 - acos((2d0*x-(left+right)) &
            /(right-left))*dble(n+1)/pi

    end function grid_Inv_Cheb_m


    !##############################################################################
    ! FUNCTION grid_Inv_Grow_m
    !
    ! Calculates inverse of gridpoints of a growing grid.
    !##############################################################################
    function grid_Inv_Grow_m(x, left, right, growth, n)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x(:)

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! growth rate
        real*8, intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! value of the inverse of the gridpoint x \in [left, right]
        real*8 :: grid_Inv_Grow_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: h


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('grid_Inv_Grow', &
            'left interval point greater than right point')

        ! check for growth
        if(growth <= 0d0)call error('grid_Inv_Grow', &
            'growth rate must be greater than zero')

        ! calculate factor
        h = (right-left)/((1d0+growth)**n-1d0)

        ! calculate grid value
        grid_Inv_Grow_m = log((x-left)/h+1d0)/log(1d0+growth)

    end function grid_Inv_Grow_m


    !##############################################################################
    ! subroutine linint_Equi_1
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
    subroutine linint_Equi_1(x, left, right, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! left interpolation point
        integer, intent(out) :: il

        ! right interpolation point
        integer, intent(out) :: ir

        ! interpolation fraction
        real*8, intent(out) :: phi

        !##### OTHER VARIABLES ####################################################

        real*8 :: h, xinv, xl, xr


        !##### ROUTINE CODE #######################################################

        ! invert the grid to get point
        xinv = grid_Inv_Equi_1(min(max(x, left), right), left, right, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        h  = (right-left)/n
        xl = h*dble(il)+left
        xr = h*dble(ir)+left

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Equi_1


    !##############################################################################
    ! subroutine linint_Equi_m
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
    subroutine linint_Equi_m(x, left, right, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x(:)

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! left interpolation point
        integer, intent(out) :: il(:)

        ! right interpolation point
        integer, intent(out) :: ir(:)

        ! interpolation fraction
        real*8, intent(out) :: phi(:)

        !##### OTHER VARIABLES ####################################################

        integer :: m
        real*8 :: h, xinv(size(x, 1)), xl(size(x, 1)), xr(size(x, 1))


        !##### ROUTINE CODE #######################################################

        ! check for sizes
        m = assert_eq(size(x, 1), size(il, 1), size(ir, 1), size(phi, 1), 'linint_Equi')
        m = m

        ! invert the grid to get point
        xinv = grid_Inv_Equi_m(min(max(x, left), right), left, right, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        h  = (right-left)/n
        xl = h*dble(il)+left
        xr = h*dble(ir)+left

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Equi_m


    !##############################################################################
    ! subroutine linint_Cheb_1
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
    subroutine linint_Cheb_1(x, left, right, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! left interpolation point
        integer, intent(out) :: il

        ! right interpolation point
        integer, intent(out) :: ir

        ! interpolation fraction
        real*8, intent(out) :: phi

        !##### OTHER VARIABLES ####################################################

        real*8 :: xinv, xl, xr
        real*8, parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! invert the grid to get point
        xinv = grid_Inv_Cheb_1(min(max(x, left), right), left, right, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        xl = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-dble(il)+0.5d0)/dble(n+1)*pi)
        xr = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-dble(ir)+0.5d0)/dble(n+1)*pi)

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Cheb_1


    !##############################################################################
    ! subroutine linint_Cheb_m
    !
    ! Calculates linear interpolant on a Chebychev grid.
    !##############################################################################
    subroutine linint_Cheb_m(x, left, right, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x(:)

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! left interpolation point
        integer, intent(out) :: il(:)

        ! right interpolation point
        integer, intent(out) :: ir(:)

        ! interpolation fraction
        real*8, intent(out) :: phi(:)

        !##### OTHER VARIABLES ####################################################

        integer :: m
        real*8 :: xinv(size(x, 1)), xl(size(x, 1)), xr(size(x, 1))
        real*8, parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! check for sizes
        m = assert_eq(size(x, 1), size(il, 1), size(ir, 1), size(phi, 1), 'linint_Equi')
        m = m

        ! invert the grid to get point
        xinv = grid_Inv_Cheb_m(min(max(x, left), right), -1d0, 1d0, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        xl = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-dble(il)+0.5d0)/dble(n+1)*pi)
        xr = (left+right)/2d0 + (right-left)/2d0* &
            cos((dble(n)-dble(ir)+0.5d0)/dble(n+1)*pi)

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Cheb_m


    !##############################################################################
    ! subroutine linint_Grow_1
    !
    ! Calculates linear interpolant on a growing grid.
    !##############################################################################
    subroutine linint_Grow_1(x, left, right, growth, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! growth rate
        real*8, intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! left interpolation point
        integer, intent(out) :: il

        ! right interpolation point
        integer, intent(out) :: ir

        ! interpolation fraction
        real*8, intent(out) :: phi

        !##### OTHER VARIABLES ####################################################

        real*8 :: h, xinv, xl, xr


        !##### ROUTINE CODE #######################################################

        ! check for left <= right
        if(left >= right)call error('linint_Grow', &
            'left interval point greater than right point')

        ! invert the grid to get point
        xinv = grid_Inv_Grow_1(min(max(x, left), right), left, right, growth, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        h = (right-left)/((1+growth)**n-1)
        xl = h*((1+growth)**dble(il)-1d0)+left
        xr = h*((1+growth)**dble(ir)-1d0)+left

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Grow_1


    !##############################################################################
    ! subroutine linint_Grow_m
    !
    ! Calculates linear interpolant on an equidistant grid.
    !##############################################################################
    subroutine linint_Grow_m(x, left, right, growth, n, il, ir, phi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point that shall be calculated
        real*8, intent(in) :: x(:)

        ! left and right interval point
        real*8, intent(in) :: left, right

        ! growth rate
        real*8, intent(in) :: growth

        ! last grid point: 0,1,...,n
        integer, intent(in) :: n

        ! left interpolation point
        integer, intent(out) :: il(:)

        ! right interpolation point
        integer, intent(out) :: ir(:)

        ! interpolation fraction
        real*8, intent(out) :: phi(:)

        !##### OTHER VARIABLES ####################################################

        integer :: m
        real*8 :: h, xinv(size(x, 1)), xl(size(x, 1)), xr(size(x, 1))


        !##### ROUTINE CODE #######################################################

        ! check for sizes
        m = assert_eq(size(x, 1), size(il, 1), size(ir, 1), size(phi, 1), 'linint_Equi')
        m = m

        ! check for left <= right
        if(left >= right)call error('linint_Grow', &
            'left interval point greater than right point')

        ! invert the grid to get point
        xinv = grid_Inv_Grow_m(min(max(x, left), right), left, right, growth, n)

        ! get left and right gridpoint
        il = min(max(floor(xinv), 0), n-1)
        ir = il+1

        ! determine left and right gridpoint
        h = (right-left)/((1+growth)**n-1)
        xl = h*((1+growth)**dble(il)-1d0)+left
        xr = h*((1+growth)**dble(ir)-1d0)+left

        ! get share on the left point
        phi = (xr-x)/(xr-xl)

    end subroutine linint_Grow_m





!##############################################################################
!##############################################################################
! MODULE matrixtools
!
! Large parts of the procedures were taken from:
!     Press, Teukolsky, Vetterling and Flannery (1992): "Numerical Recipes in
!     FORTRAN: The Art of Scientific Computing", 2nd edition, Cambridge
!     Univeristy Press, Cambridge.
!##############################################################################
!##############################################################################


    !##############################################################################
    ! FUNCTION mat_trans
    !
    ! Transposes a matrix.
    !##############################################################################
    function mat_trans(a)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix of the system
        real*8, intent(in) :: a(:, :)

        ! the inverse of the matrix
        real*8 :: mat_trans(size(a, 2), size(a, 1))


        !##### OTHER VARIABLES ####################################################

        integer :: j


        !##### ROUTINE CODE #######################################################

        ! sucessively transpose
        do j = 1, size(a,1)
            mat_trans(:, j) = a(j, :)
        enddo

    end function mat_trans


    !##############################################################################
    ! FUNCTION mat_mult2
    !
    ! Multiplies two matrices a*b.
    !##############################################################################
    function mat_mult2(a, b) result(mat)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix number one
        real*8, intent(in) :: a(:, :)

        ! matrix number two
        real*8, intent(in) :: b(:, :)

        ! the result
        real*8 :: mat(size(a, 1), size(b, 2))


        !##### OTHER VARIABLES ####################################################

        integer :: j, k, n


        !##### ROUTINE CODE #######################################################

        ! assert equality of dimensions
        n = assert_eq(size(a, 2), size(b, 1), 'mat_mult')
        n = n

        ! sucessively multiply rows and columns
        do j = 1, size(a,1)
            do k = 1, size(b, 2)
                mat(j, k) = sum(a(j, :)*b(:, k))
            enddo
        enddo

    end function mat_mult2


    !##############################################################################
    ! FUNCTION mat_mult3
    !
    ! Multiplies three matrices a*b*c.
    !##############################################################################
    function mat_mult3(a, b, c) result(mat)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix number one
        real*8, intent(in) :: a(:, :)

        ! matrix number two
        real*8, intent(in) :: b(:, :)

        ! matrix number three
        real*8, intent(in) :: c(:, :)

        ! the result
        real*8 :: mat(size(a, 1), size(c, 2))


        !##### OTHER VARIABLES ####################################################


        !##### ROUTINE CODE #######################################################

        ! sucessively multiply rows and columns
        mat = mat_mult2(a, mat_mult2(b, c))

    end function mat_mult3


    !##############################################################################
    ! FUNCTION mat_mult4
    !
    ! Multiplies four matrices a*b*c*d.
    !##############################################################################
    function mat_mult4(a, b, c, d) result(mat)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix number one
        real*8, intent(in) :: a(:, :)

        ! matrix number two
        real*8, intent(in) :: b(:, :)

        ! matrix number three
        real*8, intent(in) :: c(:, :)

        ! matrix number four
        real*8, intent(in) :: d(:, :)

        ! the result
        real*8 :: mat(size(a, 1), size(d, 2))


        !##### OTHER VARIABLES ####################################################


        !##### ROUTINE CODE #######################################################

        ! sucessively multiply rows and columns
        mat = mat_mult2(a, mat_mult3(b, c, d))

    end function mat_mult4


    !##############################################################################
    ! FUNCTION mat_mult5
    !
    ! Multiplies four matrices a*b*c*d*e.
    !##############################################################################
    function mat_mult5(a, b, c, d, e) result(mat)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix number one
        real*8, intent(in) :: a(:, :)

        ! matrix number two
        real*8, intent(in) :: b(:, :)

        ! matrix number three
        real*8, intent(in) :: c(:, :)

        ! matrix number four
        real*8, intent(in) :: d(:, :)

        ! matrix number five
        real*8, intent(in) :: e(:, :)

        ! the result
        real*8 :: mat(size(a, 1), size(e, 2))


        !##### OTHER VARIABLES ####################################################


        !##### ROUTINE CODE #######################################################

        ! sucessively multiply rows and columns
        mat = mat_mult2(a, mat_mult4(b, c, d, e))

    end function mat_mult5


    !##############################################################################
    ! FUNCTION mat_mult6
    !
    ! Multiplies four matrices a*b*c*d*e*f.
    !##############################################################################
    function mat_mult6(a, b, c, d, e, f) result(mat)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix number one
        real*8, intent(in) :: a(:, :)

        ! matrix number two
        real*8, intent(in) :: b(:, :)

        ! matrix number three
        real*8, intent(in) :: c(:, :)

        ! matrix number four
        real*8, intent(in) :: d(:, :)

        ! matrix number five
        real*8, intent(in) :: e(:, :)

        ! matrix number six
        real*8, intent(in) :: f(:, :)

        ! the result
        real*8 :: mat(size(a, 1), size(f, 2))


        !##### OTHER VARIABLES ####################################################


        !##### ROUTINE CODE #######################################################

        ! sucessively multiply rows and columns
        mat = mat_mult2(a, mat_mult5(b, c, d, e, f))

    end function mat_mult6



    !##############################################################################
    ! SUBROUTINE lu_solve
    !
    ! Solves a linear equation system by lu-decomposition.
    !##############################################################################
    subroutine lu_solve(a, b)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix of the system
        real*8, intent(in) :: a(:, :)

        ! right side of equation and solution of the system
        real*8, intent(inout) :: b(:)


        !##### OTHER VARIABLES ####################################################

        real*8 :: indx(size(b))
        real*8 :: worka(size(a, 1), size(a, 2))
        real*8 :: d
        integer :: n


        !##### ROUTINE CODE #######################################################

        ! assert size equality
        n = assert_eq(size(a,1), size(a,2), size(b), 'lu_solve')
        n = n

        ! copy matrix to working matrix
        worka = a

        ! decompose matrix
        call lu_decomp(worka, indx, d)

        ! solve system
        call lu_back(worka, indx, b)

    end subroutine lu_solve


    !##############################################################################
    ! FUNCTION lu_invert
    !
    ! Inverts a matrix by lu-decomposition.
    !##############################################################################
    function lu_invert(a)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix of the system
        real*8, intent(in) :: a(:, :)

        ! the inverse of the matrix
        real*8 :: lu_invert(size(a, 1), size(a, 2))


        !##### OTHER VARIABLES ####################################################

        integer :: j, n


        !##### ROUTINE CODE #######################################################

        ! assert size equality
        n = assert_eq(size(a,1), size(a,2), 'lu_invert')

        ! set up unity matrix
        lu_invert = 0d0
        do j = 1, n
            lu_invert(j, j) = 1d0
        enddo

        ! succesively solve the system with unity matrix
        do j = 1, n
            call lu_solve(a, lu_invert(:, j))
        enddo

    end function lu_invert


    !##############################################################################
    ! SUBROUTINE lu_decomp
    !
    ! Calculates lu-decomposition of matrices.
    !##############################################################################
    subroutine lu_decomp(a, indx, d)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix that shall be decomposed
        real*8, intent(inout) :: a(:, :)

        ! row permutation indicator due to pivoting
        real*8, intent(out) :: indx(:)

        ! indicates whether number of row permutations was even or odd
        real*8, intent(out) :: d


        !##### OTHER VARIABLES ####################################################

        real*8 :: vv(size(a,1))
        real*8, parameter :: tiny = 1.0e-20
        integer :: j, n, imax


        !##### ROUTINE CODE #######################################################

        ! check array sizes
        n = assert_eq(size(a,1), size(a,2), size(indx),'lu_decomp')

        ! initialize permutation indicator
        d = 1d0

        ! get maximum value in every row
        vv = maxval(abs(a), dim=2)

        ! if there is a zero row then matrix is singular
        if (any(abs(vv) <= 1d-100)) call error('lu_decomp', 'matrix is singular')

        ! invert v
        vv = 1d0/vv

        ! start lu-decomposition process
        do j = 1, n

            ! get index of pivot element
            imax = (j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))

            ! do pivoting if pivot element is not the first element
            if(j /= imax)then
                call swap(a(imax,:), a(j,:))
                d = -d
                vv(imax) = vv(j)
            endif

            ! indicate pivot element
            indx(j) = imax

            ! prevent division by 0
            if(abs(a(j,j))  <= 1d-100)call error('lu_decomp', 'matrix is singular')

            ! calculate new elements
            a(j+1:n, j) = a(j+1:n,j)/a(j,j)
            a(j+1:n, j+1:n) = a(j+1:n,j+1:n)-outerprod(a(j+1:n,j), a(j,j+1:n))
        enddo


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        function outerprod(a, b)

            real*8, intent(in) :: a(:), b(:)
            real*8 :: outerprod(size(a, 1),size(b, 1))

            outerprod = spread(a, dim=2, ncopies=size(b, 1)) * &
            spread(b, dim=1, ncopies=size(a, 1))

        end function outerprod


        subroutine swap(a, b)

            real*8, intent(inout) :: a(:), b(:)
            real*8 :: dum(size(a))

            dum = a
            a = b
            b = dum

        end subroutine swap


        function imaxloc(arr)

            real*8, intent(in) :: arr(:)
            integer :: imaxloc
            integer :: imax(1)

            imax = maxloc(arr(:))
            imaxloc = imax(1)

        end function imaxloc

    end subroutine lu_decomp



    !##############################################################################
    ! SUBROUTINE lu_dec
    !
    ! Calculates lu-decomposition of matrices and returns L and U matrix.
    !##############################################################################
    subroutine lu_dec(a, l, u)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! matrix that shall be decomposed
        real*8, intent(in) :: a(:, :)

        ! row permutation indicator due to pivoting
        real*8, intent(out) :: l(:, :)

        ! indicates whether number of row permutations was even of odd
        real*8, intent(out) :: u(:, :)


        !##### OTHER VARIABLES ####################################################

        real*8 :: Awork(size(a,1), size(a,2))
        real*8 :: indx(size(a,1))
        real*8 :: d
        integer :: n, j


        !##### ROUTINE CODE #######################################################

        ! check array sizes
        n = assert_eq((/size(a,1), size(a,2), size(l, 1), size(l, 2), &
        size(u, 1), size(u, 2)/), 'lu_dec')

        ! copy matrix
        Awork(:, :) = A(:, :)

        ! calculate decomposition
        call lu_decomp(Awork, indx, d)

        ! initialize matrices
        L(:, :) = 0d0
        U(:, :) = 0d0

        ! set up new matrices
        do j = 1, n

            ! diagonal element of L
            L(j, j) = 1d0

            ! other elements of L
            L(j, 1:j-1) = Awork(j, 1:j-1)

            ! elements of U
            U(j, j:n) = AWork(j, j:n)
        enddo

    end subroutine lu_dec


    !##############################################################################
    ! SUBROUTINE lu_back
    !
    ! Solves a lu decomposed matrix by backsubstitution.
    !##############################################################################
    subroutine lu_back(a, indx, b)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! lu-decomposed matrix that defines system
        real*8, intent(in) :: a(:, :)

        ! row permutation indicator due to pivoting
        real*8, intent(in) :: indx(:)

        ! right side of equation and solution of the system
        real*8, intent(inout) :: b(:)


        !##### OTHER VARIABLES ####################################################

        real*8 :: summ
        integer :: i, n, ii, ll


        !##### ROUTINE CODE #######################################################

        ! assert size equality
        n = assert_eq(size(a,1), size(a,2), size(indx), size(b), 'lu_back')

        ! start backward solving provess
        ii = 0
        do i = 1, n

            ll = indx(i)
            summ = b(ll)
            b(ll) = b(i)
            if(ii /= 0)then
                summ = summ-dot_product(a(i,ii:i-1), b(ii:i-1))
            elseif (abs(summ) >= 1d-100) then
                ii = i
            endif
            b(i)=summ
        enddo

        do i=n, 1, -1
            b(i) = (b(i)-dot_product(a(i,i+1:n), b(i+1:n)))/a(i,i)
        enddo

    end subroutine lu_back



    !##############################################################################
    ! SUBROUTINE cholesky
    !
    ! Calculates cholesky factorization of a symmetric matrix.
    !##############################################################################
    subroutine cholesky(a, l)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the matrix that should be decomposed
        real*8, intent(in) :: a(:, :)

        ! the cholesky factor
        real*8, intent(out) :: l(:, :)


        !##### OTHER VARIABLES ####################################################

        integer :: i, n
        real*8 :: summ, p(size(a,1))


        !##### ROUTINE CODE #######################################################

        ! assert equalities
        n = assert_eq(size(a,1), size(a,2), size(l, 1), size(l, 2), &
        size(p), 'normal_discrete')

        ! copy matrix
        l = a

        ! decompose matrix
        do i = 1, n
            summ = l(i,i)-dot_product(l(i,1:i-1), l(i,1:i-1))
            if(summ <= 0d0)call error('normal_discrete', &
                'Cholesky decomposition failed')
            p(i) = sqrt(summ)
            l(i+1:n,i) = (l(i,i+1:n)-matmul(l(i+1:n,1:i-1),l(i,1:i-1)))/p(i)
        enddo

        ! copy matrix
        do i = 1, n
            l(i, i) = p(i)
            l(i, i+1:n) = 0d0
        enddo

    end subroutine cholesky















!##############################################################################
!##############################################################################
! MODULE normalProb
! normal_discrete is taken from:
!     Miranda and Fackler (1992): "Applied Computational Economics and
!     Finance", MIT Press, 2002.
!##############################################################################
!##############################################################################


    !##############################################################################
    ! FUNCTION normalCDF
    !
    ! Calculates cumulated normal distribution at point p.
    !##############################################################################
    function normalCDF(p, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real*8, intent(in) :: p

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma

        ! value of the spline function at p
        real*8 :: normalCDF


        !##### OTHER VARIABLES ####################################################

        real*8 :: ptrans, ptemp, t, phi, mu_c, sigma_c
        integer :: l, r, k


        !##### ROUTINE CODE #######################################################

        ! test whether normal distribution was loaded
        if(.not. tbox_normal_loaded)call loadNormal()

        ! initialize expectation and variance
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sigma

        ! standardize evaluation point
        ptrans = (p - mu_c)/sqrt(sigma_c)

        ! map to right side if ptemp < 0
        if(ptrans < 0) ptrans = -ptrans

        ! store ptrans for other calculations
        ptemp = ptrans

        ! restrict evaluation point to interpolation interval
        ptemp = max(ptemp, 0d0)
        ptemp = min(ptemp, 5d0)

        ptemp = (ptemp - 0d0) / 0.01d0

        ! set up left and right calculation end point
        l = floor(ptemp) + 1
        r = min(l + 3, 503)

        normalCDF = 0d0

        do k = l, r

            ! set up point where to evaluate phi function
            t = abs(ptemp - k + 2)

            ! calculate spline function
            if(t < 1d0)then
                phi = 4d0 + t**2d0 * (3d0 * t - 6d0)
            elseif(t <= 2d0)then
                phi = (2d0 - t)**3d0
            else
                phi = 0d0
            endif

            ! calculate final spline value
            normalCDF = normalCDF + tbox_ncoeffs(k) * phi
        enddo

        ! extrapolate if ptrans > 5
        if(ptrans > 5d0) then
            normalCDF = normalCDF + tbox_dernormal * (ptrans - 5d0)
            normalCDF = min(normalCDF, 1d0)
        endif

        ! again invert if ptemp was < 0
        if(p - mu_c < 0d0) normalCDF = 1d0 - normalCDF

    end function normalCDF


    !##############################################################################
    ! FUNCTION normalPDF
    !
    ! Calculates normal density functions at point p.
    !##############################################################################
    function normalPDF(p, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real*8, intent(in) :: p

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma

        ! value of normal density at p
        real*8 :: normalPDF


        !##### OTHER VARIABLES ####################################################

        real*8 :: mu_c, sigma_c
        real*8, parameter :: pi = 3.1415926535897d0


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sqrt(sigma)

        normalPDF = 1/(sigma_c*sqrt(2*pi))*exp(-((p-mu_c)/sigma_c)**2/2)

    end function normalPDF


    !##############################################################################
    ! FUNCTION normalCDF_Inv
    !
    ! Calculates inverse cumulated normal distribution at point p.
    !##############################################################################
    function normalCDF_Inv(p, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point where to calculate function
        real*8 :: p

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma

        ! value of the spline function at p
        real*8 :: normalCDF_Inv


        !##### OTHER VARIABLES ####################################################

        real*8 :: ptemp, t, phi, mu_c, sigma_c
        integer :: l, r, k


        !##### ROUTINE CODE #######################################################

        ! test whether normal distribution was loaded
        if(.not. tbox_normal_loaded)call loadNormal()

        ! initialize expectation and variance
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sqrt(sigma)

        ! store for further calculations
        ptemp = p

        ! map to right side if ptemp < 0.5
        if(ptemp < 0.5d0) ptemp = 1d0 - ptemp

        ! restrict evaluation point to interpolation interval
        ptemp = max(ptemp, 0.5d0)
        ptemp = min(ptemp, 0.999999713348428d0)

        ptemp = (ptemp - 0.5d0) / (0.999999713348428d0 - 0.5d0) * 500d0

        ! set up left and right calculation end point
        l = floor(ptemp) + 1
        r = min(l + 3, 503)

        normalCDF_Inv = 0d0

        do k = l, r

            ! set up point where to evaluate phi function
            t = abs(ptemp - k + 2)

            ! calculate spline function
            if(t < 1d0)then
                phi = 4d0 + t**2d0 * (3d0 * t - 6d0)
            elseif(t <= 2d0)then
                phi = (2d0 - t)**3d0
            else
                phi = 0d0
            endif

            ! calculate final spline value
            normalCDF_Inv = normalCDF_Inv + tbox_icoeffs(k) * phi
        enddo

        ! extrapolate if ptrans > 4
        if(p > 0.999999713348428d0) then
            normalCDF_Inv = normalCDF_Inv + tbox_dernormalInv * (p - 0.999999713348428d0)
        endif

        ! again invert if ptemp was < 0
        if(p < 0.5d0) normalCDF_Inv = - normalCDF_Inv

        ! transfer to mu and sigma
        normalCDF_Inv = mu_c + sigma_c*normalCDF_Inv

    end function normalCDF_Inv


    !##############################################################################
    ! FUNCTION BinProb
    !
    ! Calculates probabilities of a normal distribution.
    !##############################################################################
    function BinProb(n, T, p)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! the number of positive draws
        integer, intent(in) :: n

        ! the total number of draws
        integer, intent(in) :: T

        ! the probability of a positive draw
        real*8, intent(in) :: p

        ! the probability of this to happen
        real*8 :: BinProb


        !##### OTHER VARIABLES ####################################################

        real*8 :: sk1, sk2, sk3
        integer :: i


        !##### ROUTINE CODE #######################################################


        if (n > T .or. p < 0d0) then
            Binprob = 0d0
            call warning('BinProb', 'Wrong input arguments!')
        else
            sk1 = 1d0
            sk2 = 1d0
            sk3 = 1d0
            do i = 2, n
                sk1 = sk1*i
            enddo
            do i = 2, T
              sk2 = sk2*i
            enddo
            do i = 2, T-n
              sk3 = sk3*i
            enddo
            Binprob = sk2/sk3/sk1*p**n*(1-p)**(T-n)
        endif

    end function


    !##############################################################################
    ! SUBROUTINE init_random_seed
    !
    ! To ensure that each ranrom draw is a new sequence (from Fortran doc)
    !##############################################################################
    subroutine init_random_seed()

        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, dt(8)
        integer, parameter :: int64 = selected_int_kind(16)
        integer(kind=int64) :: t

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(t)
        if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24_int64 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
        endif
        do i = 1, n
            seed(i) = lcg(t)
        enddo
        call random_seed(put=seed)

    contains

        ! This simple PRNG might not be good enough for real work, but is
        ! sufficient for seeding a better PRNG.
        function lcg(s)

            implicit none
            integer :: lcg
            integer(int64) :: s

            if (s == 0) then
                s = 104729
            else
                s = mod(s, 4294967296_int64)
            endif
            s = mod(s * 279470273_int64, 4294967291_int64)
            lcg = int(mod(s, int(huge(0), int64)), kind(0))

        end function lcg

    end subroutine init_random_seed


    !##############################################################################
    ! SUBROUTINE simulate_uniform_1
    !
    ! Simulates one draw from a uniform distribution.
    !##############################################################################
    subroutine simulate_uniform_1(x, a, b)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x

        ! left end of the distribution
        real*8, optional :: a

        ! right end of the distribution
        real*8, optional :: b


        !##### OTHER VARIABLES ####################################################

        real*8 :: a_c, b_c


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        a_c = 0d0
        if(present(a))a_c = a
        b_c = 1d0
        if(present(b))b_c = b

        ! initialize the random seed
        if(tbox_seed)then
            call init_random_seed()
            tbox_seed = .false.
        endif

        ! draw the random number
        call random_number(x)

        x = a_c + (b_c-a_c)*x

    end subroutine simulate_uniform_1


    !##############################################################################
    ! SUBROUTINE simulate_uniform_n
    !
    ! Simulates a series draw from a uniform distribution.
    !##############################################################################
    subroutine simulate_uniform_n(x, a, b)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x(:)

        ! left end of the distribution
        real*8, optional :: a

        ! right end of the distribution
        real*8, optional :: b


        !##### OTHER VARIABLES ####################################################

        real*8 :: a_c, b_c


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        a_c = 0d0
        if(present(a))a_c = a
        b_c = 1d0
        if(present(b))b_c = b

        ! initialize the random seed
        if(tbox_seed)then
            call init_random_seed()
            tbox_seed = .false.
        endif

        call random_number(x)

        x = a_c + (b_c-a_c)*x

    end subroutine simulate_uniform_n


    !##############################################################################
    ! Subroutine simulate_normal_1
    !
    ! Simulates one draw from a normal distribution using
    !     Box-Muller tranformation.
    !##############################################################################
    subroutine simulate_normal_1(x, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma


        !##### OTHER VARIABLES ####################################################

        real*8 :: uni1, uni2, mu_c, sigma_c
        real*8 :: pi = 3.141592653589793d0


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sigma

        ! simulate a uniform variable draw
        call simulate_uniform_1(uni1)
        call simulate_uniform_1(uni2)

        ! transform by Box-Muller transformation
        x = sqrt(-2d0*log(uni1))*cos(2d0*pi*uni2)

        ! transform to mean and variance
        x = mu_c + sqrt(sigma_c)*x

    end subroutine simulate_normal_1


    !##############################################################################
    ! SUBROUTINE simulate_normal_n
    !
    ! Simulates draws from a normal distribution using
    !     Box-Muller tranformation.
    !##############################################################################
    subroutine simulate_normal_n(x, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x(:)

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma


        !##### OTHER VARIABLES ####################################################

        real*8 :: uni1(size(x, 1)/2), uni2(size(x, 1)/2), mu_c, sigma_c
        integer :: n, in, n2
        real*8 :: pi = 3.141592653589793d0


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sigma

        ! get size of x
        n = size(x, 1)
        n2 = n/2

        ! simulate a uniform variable draw
        call simulate_uniform_n(uni1)
        call simulate_uniform_n(uni2)

        ! transform by Box-Muller transformation
        do in = 1, n2
            x(2*(in-1)+1) = sqrt(-2d0*log(uni1(in)))*cos(2d0*pi*uni2(in))
            x(2*(in-1)+2) = sqrt(-2d0*log(uni1(in)))*sin(2d0*pi*uni2(in))
        enddo

        if(2*n2 /= n)then
            call simulate_normal_1(x(n), 0d0, 1d0)
        endif

        ! transform to mean and variance
        x = mu_c + sqrt(sigma_c)*x

    end subroutine simulate_normal_n


    !##############################################################################
    ! SUBROUTINE simulate_lognormal_1
    !
    ! Simulates one draw from a log normal distribution.
    !##############################################################################
    subroutine simulate_log_normal_1(x, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma


        !##### OTHER VARIABLES ####################################################

        real*8 :: mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 1d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sigma

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! simulate normal and convert to log_normal
        call simulate_normal_1(x, mu_c, sigma_c)

        ! transform to log normal
        x = exp(x)

    end subroutine simulate_log_normal_1


    !##############################################################################
    ! SUBROUTINE simulate_log_normal_n
    !
    ! Simulates draws from a log normal distribution.
    !##############################################################################
    subroutine simulate_log_normal_n(x, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x(:)

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma


        !##### OTHER VARIABLES ####################################################

        real*8 :: mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 1d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sigma

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! simulate normal and convert to log_normal
        call simulate_normal_n(x, mu_c, sigma_c)

        ! transform to log normal
        x = exp(x)

    end subroutine simulate_log_normal_n


    !##############################################################################
    ! SUBROUTINE simulate_Gamma_1
    !
    ! Simulates draws from a gamma distribution using the
    !     Marsaglia/Tsang method.
    !##############################################################################
    subroutine simulate_Gamma_1(x, alpha, beta)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x

        ! shape parameter
        real*8, optional :: alpha

        ! rate parameter
        real*8, optional :: beta


        !##### OTHER VARIABLES ####################################################

        real*8 :: alpha_c, beta_c, d, c, uni


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        alpha_c = 1d0
        if(present(alpha))alpha_c = alpha
        beta_c = 1d0
        if(present(beta))beta_c = beta

        ! check for validity of beta
        if(beta_c <= 0d0)then
            call error('simulate_Gamma','beta has a non-positive value')
        endif

        ! check whether alpha >= 1
        if(alpha_c >= 1d0)then

            ! calculate c and d
            d = alpha_c - 1d0/3d0
            c = 1d0/sqrt(9d0*d)

            ! simulate the gammas step by step
            x = simulate_Gamma_plain(c, d)
        else

            ! generate alpha+1 variables
            d = alpha_c + 2d0/3d0
            c = 1d0/sqrt(9d0*d)

            ! simulate the gammas step by step
            x = simulate_Gamma_plain(c, d)

            ! add uniform transformation
            call simulate_uniform_1(uni)
            x = x*uni**(1d0/alpha_c)

        endif

        ! apply scaling parameter
        x = x/beta_c


    end subroutine simulate_Gamma_1


    !##############################################################################
    ! SUBROUTINE simulate_Gamma_n
    !
    ! Simulates draws from a gamma distribution using the
    !     Marsaglia/Tsang method.
    !##############################################################################
    subroutine simulate_Gamma_n(x, alpha, beta)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x(:)

        ! shape parameter
        real*8, optional :: alpha

        ! rate parameter
        real*8, optional :: beta


        !##### OTHER VARIABLES ####################################################

        real*8 :: alpha_c, beta_c, d, c, uni(size(x, 1))
        integer :: n, in


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        alpha_c = 1d0
        if(present(alpha))alpha_c = alpha
        beta_c = 1d0
        if(present(beta))beta_c = beta

        ! check for validity of beta
        if(beta_c <= 0d0)then
            call error('simulate_Gamma','beta has a non-positive value')
        endif

        ! get size of x
        n = size(x, 1)

        ! check whether alpha >= 1
        if(alpha_c >= 1d0)then

            ! calculate c and d
            d = alpha_c - 1d0/3d0
            c = 1d0/sqrt(9d0*d)

            ! simulate the gammas step by step
            do in = 1, n
                x(in) = simulate_Gamma_plain(c, d)
            enddo
        else

            ! generate alpha+1 variables
            d = alpha_c + 2d0/3d0
            c = 1d0/sqrt(9d0*d)

            ! simulate the gammas step by step
            do in = 1, n
                x(in) = simulate_Gamma_plain(c, d)
            enddo

            ! add uniform transformation
            call simulate_uniform_n(uni)
            x = x*uni**(1d0/alpha_c)

        endif

        ! apply scaling parameter
        x = x/beta_c


    end subroutine simulate_Gamma_n


    !##############################################################################
    ! FUNCTION simulate_Gamma_plain
    !
    ! Simulates one draw from gamma distribution without any testing.
    !##############################################################################
    function simulate_Gamma_plain(c, d) result(res)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! c parameter
        real*8 :: c

        ! d parameter
        real*8 :: d

        ! a gamma realization
        real*8 :: res


        !##### OTHER VARIABLES ####################################################

        real*8 :: v, x, u


        !##### ROUTINE CODE #######################################################

        do
            ! generate v
            call simulate_normal_1(x)
            v = (1d0+c*x)
            v = v*v*v

            ! test v
            if(v <= 0d0)cycle

            ! generate u
            call simulate_uniform_1(u)

            ! check squeeze
            if(u < 1d0 - 0.0331d0*(x*x)*(x*x))then
                res = d*v
                return
            endif

            ! check real condition
            if(log(u) < 0.5d0*x*x + d*(1d0 - v + log(v)))then
                res = d*v
                return
            endif
        enddo

    end function simulate_Gamma_plain


    !##############################################################################
    ! SUBROUTINE simulate_beta_1
    !
    ! Simulates one draw from a beta distribution.
    !##############################################################################
    subroutine simulate_beta_1(x, p, q)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x

        ! shape parameter
        real*8, optional :: p

        ! rate parameter
        real*8, optional :: q


        !##### OTHER VARIABLES ####################################################

        real*8 :: p_c, q_c, gamma1, gamma2, bern


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        p_c = 1d0
        if(present(p))p_c = p
        q_c = 1d0
        if(present(q))q_c = q

        ! check for validity of parameters
        if(p_c <= 0d0)then
            call error('simulate_beta','p has a non-positive value')
        endif
        if(q_c <= 0d0)then
            call error('simulate_beta','q has a non-positive value')
        endif

        ! get bernoulli parametzer for small p and q
        bern = p/(p+q)

        call simulate_Gamma_1(gamma1, p)
        call simulate_Gamma_1(gamma2, q)

        ! check whether gammas both greaters 0
        if(gamma1 > 0d0 .or. gamma2 > 0d0)then
            x = gamma1/(gamma1+gamma2)
        else
            x = simulate_bernoulli_plain(bern)
        endif

    end subroutine simulate_beta_1


    !##############################################################################
    ! SUBROUTINE simulate_beta_n
    !
    ! Simulates draws from a beta distribution.
    !##############################################################################
    subroutine simulate_beta_n(x, p, q)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! point into which the draw should be saved
        real*8, intent(out) :: x(:)

        ! shape parameter
        real*8, optional :: p

        ! rate parameter
        real*8, optional :: q


        !##### OTHER VARIABLES ####################################################

        real*8 :: p_c, q_c, gamma1(size(x,1)), gamma2(size(x,1)), bern
        integer :: n, in


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        p_c = 1d0
        if(present(p))p_c = p
        q_c = 1d0
        if(present(q))q_c = q

        ! check for validity of parameters
        if(p_c <= 0d0)then
            call error('simulate_beta','p has a non-positive value')
        endif
        if(q_c <= 0d0)then
            call error('simulate_beta','q has a non-positive value')
        endif

        ! get size of x
        n = size(x, 1)

        ! get bernoulli parametzer for small p and q
        bern = p/(p+q)

        call simulate_Gamma_n(gamma1, p)
        call simulate_Gamma_n(gamma2, q)

        do in = 1, n

            ! check whether gammas both greaters 0
            if(gamma1(in) > 0d0 .or. gamma2(in) > 0d0)then
                x(in) = gamma1(in)/(gamma1(in)+gamma2(in))
            else
                x(in) = simulate_bernoulli_plain(bern)
            endif
        enddo

    end subroutine simulate_beta_n


    !##############################################################################
    ! FUNCTION simulate_bernoulli_plain
    !
    ! Simulates one draw from bernoulli distribution
    !##############################################################################
    function simulate_bernoulli_plain(p) result(res)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! bernoulli parameter
        real*8 :: p

        ! a bernoulli realization
        real*8 :: res


        !##### OTHER VARIABLES ####################################################

        real*8 :: x


        !##### ROUTINE CODE #######################################################

        call simulate_uniform_1(x)

        if(x <= p)then
            res = 1d0
        else
            res = 0d0
        endif

    end function simulate_bernoulli_plain



    !##############################################################################
    ! SUBROUTINE normal_discrete_1
    !
    ! Creates n points and probabilities for a normal distribution.
    !
    ! Taken from Miranda and Fackler's CompEcon Toolkit
    !##############################################################################
    subroutine normal_discrete_1(x, prob, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! discrete points of normal distribution
        real*8, intent(out) :: x(:)

        ! probability weights
        real*8, intent(out) :: prob(:)

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma


        !##### OTHER VARIABLES ####################################################

        real*8 :: mu_c, sigma_c, pim4, z, z1, p1, p2, p3, pp
        integer :: n, m, i, j, its
        integer, parameter :: maxit = 200
        real*8, parameter :: pi = 3.1415926535897932d0


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sqrt(sigma)

        ! check for right array sizes
        n = assert_eq(size(x,1), size(prob,1), 'normal_discrete')

        ! calculate 1/pi^0.25
        pim4 = 1d0/pi**0.25d0

        ! get number of points
        m = (n+1)/2

        ! initialize x and prob
        x = 0d0
        prob = 0d0

        ! start iteration
        do i = 1, m

            ! set reasonable starting values
            if(i == 1)then
                z = sqrt(dble(2*n+1))-1.85575d0*(dble(2*n+1)**(-1d0/6d0))
            elseif(i == 2)then
                z = z - 1.14d0*(dble(n)**0.426d0)/z
            elseif(i == 3)then
                z = 1.86d0*z+0.86d0*x(1)
            elseif(i == 4)then
                z = 1.91d0*z+0.91d0*x(2);
            else
                z = 2d0*z+x(i-2);
            endif

            ! root finding iterations
            its = 0
            do while(its < maxit)
                its = its+1
                p1 = pim4
                p2 = 0d0
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = z*sqrt(2d0/dble(j))*p2-sqrt(dble(j-1)/dble(j))*p3
                enddo
                pp = sqrt(2d0*dble(n))*p2
                z1 = z
                z  = z1-p1/pp
                if(abs(z-z1) < 1e-14)exit
            enddo
            if(its >= maxit)then
                call error('normal_discrete', &
                    'Could not discretize normal distribution')
            endif
            x(n+1-i) = z
            x(i) = -z
            prob(i) = 2d0/pp**2
            prob(n+1-i) = prob(i)
        enddo

        ! set output data
        prob = prob/sqrt(pi)
        x = x*sqrt(2d0)*sigma_c + mu_c

    end subroutine normal_discrete_1


    !##############################################################################
    ! SUBROUTINE normal_discrete_2
    !
    ! Creates n1*n2 points and probabilities for a two-dimensional normal
    !     distribution.
    !
    ! Taken from Miranda and Fackler's CompEcon Toolkit
    !##############################################################################
    subroutine normal_discrete_2(n, x, prob, mu, sigma, rho)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! number of points in every direction
        integer, intent(in) :: n(2)

        ! discrete points of normal distribution
        real*8, intent(out) :: x(:, :)

        ! probability weights
        real*8, intent(out) :: prob(:)

        ! expectation of distribution
        real*8, optional :: mu(2)

        ! variance of distribution
        real*8, optional :: sigma(2)

        ! correlation of distribution
        real*8, optional :: rho


        !##### OTHER VARIABLES ####################################################

        real*8 :: mu_c(2), sig_c(2), rho_c, sigma_c(2,2), l(2,2)
        real*8 :: x1(n(1)), x2(n(2)), p1(n(1)), p2(n(2))
        integer :: m, j, k


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 0d0
        if(present(mu))mu_c = mu
        sig_c(1) = 1d0
        if(present(sigma))sig_c = sigma
        rho_c = 0d0
        if(present(rho))rho_c = rho

        ! set up variance covariance matrix
        sigma_c(1, 1) = sig_c(1)
        sigma_c(2, 2) = sig_c(2)
        sigma_c(1, 2) = rho_c*sqrt(sig_c(1)*sig_c(2))
        sigma_c(2, 1) = sigma_c(1, 2)

        ! check for right array sizes
        m = assert_eq(size(x,1), size(prob,1), n(1)*n(2), 'normal_discrete')
        m = assert_eq(size(x,2), 2, 'normal_discrete')

        ! check whether sigma is symmetric
        if(any(abs(transpose(sigma_c) - sigma_c) > 1d-20)) &
            call error('normal_discrete', &
            'Variance-Covariance matrix is not symmetric')

        ! get standard normal distributed random variables
        call normal_discrete(x1, p1, 0d0, 1d0)
        call normal_discrete(x2, p2, 0d0, 1d0)

        ! get joint distribution
        m = 1
        do k = 1, n(2)
            do j = 1, n(1)
                prob(m) = p1(j)*p2(k)
                x(m, :) = (/x1(j), x2(k)/)
                m = m+1
            enddo
        enddo

        ! decompose var-cov matrix
        if(.not.any(abs(sig_c)  <= 1d-100))then
            call cholesky(sigma_c, l)
        else
            l = 0d0
            l(1,1) = sqrt(sig_c(1))
            l(2,2) = sqrt(sig_c(2))
        endif

        ! calculate distribution
        x = matmul(x, transpose(l))
        x(:, 1) = x(:, 1) + mu_c(1)
        x(:, 2) = x(:, 2) + mu_c(2)

    end subroutine normal_discrete_2


    !##############################################################################
    ! SUBROUTINE log_normal_discrete_1
    !
    ! Creates n points and probabilities for a log-normal distribution.
    !
    ! Taken from Miranda and Fackler's CompEcon Toolkit
    !##############################################################################
    subroutine log_normal_discrete_1(x, prob, mu, sigma)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! discrete points of normal distribution
        real*8, intent(out) :: x(:)

        ! probability weights
        real*8, intent(out) :: prob(:)

        ! expectation of distribution
        real*8, optional :: mu

        ! variance of distribution
        real*8, optional :: sigma


        !##### OTHER VARIABLES ####################################################

        real*8 :: mu_c, sigma_c
        integer :: n

        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 1d0
        if(present(mu))mu_c = mu
        sigma_c = 1d0
        if(present(sigma))sigma_c = sigma

        ! get expectation and variance
        sigma_c = log(1d0+sigma_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sigma_c

        ! check for right array sizes
        n = assert_eq(size(x,1), size(prob,1), 'normal_discrete')
        n = n

        call normal_discrete(x, prob, mu_c, sigma_c)

        x = exp(x)

    end subroutine log_normal_discrete_1


    !##############################################################################
    ! SUBROUTINE log_normal_discrete_2
    !
    ! Creates n1*n2 points and probabilities for a two-dimensional log-normal
    !     distribution.
    !
    ! Taken from Miranda and Fackler's CompEcon Toolkit
    !##############################################################################
    subroutine log_normal_discrete_2(n, x, prob, mu, sigma, rho)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! number of points in every direction
        integer, intent(in) :: n(2)

        ! discrete points of normal distribution
        real*8, intent(out) :: x(:, :)

        ! probability weights
        real*8, intent(out) :: prob(:)

        ! expectation of distribution
        real*8, optional :: mu(2)

        ! variance of distribution
        real*8, optional :: sigma(2)

        ! correlation of distribution
        real*8, optional :: rho


        !##### OTHER VARIABLES ####################################################

        real*8 :: mu_c(2), sig_c(2), rho_c, sigma_c(2,2), l(2,2)
        real*8 :: x1(n(1)), x2(n(2)), p1(n(1)), p2(n(2))
        integer :: m, j, k


        !##### ROUTINE CODE #######################################################

        ! initialize expectation and variance
        mu_c = 1d0
        if(present(mu))mu_c = mu
        sig_c(:) = 1d0
        if(present(sigma))sig_c = sigma
        rho_c = 0d0
        if(present(rho))rho_c = rho

        ! get expectation and variances
        sig_c = log(1d0+sig_c/mu_c**2)
        mu_c  = log(mu_c)-0.5d0*sig_c

        ! set up covariance matrix
        sigma_c(1, 1) = sig_c(1)
        sigma_c(2, 2) = sig_c(2)
        sigma_c(1, 2) = log(rho_c*sqrt(exp(sig_c(1))-1d0)* &
            sqrt(exp(sig_c(2))-1d0)+1d0)
        sigma_c(2, 1) = sigma_c(1, 2)

        ! check for right array sizes
        m = assert_eq(size(x,1), size(prob,1), n(1)*n(2), 'normal_discrete')
        m = assert_eq(size(x,2), 2, 'normal_discrete')

        ! check whether sigma is symmetric
        if(any(abs(transpose(sigma_c) - sigma_c) > 1d-20)) &
            call error('normal_discrete', &
            'Variance-Covariance matrix is not symmetric')

        ! get standard normal distributed random variables
        call normal_discrete(x1, p1, 0d0, 1d0)
        call normal_discrete(x2, p2, 0d0, 1d0)

        ! get joint distribution
        m = 1
        do k = 1, n(2)
            do j = 1, n(1)
                prob(m) = p1(j)*p2(k)
                x(m, :) = (/x1(j), x2(k)/)
                m = m+1
            enddo
        enddo

        ! decompose var-cov matrix
        if(.not.any(abs(sig_c) <= 1d-100))then
            call cholesky(sigma_c, l)
        else
            l = 0d0
            l(1,1) = sqrt(sig_c(1))
            l(2,2) = sqrt(sig_c(2))
        endif

        ! calculate distribution
        x = matmul(x, transpose(l))
        x(:, 1) = x(:, 1) + mu_c(1)
        x(:, 2) = x(:, 2) + mu_c(2)
        x = exp(x)

    end subroutine log_normal_discrete_2


    !##############################################################################
    ! SUBROUTINE loadNormal
    !
    ! Loads coefficients of standard normal distribution. No Inputs no Outputs.
    !##############################################################################
    subroutine loadNormal()

        implicit none


        ! coefficients for normal distribution
        tbox_ncoeffs(1  ) = 0.08266842952037d0
        tbox_ncoeffs(2  ) = 0.08333333333666d0
        tbox_ncoeffs(3  ) = 0.08399823713300d0
        tbox_ncoeffs(4  ) = 0.08466307444597d0
        tbox_ncoeffs(5  ) = 0.08532777880002d0
        tbox_ncoeffs(6  ) = 0.08599228376807d0
        tbox_ncoeffs(7  ) = 0.08665652298054d0
        tbox_ncoeffs(8  ) = 0.08732043014815d0
        tbox_ncoeffs(9  ) = 0.08798393908098d0
        tbox_ncoeffs(10 ) = 0.08864698370847d0
        tbox_ncoeffs(11 ) = 0.08930949809913d0
        tbox_ncoeffs(12 ) = 0.08997141648018d0
        tbox_ncoeffs(13 ) = 0.09063267325717d0
        tbox_ncoeffs(14 ) = 0.09129320303345d0
        tbox_ncoeffs(15 ) = 0.09195294062960d0
        tbox_ncoeffs(16 ) = 0.09261182110270d0
        tbox_ncoeffs(17 ) = 0.09326977976551d0
        tbox_ncoeffs(18 ) = 0.09392675220552d0
        tbox_ncoeffs(19 ) = 0.09458267430386d0
        tbox_ncoeffs(20 ) = 0.09523748225408d0
        tbox_ncoeffs(21 ) = 0.09589111258072d0
        tbox_ncoeffs(22 ) = 0.09654350215782d0
        tbox_ncoeffs(23 ) = 0.09719458822711d0
        tbox_ncoeffs(24 ) = 0.09784430841619d0
        tbox_ncoeffs(25 ) = 0.09849260075636d0
        tbox_ncoeffs(26 ) = 0.09913940370038d0
        tbox_ncoeffs(27 ) = 0.09978465613992d0
        tbox_ncoeffs(28 ) = 0.10042829742287d0
        tbox_ncoeffs(29 ) = 0.10107026737037d0
        tbox_ncoeffs(30 ) = 0.10171050629368d0
        tbox_ncoeffs(31 ) = 0.10234895501070d0
        tbox_ncoeffs(32 ) = 0.10298555486238d0
        tbox_ncoeffs(33 ) = 0.10362024772873d0
        tbox_ncoeffs(34 ) = 0.10425297604471d0
        tbox_ncoeffs(35 ) = 0.10488368281574d0
        tbox_ncoeffs(36 ) = 0.10551231163298d0
        tbox_ncoeffs(37 ) = 0.10613880668836d0
        tbox_ncoeffs(38 ) = 0.10676311278921d0
        tbox_ncoeffs(39 ) = 0.10738517537278d0
        tbox_ncoeffs(40 ) = 0.10800494052022d0
        tbox_ncoeffs(41 ) = 0.10862235497049d0
        tbox_ncoeffs(42 ) = 0.10923736613378d0
        tbox_ncoeffs(43 ) = 0.10984992210470d0
        tbox_ncoeffs(44 ) = 0.11045997167510d0
        tbox_ncoeffs(45 ) = 0.11106746434664d0
        tbox_ncoeffs(46 ) = 0.11167235034289d0
        tbox_ncoeffs(47 ) = 0.11227458062122d0
        tbox_ncoeffs(48 ) = 0.11287410688430d0
        tbox_ncoeffs(49 ) = 0.11347088159122d0
        tbox_ncoeffs(50 ) = 0.11406485796828d0
        tbox_ncoeffs(51 ) = 0.11465599001946d0
        tbox_ncoeffs(52 ) = 0.11524423253650d0
        tbox_ncoeffs(53 ) = 0.11582954110857d0
        tbox_ncoeffs(54 ) = 0.11641187213170d0
        tbox_ncoeffs(55 ) = 0.11699118281768d0
        tbox_ncoeffs(56 ) = 0.11756743120273d0
        tbox_ncoeffs(57 ) = 0.11814057615572d0
        tbox_ncoeffs(58 ) = 0.11871057738604d0
        tbox_ncoeffs(59 ) = 0.11927739545108d0
        tbox_ncoeffs(60 ) = 0.11984099176332d0
        tbox_ncoeffs(61 ) = 0.12040132859708d0
        tbox_ncoeffs(62 ) = 0.12095836909488d0
        tbox_ncoeffs(63 ) = 0.12151207727333d0
        tbox_ncoeffs(64 ) = 0.12206241802879d0
        tbox_ncoeffs(65 ) = 0.12260935714252d0
        tbox_ncoeffs(66 ) = 0.12315286128547d0
        tbox_ncoeffs(67 ) = 0.12369289802275d0
        tbox_ncoeffs(68 ) = 0.12422943581766d0
        tbox_ncoeffs(69 ) = 0.12476244403529d0
        tbox_ncoeffs(70 ) = 0.12529189294588d0
        tbox_ncoeffs(71 ) = 0.12581775372763d0
        tbox_ncoeffs(72 ) = 0.12633999846928d0
        tbox_ncoeffs(73 ) = 0.12685860017217d0
        tbox_ncoeffs(74 ) = 0.12737353275205d0
        tbox_ncoeffs(75 ) = 0.12788477104040d0
        tbox_ncoeffs(76 ) = 0.12839229078547d0
        tbox_ncoeffs(77 ) = 0.12889606865292d0
        tbox_ncoeffs(78 ) = 0.12939608222599d0
        tbox_ncoeffs(79 ) = 0.12989231000551d0
        tbox_ncoeffs(80 ) = 0.13038473140932d0
        tbox_ncoeffs(81 ) = 0.13087332677149d0
        tbox_ncoeffs(82 ) = 0.13135807734110d0
        tbox_ncoeffs(83 ) = 0.13183896528071d0
        tbox_ncoeffs(84 ) = 0.13231597366444d0
        tbox_ncoeffs(85 ) = 0.13278908647571d0
        tbox_ncoeffs(86 ) = 0.13325828860466d0
        tbox_ncoeffs(87 ) = 0.13372356584521d0
        tbox_ncoeffs(88 ) = 0.13418490489180d0
        tbox_ncoeffs(89 ) = 0.13464229333577d0
        tbox_ncoeffs(90 ) = 0.13509571966143d0
        tbox_ncoeffs(91 ) = 0.13554517324182d0
        tbox_ncoeffs(92 ) = 0.13599064433414d0
        tbox_ncoeffs(93 ) = 0.13643212407487d0
        tbox_ncoeffs(94 ) = 0.13686960447458d0
        tbox_ncoeffs(95 ) = 0.13730307841244d0
        tbox_ncoeffs(96 ) = 0.13773253963042d0
        tbox_ncoeffs(97 ) = 0.13815798272726d0
        tbox_ncoeffs(98 ) = 0.13857940315206d0
        tbox_ncoeffs(99 ) = 0.13899679719767d0
        tbox_ncoeffs(100) = 0.13941016199375d0
        tbox_ncoeffs(101) = 0.13981949549964d0
        tbox_ncoeffs(102) = 0.14022479649686d0
        tbox_ncoeffs(103) = 0.14062606458146d0
        tbox_ncoeffs(104) = 0.14102330015605d0
        tbox_ncoeffs(105) = 0.14141650442162d0
        tbox_ncoeffs(106) = 0.14180567936913d0
        tbox_ncoeffs(107) = 0.14219082777086d0
        tbox_ncoeffs(108) = 0.14257195317152d0
        tbox_ncoeffs(109) = 0.14294905987916d0
        tbox_ncoeffs(110) = 0.14332215295591d0
        tbox_ncoeffs(111) = 0.14369123820844d0
        tbox_ncoeffs(112) = 0.14405632217830d0
        tbox_ncoeffs(113) = 0.14441741213199d0
        tbox_ncoeffs(114) = 0.14477451605098d0
        tbox_ncoeffs(115) = 0.14512764262136d0
        tbox_ncoeffs(116) = 0.14547680122355d0
        tbox_ncoeffs(117) = 0.14582200192164d0
        tbox_ncoeffs(118) = 0.14616325545274d0
        tbox_ncoeffs(119) = 0.14650057321607d0
        tbox_ncoeffs(120) = 0.14683396726198d0
        tbox_ncoeffs(121) = 0.14716345028082d0
        tbox_ncoeffs(122) = 0.14748903559165d0
        tbox_ncoeffs(123) = 0.14781073713089d0
        tbox_ncoeffs(124) = 0.14812856944082d0
        tbox_ncoeffs(125) = 0.14844254765799d0
        tbox_ncoeffs(126) = 0.14875268750152d0
        tbox_ncoeffs(127) = 0.14905900526133d0
        tbox_ncoeffs(128) = 0.14936151778628d0
        tbox_ncoeffs(129) = 0.14966024247223d0
        tbox_ncoeffs(130) = 0.14995519725000d0
        tbox_ncoeffs(131) = 0.15024640057335d0
        tbox_ncoeffs(132) = 0.15053387140685d0
        tbox_ncoeffs(133) = 0.15081762921365d0
        tbox_ncoeffs(134) = 0.15109769394332d0
        tbox_ncoeffs(135) = 0.15137408601959d0
        tbox_ncoeffs(136) = 0.15164682632804d0
        tbox_ncoeffs(137) = 0.15191593620381d0
        tbox_ncoeffs(138) = 0.15218143741931d0
        tbox_ncoeffs(139) = 0.15244335217185d0
        tbox_ncoeffs(140) = 0.15270170307132d0
        tbox_ncoeffs(141) = 0.15295651312786d0
        tbox_ncoeffs(142) = 0.15320780573956d0
        tbox_ncoeffs(143) = 0.15345560468012d0
        tbox_ncoeffs(144) = 0.15369993408657d0
        tbox_ncoeffs(145) = 0.15394081844705d0
        tbox_ncoeffs(146) = 0.15417828258851d0
        tbox_ncoeffs(147) = 0.15441235166460d0
        tbox_ncoeffs(148) = 0.15464305114346d0
        tbox_ncoeffs(149) = 0.15487040679567d0
        tbox_ncoeffs(150) = 0.15509444468217d0
        tbox_ncoeffs(151) = 0.15531519114231d0
        tbox_ncoeffs(152) = 0.15553267278188d0
        tbox_ncoeffs(153) = 0.15574691646132d0
        tbox_ncoeffs(154) = 0.15595794928392d0
        tbox_ncoeffs(155) = 0.15616579858408d0
        tbox_ncoeffs(156) = 0.15637049191578d0
        tbox_ncoeffs(157) = 0.15657205704098d0
        tbox_ncoeffs(158) = 0.15677052191823d0
        tbox_ncoeffs(159) = 0.15696591469131d0
        tbox_ncoeffs(160) = 0.15715826367800d0
        tbox_ncoeffs(161) = 0.15734759735896d0
        tbox_ncoeffs(162) = 0.15753394436670d0
        tbox_ncoeffs(163) = 0.15771733347468d0
        tbox_ncoeffs(164) = 0.15789779358648d0
        tbox_ncoeffs(165) = 0.15807535372517d0
        tbox_ncoeffs(166) = 0.15825004302275d0
        tbox_ncoeffs(167) = 0.15842189070971d0
        tbox_ncoeffs(168) = 0.15859092610475d0
        tbox_ncoeffs(169) = 0.15875717860458d0
        tbox_ncoeffs(170) = 0.15892067767398d0
        tbox_ncoeffs(171) = 0.15908145283580d0
        tbox_ncoeffs(172) = 0.15923953366129d0
        tbox_ncoeffs(173) = 0.15939494976049d0
        tbox_ncoeffs(174) = 0.15954773077272d0
        tbox_ncoeffs(175) = 0.15969790635731d0
        tbox_ncoeffs(176) = 0.15984550618444d0
        tbox_ncoeffs(177) = 0.15999055992611d0
        tbox_ncoeffs(178) = 0.16013309724729d0
        tbox_ncoeffs(179) = 0.16027314779723d0
        tbox_ncoeffs(180) = 0.16041074120091d0
        tbox_ncoeffs(181) = 0.16054590705063d0
        tbox_ncoeffs(182) = 0.16067867489785d0
        tbox_ncoeffs(183) = 0.16080907424505d0
        tbox_ncoeffs(184) = 0.16093713453790d0
        tbox_ncoeffs(185) = 0.16106288515746d0
        tbox_ncoeffs(186) = 0.16118635541264d0
        tbox_ncoeffs(187) = 0.16130757453281d0
        tbox_ncoeffs(188) = 0.16142657166051d0
        tbox_ncoeffs(189) = 0.16154337584441d0
        tbox_ncoeffs(190) = 0.16165801603239d0
        tbox_ncoeffs(191) = 0.16177052106482d0
        tbox_ncoeffs(192) = 0.16188091966793d0
        tbox_ncoeffs(193) = 0.16198924044747d0
        tbox_ncoeffs(194) = 0.16209551188243d0
        tbox_ncoeffs(195) = 0.16219976231898d0
        tbox_ncoeffs(196) = 0.16230201996458d0
        tbox_ncoeffs(197) = 0.16240231288223d0
        tbox_ncoeffs(198) = 0.16250066898487d0
        tbox_ncoeffs(199) = 0.16259711603006d0
        tbox_ncoeffs(200) = 0.16269168161466d0
        tbox_ncoeffs(201) = 0.16278439316979d0
        tbox_ncoeffs(202) = 0.16287527795595d0
        tbox_ncoeffs(203) = 0.16296436305822d0
        tbox_ncoeffs(204) = 0.16305167538173d0
        tbox_ncoeffs(205) = 0.16313724164721d0
        tbox_ncoeffs(206) = 0.16322108838675d0
        tbox_ncoeffs(207) = 0.16330324193970d0
        tbox_ncoeffs(208) = 0.16338372844873d0
        tbox_ncoeffs(209) = 0.16346257385602d0
        tbox_ncoeffs(210) = 0.16353980389969d0
        tbox_ncoeffs(211) = 0.16361544411029d0
        tbox_ncoeffs(212) = 0.16368951980749d0
        tbox_ncoeffs(213) = 0.16376205609692d0
        tbox_ncoeffs(214) = 0.16383307786715d0
        tbox_ncoeffs(215) = 0.16390260978683d0
        tbox_ncoeffs(216) = 0.16397067630194d0
        tbox_ncoeffs(217) = 0.16403730163326d0
        tbox_ncoeffs(218) = 0.16410250977393d0
        tbox_ncoeffs(219) = 0.16416632448711d0
        tbox_ncoeffs(220) = 0.16422876930390d0
        tbox_ncoeffs(221) = 0.16428986752129d0
        tbox_ncoeffs(222) = 0.16434964220027d0
        tbox_ncoeffs(223) = 0.16440811616413d0
        tbox_ncoeffs(224) = 0.16446531199680d0
        tbox_ncoeffs(225) = 0.16452125204142d0
        tbox_ncoeffs(226) = 0.16457595839893d0
        tbox_ncoeffs(227) = 0.16462945292691d0
        tbox_ncoeffs(228) = 0.16468175723840d0
        tbox_ncoeffs(229) = 0.16473289270097d0
        tbox_ncoeffs(230) = 0.16478288043584d0
        tbox_ncoeffs(231) = 0.16483174131713d0
        tbox_ncoeffs(232) = 0.16487949597123d0
        tbox_ncoeffs(233) = 0.16492616477627d0
        tbox_ncoeffs(234) = 0.16497176786173d0
        tbox_ncoeffs(235) = 0.16501632510810d0
        tbox_ncoeffs(236) = 0.16505985614672d0
        tbox_ncoeffs(237) = 0.16510238035967d0
        tbox_ncoeffs(238) = 0.16514391687977d0
        tbox_ncoeffs(239) = 0.16518448459069d0
        tbox_ncoeffs(240) = 0.16522410212715d0
        tbox_ncoeffs(241) = 0.16526278787521d0
        tbox_ncoeffs(242) = 0.16530055997267d0
        tbox_ncoeffs(243) = 0.16533743630952d0
        tbox_ncoeffs(244) = 0.16537343452852d0
        tbox_ncoeffs(245) = 0.16540857202584d0
        tbox_ncoeffs(246) = 0.16544286595180d0
        tbox_ncoeffs(247) = 0.16547633321163d0
        tbox_ncoeffs(248) = 0.16550899046642d0
        tbox_ncoeffs(249) = 0.16554085413405d0
        tbox_ncoeffs(250) = 0.16557194039022d0
        tbox_ncoeffs(251) = 0.16560226516953d0
        tbox_ncoeffs(252) = 0.16563184416672d0
        tbox_ncoeffs(253) = 0.16566069283783d0
        tbox_ncoeffs(254) = 0.16568882640156d0
        tbox_ncoeffs(255) = 0.16571625984060d0
        tbox_ncoeffs(256) = 0.16574300790308d0
        tbox_ncoeffs(257) = 0.16576908510400d0
        tbox_ncoeffs(258) = 0.16579450572684d0
        tbox_ncoeffs(259) = 0.16581928382508d0
        tbox_ncoeffs(260) = 0.16584343322386d0
        tbox_ncoeffs(261) = 0.16586696752170d0
        tbox_ncoeffs(262) = 0.16588990009220d0
        tbox_ncoeffs(263) = 0.16591224408580d0
        tbox_ncoeffs(264) = 0.16593401243165d0
        tbox_ncoeffs(265) = 0.16595521783947d0
        tbox_ncoeffs(266) = 0.16597587280139d0
        tbox_ncoeffs(267) = 0.16599598959395d0
        tbox_ncoeffs(268) = 0.16601558028006d0
        tbox_ncoeffs(269) = 0.16603465671097d0
        tbox_ncoeffs(270) = 0.16605323052837d0
        tbox_ncoeffs(271) = 0.16607131316638d0
        tbox_ncoeffs(272) = 0.16608891585371d0
        tbox_ncoeffs(273) = 0.16610604961574d0
        tbox_ncoeffs(274) = 0.16612272527667d0
        tbox_ncoeffs(275) = 0.16613895346170d0
        tbox_ncoeffs(276) = 0.16615474459919d0
        tbox_ncoeffs(277) = 0.16617010892290d0
        tbox_ncoeffs(278) = 0.16618505647417d0
        tbox_ncoeffs(279) = 0.16619959710420d0
        tbox_ncoeffs(280) = 0.16621374047627d0
        tbox_ncoeffs(281) = 0.16622749606802d0
        tbox_ncoeffs(282) = 0.16624087317374d0
        tbox_ncoeffs(283) = 0.16625388090661d0
        tbox_ncoeffs(284) = 0.16626652820105d0
        tbox_ncoeffs(285) = 0.16627882381500d0
        tbox_ncoeffs(286) = 0.16629077633222d0
        tbox_ncoeffs(287) = 0.16630239416459d0
        tbox_ncoeffs(288) = 0.16631368555449d0
        tbox_ncoeffs(289) = 0.16632465857703d0
        tbox_ncoeffs(290) = 0.16633532114244d0
        tbox_ncoeffs(291) = 0.16634568099833d0
        tbox_ncoeffs(292) = 0.16635574573206d0
        tbox_ncoeffs(293) = 0.16636552277304d0
        tbox_ncoeffs(294) = 0.16637501939499d0
        tbox_ncoeffs(295) = 0.16638424271833d0
        tbox_ncoeffs(296) = 0.16639319971241d0
        tbox_ncoeffs(297) = 0.16640189719786d0
        tbox_ncoeffs(298) = 0.16641034184880d0
        tbox_ncoeffs(299) = 0.16641854019521d0
        tbox_ncoeffs(300) = 0.16642649862513d0
        tbox_ncoeffs(301) = 0.16643422338694d0
        tbox_ncoeffs(302) = 0.16644172059162d0
        tbox_ncoeffs(303) = 0.16644899621496d0
        tbox_ncoeffs(304) = 0.16645605609979d0
        tbox_ncoeffs(305) = 0.16646290595821d0
        tbox_ncoeffs(306) = 0.16646955137377d0
        tbox_ncoeffs(307) = 0.16647599780362d0
        tbox_ncoeffs(308) = 0.16648225058074d0
        tbox_ncoeffs(309) = 0.16648831491603d0
        tbox_ncoeffs(310) = 0.16649419590048d0
        tbox_ncoeffs(311) = 0.16649989850726d0
        tbox_ncoeffs(312) = 0.16650542759386d0
        tbox_ncoeffs(313) = 0.16651078790409d0
        tbox_ncoeffs(314) = 0.16651598407025d0
        tbox_ncoeffs(315) = 0.16652102061508d0
        tbox_ncoeffs(316) = 0.16652590195381d0
        tbox_ncoeffs(317) = 0.16653063239621d0
        tbox_ncoeffs(318) = 0.16653521614851d0
        tbox_ncoeffs(319) = 0.16653965731538d0
        tbox_ncoeffs(320) = 0.16654395990190d0
        tbox_ncoeffs(321) = 0.16654812781546d0
        tbox_ncoeffs(322) = 0.16655216486763d0
        tbox_ncoeffs(323) = 0.16655607477611d0
        tbox_ncoeffs(324) = 0.16655986116650d0
        tbox_ncoeffs(325) = 0.16656352757422d0
        tbox_ncoeffs(326) = 0.16656707744623d0
        tbox_ncoeffs(327) = 0.16657051414291d0
        tbox_ncoeffs(328) = 0.16657384093974d0
        tbox_ncoeffs(329) = 0.16657706102911d0
        tbox_ncoeffs(330) = 0.16658017752200d0
        tbox_ncoeffs(331) = 0.16658319344969d0
        tbox_ncoeffs(332) = 0.16658611176545d0
        tbox_ncoeffs(333) = 0.16658893534614d0
        tbox_ncoeffs(334) = 0.16659166699389d0
        tbox_ncoeffs(335) = 0.16659430943769d0
        tbox_ncoeffs(336) = 0.16659686533495d0
        tbox_ncoeffs(337) = 0.16659933727306d0
        tbox_ncoeffs(338) = 0.16660172777095d0
        tbox_ncoeffs(339) = 0.16660403928057d0
        tbox_ncoeffs(340) = 0.16660627418839d0
        tbox_ncoeffs(341) = 0.16660843481686d0
        tbox_ncoeffs(342) = 0.16661052342585d0
        tbox_ncoeffs(343) = 0.16661254221407d0
        tbox_ncoeffs(344) = 0.16661449332045d0
        tbox_ncoeffs(345) = 0.16661637882554d0
        tbox_ncoeffs(346) = 0.16661820075280d0
        tbox_ncoeffs(347) = 0.16661996106999d0
        tbox_ncoeffs(348) = 0.16662166169043d0
        tbox_ncoeffs(349) = 0.16662330447431d0
        tbox_ncoeffs(350) = 0.16662489122990d0
        tbox_ncoeffs(351) = 0.16662642371482d0
        tbox_ncoeffs(352) = 0.16662790363726d0
        tbox_ncoeffs(353) = 0.16662933265712d0
        tbox_ncoeffs(354) = 0.16663071238726d0
        tbox_ncoeffs(355) = 0.16663204439455d0
        tbox_ncoeffs(356) = 0.16663333020107d0
        tbox_ncoeffs(357) = 0.16663457128517d0
        tbox_ncoeffs(358) = 0.16663576908260d0
        tbox_ncoeffs(359) = 0.16663692498750d0
        tbox_ncoeffs(360) = 0.16663804035350d0
        tbox_ncoeffs(361) = 0.16663911649474d0
        tbox_ncoeffs(362) = 0.16664015468682d0
        tbox_ncoeffs(363) = 0.16664115616783d0
        tbox_ncoeffs(364) = 0.16664212213929d0
        tbox_ncoeffs(365) = 0.16664305376710d0
        tbox_ncoeffs(366) = 0.16664395218244d0
        tbox_ncoeffs(367) = 0.16664481848270d0
        tbox_ncoeffs(368) = 0.16664565373234d0
        tbox_ncoeffs(369) = 0.16664645896379d0
        tbox_ncoeffs(370) = 0.16664723517823d0
        tbox_ncoeffs(371) = 0.16664798334649d0
        tbox_ncoeffs(372) = 0.16664870440983d0
        tbox_ncoeffs(373) = 0.16664939928072d0
        tbox_ncoeffs(374) = 0.16665006884363d0
        tbox_ncoeffs(375) = 0.16665071395579d0
        tbox_ncoeffs(376) = 0.16665133544793d0
        tbox_ncoeffs(377) = 0.16665193412500d0
        tbox_ncoeffs(378) = 0.16665251076687d0
        tbox_ncoeffs(379) = 0.16665306612903d0
        tbox_ncoeffs(380) = 0.16665360094329d0
        tbox_ncoeffs(381) = 0.16665411591840d0
        tbox_ncoeffs(382) = 0.16665461174071d0
        tbox_ncoeffs(383) = 0.16665508907482d0
        tbox_ncoeffs(384) = 0.16665554856414d0
        tbox_ncoeffs(385) = 0.16665599083158d0
        tbox_ncoeffs(386) = 0.16665641648004d0
        tbox_ncoeffs(387) = 0.16665682609306d0
        tbox_ncoeffs(388) = 0.16665722023531d0
        tbox_ncoeffs(389) = 0.16665759945318d0
        tbox_ncoeffs(390) = 0.16665796427532d0
        tbox_ncoeffs(391) = 0.16665831521312d0
        tbox_ncoeffs(392) = 0.16665865276121d0
        tbox_ncoeffs(393) = 0.16665897739802d0
        tbox_ncoeffs(394) = 0.16665928958617d0
        tbox_ncoeffs(395) = 0.16665958977300d0
        tbox_ncoeffs(396) = 0.16665987839103d0
        tbox_ncoeffs(397) = 0.16666015585833d0
        tbox_ncoeffs(398) = 0.16666042257905d0
        tbox_ncoeffs(399) = 0.16666067894377d0
        tbox_ncoeffs(400) = 0.16666092532995d0
        tbox_ncoeffs(401) = 0.16666116210230d0
        tbox_ncoeffs(402) = 0.16666138961320d0
        tbox_ncoeffs(403) = 0.16666160820305d0
        tbox_ncoeffs(404) = 0.16666181820066d0
        tbox_ncoeffs(405) = 0.16666201992359d0
        tbox_ncoeffs(406) = 0.16666221367853d0
        tbox_ncoeffs(407) = 0.16666239976158d0
        tbox_ncoeffs(408) = 0.16666257845867d0
        tbox_ncoeffs(409) = 0.16666275004578d0
        tbox_ncoeffs(410) = 0.16666291478934d0
        tbox_ncoeffs(411) = 0.16666307294647d0
        tbox_ncoeffs(412) = 0.16666322476532d0
        tbox_ncoeffs(413) = 0.16666337048533d0
        tbox_ncoeffs(414) = 0.16666351033750d0
        tbox_ncoeffs(415) = 0.16666364454471d0
        tbox_ncoeffs(416) = 0.16666377332193d0
        tbox_ncoeffs(417) = 0.16666389687649d0
        tbox_ncoeffs(418) = 0.16666401540836d0
        tbox_ncoeffs(419) = 0.16666412911034d0
        tbox_ncoeffs(420) = 0.16666423816833d0
        tbox_ncoeffs(421) = 0.16666434276155d0
        tbox_ncoeffs(422) = 0.16666444306275d0
        tbox_ncoeffs(423) = 0.16666453923844d0
        tbox_ncoeffs(424) = 0.16666463144908d0
        tbox_ncoeffs(425) = 0.16666471984931d0
        tbox_ncoeffs(426) = 0.16666480458810d0
        tbox_ncoeffs(427) = 0.16666488580898d0
        tbox_ncoeffs(428) = 0.16666496365022d0
        tbox_ncoeffs(429) = 0.16666503824498d0
        tbox_ncoeffs(430) = 0.16666510972151d0
        tbox_ncoeffs(431) = 0.16666517820332d0
        tbox_ncoeffs(432) = 0.16666524380931d0
        tbox_ncoeffs(433) = 0.16666530665397d0
        tbox_ncoeffs(434) = 0.16666536684750d0
        tbox_ncoeffs(435) = 0.16666542449598d0
        tbox_ncoeffs(436) = 0.16666547970149d0
        tbox_ncoeffs(437) = 0.16666553256227d0
        tbox_ncoeffs(438) = 0.16666558317284d0
        tbox_ncoeffs(439) = 0.16666563162416d0
        tbox_ncoeffs(440) = 0.16666567800370d0
        tbox_ncoeffs(441) = 0.16666572239561d0
        tbox_ncoeffs(442) = 0.16666576488083d0
        tbox_ncoeffs(443) = 0.16666580553718d0
        tbox_ncoeffs(444) = 0.16666584443950d0
        tbox_ncoeffs(445) = 0.16666588165975d0
        tbox_ncoeffs(446) = 0.16666591726709d0
        tbox_ncoeffs(447) = 0.16666595132801d0
        tbox_ncoeffs(448) = 0.16666598390641d0
        tbox_ncoeffs(449) = 0.16666601506371d0
        tbox_ncoeffs(450) = 0.16666604485891d0
        tbox_ncoeffs(451) = 0.16666607334871d0
        tbox_ncoeffs(452) = 0.16666610058758d0
        tbox_ncoeffs(453) = 0.16666612662784d0
        tbox_ncoeffs(454) = 0.16666615151975d0
        tbox_ncoeffs(455) = 0.16666617531157d0
        tbox_ncoeffs(456) = 0.16666619804964d0
        tbox_ncoeffs(457) = 0.16666621977846d0
        tbox_ncoeffs(458) = 0.16666624054075d0
        tbox_ncoeffs(459) = 0.16666626037751d0
        tbox_ncoeffs(460) = 0.16666627932812d0
        tbox_ncoeffs(461) = 0.16666629743035d0
        tbox_ncoeffs(462) = 0.16666631472044d0
        tbox_ncoeffs(463) = 0.16666633123319d0
        tbox_ncoeffs(464) = 0.16666634700196d0
        tbox_ncoeffs(465) = 0.16666636205877d0
        tbox_ncoeffs(466) = 0.16666637643432d0
        tbox_ncoeffs(467) = 0.16666639015807d0
        tbox_ncoeffs(468) = 0.16666640325826d0
        tbox_ncoeffs(469) = 0.16666641576197d0
        tbox_ncoeffs(470) = 0.16666642769517d0
        tbox_ncoeffs(471) = 0.16666643908275d0
        tbox_ncoeffs(472) = 0.16666644994857d0
        tbox_ncoeffs(473) = 0.16666646031550d0
        tbox_ncoeffs(474) = 0.16666647020546d0
        tbox_ncoeffs(475) = 0.16666647963944d0
        tbox_ncoeffs(476) = 0.16666648863758d0
        tbox_ncoeffs(477) = 0.16666649721914d0
        tbox_ncoeffs(478) = 0.16666650540260d0
        tbox_ncoeffs(479) = 0.16666651320565d0
        tbox_ncoeffs(480) = 0.16666652064521d0
        tbox_ncoeffs(481) = 0.16666652773753d0
        tbox_ncoeffs(482) = 0.16666653449812d0
        tbox_ncoeffs(483) = 0.16666654094186d0
        tbox_ncoeffs(484) = 0.16666654708298d0
        tbox_ncoeffs(485) = 0.16666655293512d0
        tbox_ncoeffs(486) = 0.16666655851130d0
        tbox_ncoeffs(487) = 0.16666656382403d0
        tbox_ncoeffs(488) = 0.16666656888522d0
        tbox_ncoeffs(489) = 0.16666657370632d0
        tbox_ncoeffs(490) = 0.16666657829825d0
        tbox_ncoeffs(491) = 0.16666658267147d0
        tbox_ncoeffs(492) = 0.16666658683598d0
        tbox_ncoeffs(493) = 0.16666659080134d0
        tbox_ncoeffs(494) = 0.16666659457670d0
        tbox_ncoeffs(495) = 0.16666659817080d0
        tbox_ncoeffs(496) = 0.16666660159201d0
        tbox_ncoeffs(497) = 0.16666660484831d0
        tbox_ncoeffs(498) = 0.16666660794734d0
        tbox_ncoeffs(499) = 0.16666661089641d0
        tbox_ncoeffs(500) = 0.16666661370249d0
        tbox_ncoeffs(501) = 0.16666661637226d0
        tbox_ncoeffs(502) = 0.16666661891201d0
        tbox_ncoeffs(503) = 0.16666662132813d0

        ! coefficients for inverse distribution
        tbox_icoeffs(1  ) = -0.00041777129562d0
        tbox_icoeffs(2  ) = 0.00000000000002d0
        tbox_icoeffs(3  ) = 0.00041777129555d0
        tbox_icoeffs(4  ) = 0.00083554521604d0
        tbox_icoeffs(5  ) = 0.00125332438654d0
        tbox_icoeffs(6  ) = 0.00167111143232d0
        tbox_icoeffs(7  ) = 0.00208890897904d0
        tbox_icoeffs(8  ) = 0.00250671965277d0
        tbox_icoeffs(9  ) = 0.00292454608021d0
        tbox_icoeffs(10 ) = 0.00334239088872d0
        tbox_icoeffs(11 ) = 0.00376025670647d0
        tbox_icoeffs(12 ) = 0.00417814616256d0
        tbox_icoeffs(13 ) = 0.00459606188715d0
        tbox_icoeffs(14 ) = 0.00501400651153d0
        tbox_icoeffs(15 ) = 0.00543198266827d0
        tbox_icoeffs(16 ) = 0.00584999299134d0
        tbox_icoeffs(17 ) = 0.00626804011620d0
        tbox_icoeffs(18 ) = 0.00668612667993d0
        tbox_icoeffs(19 ) = 0.00710425532138d0
        tbox_icoeffs(20 ) = 0.00752242868122d0
        tbox_icoeffs(21 ) = 0.00794064940211d0
        tbox_icoeffs(22 ) = 0.00835892012881d0
        tbox_icoeffs(23 ) = 0.00877724350826d0
        tbox_icoeffs(24 ) = 0.00919562218976d0
        tbox_icoeffs(25 ) = 0.00961405882504d0
        tbox_icoeffs(26 ) = 0.01003255606838d0
        tbox_icoeffs(27 ) = 0.01045111657676d0
        tbox_icoeffs(28 ) = 0.01086974300997d0
        tbox_icoeffs(29 ) = 0.01128843803069d0
        tbox_icoeffs(30 ) = 0.01170720430466d0
        tbox_icoeffs(31 ) = 0.01212604450079d0
        tbox_icoeffs(32 ) = 0.01254496129125d0
        tbox_icoeffs(33 ) = 0.01296395735162d0
        tbox_icoeffs(34 ) = 0.01338303536101d0
        tbox_icoeffs(35 ) = 0.01380219800217d0
        tbox_icoeffs(36 ) = 0.01422144796162d0
        tbox_icoeffs(37 ) = 0.01464078792976d0
        tbox_icoeffs(38 ) = 0.01506022060101d0
        tbox_icoeffs(39 ) = 0.01547974867392d0
        tbox_icoeffs(40 ) = 0.01589937485133d0
        tbox_icoeffs(41 ) = 0.01631910184041d0
        tbox_icoeffs(42 ) = 0.01673893235289d0
        tbox_icoeffs(43 ) = 0.01715886910510d0
        tbox_icoeffs(44 ) = 0.01757891481815d0
        tbox_icoeffs(45 ) = 0.01799907221804d0
        tbox_icoeffs(46 ) = 0.01841934403577d0
        tbox_icoeffs(47 ) = 0.01883973300748d0
        tbox_icoeffs(48 ) = 0.01926024187461d0
        tbox_icoeffs(49 ) = 0.01968087338397d0
        tbox_icoeffs(50 ) = 0.02010163028790d0
        tbox_icoeffs(51 ) = 0.02052251534443d0
        tbox_icoeffs(52 ) = 0.02094353131733d0
        tbox_icoeffs(53 ) = 0.02136468097634d0
        tbox_icoeffs(54 ) = 0.02178596709724d0
        tbox_icoeffs(55 ) = 0.02220739246197d0
        tbox_icoeffs(56 ) = 0.02262895985883d0
        tbox_icoeffs(57 ) = 0.02305067208255d0
        tbox_icoeffs(58 ) = 0.02347253193447d0
        tbox_icoeffs(59 ) = 0.02389454222264d0
        tbox_icoeffs(60 ) = 0.02431670576198d0
        tbox_icoeffs(61 ) = 0.02473902537441d0
        tbox_icoeffs(62 ) = 0.02516150388900d0
        tbox_icoeffs(63 ) = 0.02558414414208d0
        tbox_icoeffs(64 ) = 0.02600694897743d0
        tbox_icoeffs(65 ) = 0.02642992124635d0
        tbox_icoeffs(66 ) = 0.02685306380788d0
        tbox_icoeffs(67 ) = 0.02727637952888d0
        tbox_icoeffs(68 ) = 0.02769987128423d0
        tbox_icoeffs(69 ) = 0.02812354195691d0
        tbox_icoeffs(70 ) = 0.02854739443821d0
        tbox_icoeffs(71 ) = 0.02897143162783d0
        tbox_icoeffs(72 ) = 0.02939565643406d0
        tbox_icoeffs(73 ) = 0.02982007177390d0
        tbox_icoeffs(74 ) = 0.03024468057323d0
        tbox_icoeffs(75 ) = 0.03066948576697d0
        tbox_icoeffs(76 ) = 0.03109449029920d0
        tbox_icoeffs(77 ) = 0.03151969712335d0
        tbox_icoeffs(78 ) = 0.03194510920233d0
        tbox_icoeffs(79 ) = 0.03237072950868d0
        tbox_icoeffs(80 ) = 0.03279656102476d0
        tbox_icoeffs(81 ) = 0.03322260674288d0
        tbox_icoeffs(82 ) = 0.03364886966548d0
        tbox_icoeffs(83 ) = 0.03407535280525d0
        tbox_icoeffs(84 ) = 0.03450205918534d0
        tbox_icoeffs(85 ) = 0.03492899183952d0
        tbox_icoeffs(86 ) = 0.03535615381230d0
        tbox_icoeffs(87 ) = 0.03578354815915d0
        tbox_icoeffs(88 ) = 0.03621117794662d0
        tbox_icoeffs(89 ) = 0.03663904625254d0
        tbox_icoeffs(90 ) = 0.03706715616621d0
        tbox_icoeffs(91 ) = 0.03749551078850d0
        tbox_icoeffs(92 ) = 0.03792411323210d0
        tbox_icoeffs(93 ) = 0.03835296662164d0
        tbox_icoeffs(94 ) = 0.03878207409392d0
        tbox_icoeffs(95 ) = 0.03921143879803d0
        tbox_icoeffs(96 ) = 0.03964106389558d0
        tbox_icoeffs(97 ) = 0.04007095256082d0
        tbox_icoeffs(98 ) = 0.04050110798091d0
        tbox_icoeffs(99 ) = 0.04093153335604d0
        tbox_icoeffs(100) = 0.04136223189960d0
        tbox_icoeffs(101) = 0.04179320683846d0
        tbox_icoeffs(102) = 0.04222446141305d0
        tbox_icoeffs(103) = 0.04265599887763d0
        tbox_icoeffs(104) = 0.04308782250045d0
        tbox_icoeffs(105) = 0.04351993556395d0
        tbox_icoeffs(106) = 0.04395234136498d0
        tbox_icoeffs(107) = 0.04438504321494d0
        tbox_icoeffs(108) = 0.04481804444005d0
        tbox_icoeffs(109) = 0.04525134838153d0
        tbox_icoeffs(110) = 0.04568495839577d0
        tbox_icoeffs(111) = 0.04611887785460d0
        tbox_icoeffs(112) = 0.04655311014544d0
        tbox_icoeffs(113) = 0.04698765867155d0
        tbox_icoeffs(114) = 0.04742252685225d0
        tbox_icoeffs(115) = 0.04785771812308d0
        tbox_icoeffs(116) = 0.04829323593610d0
        tbox_icoeffs(117) = 0.04872908376003d0
        tbox_icoeffs(118) = 0.04916526508055d0
        tbox_icoeffs(119) = 0.04960178340047d0
        tbox_icoeffs(120) = 0.05003864223996d0
        tbox_icoeffs(121) = 0.05047584513682d0
        tbox_icoeffs(122) = 0.05091339564669d0
        tbox_icoeffs(123) = 0.05135129734326d0
        tbox_icoeffs(124) = 0.05178955381855d0
        tbox_icoeffs(125) = 0.05222816868314d0
        tbox_icoeffs(126) = 0.05266714556640d0
        tbox_icoeffs(127) = 0.05310648811675d0
        tbox_icoeffs(128) = 0.05354620000188d0
        tbox_icoeffs(129) = 0.05398628490906d0
        tbox_icoeffs(130) = 0.05442674654534d0
        tbox_icoeffs(131) = 0.05486758863782d0
        tbox_icoeffs(132) = 0.05530881493393d0
        tbox_icoeffs(133) = 0.05575042920169d0
        tbox_icoeffs(134) = 0.05619243522995d0
        tbox_icoeffs(135) = 0.05663483682869d0
        tbox_icoeffs(136) = 0.05707763782928d0
        tbox_icoeffs(137) = 0.05752084208475d0
        tbox_icoeffs(138) = 0.05796445347010d0
        tbox_icoeffs(139) = 0.05840847588254d0
        tbox_icoeffs(140) = 0.05885291324180d0
        tbox_icoeffs(141) = 0.05929776949043d0
        tbox_icoeffs(142) = 0.05974304859406d0
        tbox_icoeffs(143) = 0.06018875454175d0
        tbox_icoeffs(144) = 0.06063489134623d0
        tbox_icoeffs(145) = 0.06108146304424d0
        tbox_icoeffs(146) = 0.06152847369685d0
        tbox_icoeffs(147) = 0.06197592738972d0
        tbox_icoeffs(148) = 0.06242382823348d0
        tbox_icoeffs(149) = 0.06287218036401d0
        tbox_icoeffs(150) = 0.06332098794277d0
        tbox_icoeffs(151) = 0.06377025515714d0
        tbox_icoeffs(152) = 0.06421998622073d0
        tbox_icoeffs(153) = 0.06467018537375d0
        tbox_icoeffs(154) = 0.06512085688333d0
        tbox_icoeffs(155) = 0.06557200504384d0
        tbox_icoeffs(156) = 0.06602363417731d0
        tbox_icoeffs(157) = 0.06647574863371d0
        tbox_icoeffs(158) = 0.06692835279135d0
        tbox_icoeffs(159) = 0.06738145105726d0
        tbox_icoeffs(160) = 0.06783504786749d0
        tbox_icoeffs(161) = 0.06828914768757d0
        tbox_icoeffs(162) = 0.06874375501282d0
        tbox_icoeffs(163) = 0.06919887436877d0
        tbox_icoeffs(164) = 0.06965451031154d0
        tbox_icoeffs(165) = 0.07011066742824d0
        tbox_icoeffs(166) = 0.07056735033735d0
        tbox_icoeffs(167) = 0.07102456368915d0
        tbox_icoeffs(168) = 0.07148231216611d0
        tbox_icoeffs(169) = 0.07194060048333d0
        tbox_icoeffs(170) = 0.07239943338895d0
        tbox_icoeffs(171) = 0.07285881566455d0
        tbox_icoeffs(172) = 0.07331875212564d0
        tbox_icoeffs(173) = 0.07377924762205d0
        tbox_icoeffs(174) = 0.07424030703842d0
        tbox_icoeffs(175) = 0.07470193529460d0
        tbox_icoeffs(176) = 0.07516413734617d0
        tbox_icoeffs(177) = 0.07562691818485d0
        tbox_icoeffs(178) = 0.07609028283901d0
        tbox_icoeffs(179) = 0.07655423637413d0
        tbox_icoeffs(180) = 0.07701878389330d0
        tbox_icoeffs(181) = 0.07748393053769d0
        tbox_icoeffs(182) = 0.07794968148709d0
        tbox_icoeffs(183) = 0.07841604196038d0
        tbox_icoeffs(184) = 0.07888301721605d0
        tbox_icoeffs(185) = 0.07935061255274d0
        tbox_icoeffs(186) = 0.07981883330977d0
        tbox_icoeffs(187) = 0.08028768486765d0
        tbox_icoeffs(188) = 0.08075717264867d0
        tbox_icoeffs(189) = 0.08122730211743d0
        tbox_icoeffs(190) = 0.08169807878137d0
        tbox_icoeffs(191) = 0.08216950819143d0
        tbox_icoeffs(192) = 0.08264159594252d0
        tbox_icoeffs(193) = 0.08311434767422d0
        tbox_icoeffs(194) = 0.08358776907127d0
        tbox_icoeffs(195) = 0.08406186586428d0
        tbox_icoeffs(196) = 0.08453664383026d0
        tbox_icoeffs(197) = 0.08501210879331d0
        tbox_icoeffs(198) = 0.08548826662520d0
        tbox_icoeffs(199) = 0.08596512324608d0
        tbox_icoeffs(200) = 0.08644268462507d0
        tbox_icoeffs(201) = 0.08692095678098d0
        tbox_icoeffs(202) = 0.08739994578294d0
        tbox_icoeffs(203) = 0.08787965775115d0
        tbox_icoeffs(204) = 0.08836009885749d0
        tbox_icoeffs(205) = 0.08884127532633d0
        tbox_icoeffs(206) = 0.08932319343516d0
        tbox_icoeffs(207) = 0.08980585951539d0
        tbox_icoeffs(208) = 0.09028927995306d0
        tbox_icoeffs(209) = 0.09077346118960d0
        tbox_icoeffs(210) = 0.09125840972261d0
        tbox_icoeffs(211) = 0.09174413210662d0
        tbox_icoeffs(212) = 0.09223063495394d0
        tbox_icoeffs(213) = 0.09271792493537d0
        tbox_icoeffs(214) = 0.09320600878113d0
        tbox_icoeffs(215) = 0.09369489328161d0
        tbox_icoeffs(216) = 0.09418458528826d0
        tbox_icoeffs(217) = 0.09467509171444d0
        tbox_icoeffs(218) = 0.09516641953632d0
        tbox_icoeffs(219) = 0.09565857579372d0
        tbox_icoeffs(220) = 0.09615156759108d0
        tbox_icoeffs(221) = 0.09664540209833d0
        tbox_icoeffs(222) = 0.09714008655187d0
        tbox_icoeffs(223) = 0.09763562825550d0
        tbox_icoeffs(224) = 0.09813203458138d0
        tbox_icoeffs(225) = 0.09862931297106d0
        tbox_icoeffs(226) = 0.09912747093643d0
        tbox_icoeffs(227) = 0.09962651606079d0
        tbox_icoeffs(228) = 0.10012645599988d0
        tbox_icoeffs(229) = 0.10062729848290d0
        tbox_icoeffs(230) = 0.10112905131365d0
        tbox_icoeffs(231) = 0.10163172237156d0
        tbox_icoeffs(232) = 0.10213531961288d0
        tbox_icoeffs(233) = 0.10263985107174d0
        tbox_icoeffs(234) = 0.10314532486136d0
        tbox_icoeffs(235) = 0.10365174917520d0
        tbox_icoeffs(236) = 0.10415913228816d0
        tbox_icoeffs(237) = 0.10466748255783d0
        tbox_icoeffs(238) = 0.10517680842570d0
        tbox_icoeffs(239) = 0.10568711841841d0
        tbox_icoeffs(240) = 0.10619842114908d0
        tbox_icoeffs(241) = 0.10671072531863d0
        tbox_icoeffs(242) = 0.10722403971704d0
        tbox_icoeffs(243) = 0.10773837322480d0
        tbox_icoeffs(244) = 0.10825373481423d0
        tbox_icoeffs(245) = 0.10877013355096d0
        tbox_icoeffs(246) = 0.10928757859528d0
        tbox_icoeffs(247) = 0.10980607920370d0
        tbox_icoeffs(248) = 0.11032564473037d0
        tbox_icoeffs(249) = 0.11084628462865d0
        tbox_icoeffs(250) = 0.11136800845264d0
        tbox_icoeffs(251) = 0.11189082585876d0
        tbox_icoeffs(252) = 0.11241474660740d0
        tbox_icoeffs(253) = 0.11293978056449d0
        tbox_icoeffs(254) = 0.11346593770324d0
        tbox_icoeffs(255) = 0.11399322810583d0
        tbox_icoeffs(256) = 0.11452166196514d0
        tbox_icoeffs(257) = 0.11505124958655d0
        tbox_icoeffs(258) = 0.11558200138971d0
        tbox_icoeffs(259) = 0.11611392791045d0
        tbox_icoeffs(260) = 0.11664703980258d0
        tbox_icoeffs(261) = 0.11718134783991d0
        tbox_icoeffs(262) = 0.11771686291814d0
        tbox_icoeffs(263) = 0.11825359605686d0
        tbox_icoeffs(264) = 0.11879155840166d0
        tbox_icoeffs(265) = 0.11933076122613d0
        tbox_icoeffs(266) = 0.11987121593403d0
        tbox_icoeffs(267) = 0.12041293406146d0
        tbox_icoeffs(268) = 0.12095592727905d0
        tbox_icoeffs(269) = 0.12150020739426d0
        tbox_icoeffs(270) = 0.12204578635364d0
        tbox_icoeffs(271) = 0.12259267624522d0
        tbox_icoeffs(272) = 0.12314088930091d0
        tbox_icoeffs(273) = 0.12369043789893d0
        tbox_icoeffs(274) = 0.12424133456635d0
        tbox_icoeffs(275) = 0.12479359198166d0
        tbox_icoeffs(276) = 0.12534722297734d0
        tbox_icoeffs(277) = 0.12590224054258d0
        tbox_icoeffs(278) = 0.12645865782599d0
        tbox_icoeffs(279) = 0.12701648813838d0
        tbox_icoeffs(280) = 0.12757574495563d0
        tbox_icoeffs(281) = 0.12813644192162d0
        tbox_icoeffs(282) = 0.12869859285116d0
        tbox_icoeffs(283) = 0.12926221173306d0
        tbox_icoeffs(284) = 0.12982731273326d0
        tbox_icoeffs(285) = 0.13039391019800d0
        tbox_icoeffs(286) = 0.13096201865702d0
        tbox_icoeffs(287) = 0.13153165282700d0
        tbox_icoeffs(288) = 0.13210282761483d0
        tbox_icoeffs(289) = 0.13267555812120d0
        tbox_icoeffs(290) = 0.13324985964408d0
        tbox_icoeffs(291) = 0.13382574768241d0
        tbox_icoeffs(292) = 0.13440323793981d0
        tbox_icoeffs(293) = 0.13498234632837d0
        tbox_icoeffs(294) = 0.13556308897257d0
        tbox_icoeffs(295) = 0.13614548221329d0
        tbox_icoeffs(296) = 0.13672954261184d0
        tbox_icoeffs(297) = 0.13731528695420d0
        tbox_icoeffs(298) = 0.13790273225528d0
        tbox_icoeffs(299) = 0.13849189576328d0
        tbox_icoeffs(300) = 0.13908279496421d0
        tbox_icoeffs(301) = 0.13967544758649d0
        tbox_icoeffs(302) = 0.14026987160563d0
        tbox_icoeffs(303) = 0.14086608524907d0
        tbox_icoeffs(304) = 0.14146410700115d0
        tbox_icoeffs(305) = 0.14206395560812d0
        tbox_icoeffs(306) = 0.14266565008336d0
        tbox_icoeffs(307) = 0.14326920971271d0
        tbox_icoeffs(308) = 0.14387465405991d0
        tbox_icoeffs(309) = 0.14448200297217d0
        tbox_icoeffs(310) = 0.14509127658593d0
        tbox_icoeffs(311) = 0.14570249533272d0
        tbox_icoeffs(312) = 0.14631567994521d0
        tbox_icoeffs(313) = 0.14693085146335d0
        tbox_icoeffs(314) = 0.14754803124074d0
        tbox_icoeffs(315) = 0.14816724095114d0
        tbox_icoeffs(316) = 0.14878850259513d0
        tbox_icoeffs(317) = 0.14941183850696d0
        tbox_icoeffs(318) = 0.15003727136160d0
        tbox_icoeffs(319) = 0.15066482418194d0
        tbox_icoeffs(320) = 0.15129452034622d0
        tbox_icoeffs(321) = 0.15192638359562d0
        tbox_icoeffs(322) = 0.15256043804210d0
        tbox_icoeffs(323) = 0.15319670817639d0
        tbox_icoeffs(324) = 0.15383521887628d0
        tbox_icoeffs(325) = 0.15447599541508d0
        tbox_icoeffs(326) = 0.15511906347029d0
        tbox_icoeffs(327) = 0.15576444913262d0
        tbox_icoeffs(328) = 0.15641217891512d0
        tbox_icoeffs(329) = 0.15706227976270d0
        tbox_icoeffs(330) = 0.15771477906179d0
        tbox_icoeffs(331) = 0.15836970465040d0
        tbox_icoeffs(332) = 0.15902708482837d0
        tbox_icoeffs(333) = 0.15968694836797d0
        tbox_icoeffs(334) = 0.16034932452477d0
        tbox_icoeffs(335) = 0.16101424304889d0
        tbox_icoeffs(336) = 0.16168173419649d0
        tbox_icoeffs(337) = 0.16235182874164d0
        tbox_icoeffs(338) = 0.16302455798860d0
        tbox_icoeffs(339) = 0.16369995378435d0
        tbox_icoeffs(340) = 0.16437804853160d0
        tbox_icoeffs(341) = 0.16505887520213d0
        tbox_icoeffs(342) = 0.16574246735056d0
        tbox_icoeffs(343) = 0.16642885912855d0
        tbox_icoeffs(344) = 0.16711808529938d0
        tbox_icoeffs(345) = 0.16781018125309d0
        tbox_icoeffs(346) = 0.16850518302198d0
        tbox_icoeffs(347) = 0.16920312729665d0
        tbox_icoeffs(348) = 0.16990405144257d0
        tbox_icoeffs(349) = 0.17060799351710d0
        tbox_icoeffs(350) = 0.17131499228713d0
        tbox_icoeffs(351) = 0.17202508724725d0
        tbox_icoeffs(352) = 0.17273831863853d0
        tbox_icoeffs(353) = 0.17345472746787d0
        tbox_icoeffs(354) = 0.17417435552805d0
        tbox_icoeffs(355) = 0.17489724541840d0
        tbox_icoeffs(356) = 0.17562344056613d0
        tbox_icoeffs(357) = 0.17635298524848d0
        tbox_icoeffs(358) = 0.17708592461550d0
        tbox_icoeffs(359) = 0.17782230471365d0
        tbox_icoeffs(360) = 0.17856217251026d0
        tbox_icoeffs(361) = 0.17930557591875d0
        tbox_icoeffs(362) = 0.18005256382476d0
        tbox_icoeffs(363) = 0.18080318611323d0
        tbox_icoeffs(364) = 0.18155749369635d0
        tbox_icoeffs(365) = 0.18231553854259d0
        tbox_icoeffs(366) = 0.18307737370674d0
        tbox_icoeffs(367) = 0.18384305336097d0
        tbox_icoeffs(368) = 0.18461263282712d0
        tbox_icoeffs(369) = 0.18538616861011d0
        tbox_icoeffs(370) = 0.18616371843257d0
        tbox_icoeffs(371) = 0.18694534127078d0
        tbox_icoeffs(372) = 0.18773109739197d0
        tbox_icoeffs(373) = 0.18852104839301d0
        tbox_icoeffs(374) = 0.18931525724051d0
        tbox_icoeffs(375) = 0.19011378831258d0
        tbox_icoeffs(376) = 0.19091670744210d0
        tbox_icoeffs(377) = 0.19172408196169d0
        tbox_icoeffs(378) = 0.19253598075050d0
        tbox_icoeffs(379) = 0.19335247428278d0
        tbox_icoeffs(380) = 0.19417363467841d0
        tbox_icoeffs(381) = 0.19499953575548d0
        tbox_icoeffs(382) = 0.19583025308497d0
        tbox_icoeffs(383) = 0.19666586404771d0
        tbox_icoeffs(384) = 0.19750644789362d0
        tbox_icoeffs(385) = 0.19835208580345d0
        tbox_icoeffs(386) = 0.19920286095313d0
        tbox_icoeffs(387) = 0.20005885858082d0
        tbox_icoeffs(388) = 0.20092016605678d0
        tbox_icoeffs(389) = 0.20178687295628d0
        tbox_icoeffs(390) = 0.20265907113575d0
        tbox_icoeffs(391) = 0.20353685481209d0
        tbox_icoeffs(392) = 0.20442032064566d0
        tbox_icoeffs(393) = 0.20530956782687d0
        tbox_icoeffs(394) = 0.20620469816670d0
        tbox_icoeffs(395) = 0.20710581619134d0
        tbox_icoeffs(396) = 0.20801302924114d0
        tbox_icoeffs(397) = 0.20892644757417d0
        tbox_icoeffs(398) = 0.20984618447462d0
        tbox_icoeffs(399) = 0.21077235636628d0
        tbox_icoeffs(400) = 0.21170508293150d0
        tbox_icoeffs(401) = 0.21264448723579d0
        tbox_icoeffs(402) = 0.21359069585852d0
        tbox_icoeffs(403) = 0.21454383903007d0
        tbox_icoeffs(404) = 0.21550405077571d0
        tbox_icoeffs(405) = 0.21647146906673d0
        tbox_icoeffs(406) = 0.21744623597924d0
        tbox_icoeffs(407) = 0.21842849786106d0
        tbox_icoeffs(408) = 0.21941840550722d0
        tbox_icoeffs(409) = 0.22041611434471d0
        tbox_icoeffs(410) = 0.22142178462693d0
        tbox_icoeffs(411) = 0.22243558163851d0
        tbox_icoeffs(412) = 0.22345767591123d0
        tbox_icoeffs(413) = 0.22448824345171d0
        tbox_icoeffs(414) = 0.22552746598169d0
        tbox_icoeffs(415) = 0.22657553119170d0
        tbox_icoeffs(416) = 0.22763263300908d0
        tbox_icoeffs(417) = 0.22869897188130d0
        tbox_icoeffs(418) = 0.22977475507574d0
        tbox_icoeffs(419) = 0.23086019699690d0
        tbox_icoeffs(420) = 0.23195551952250d0
        tbox_icoeffs(421) = 0.23306095235967d0
        tbox_icoeffs(422) = 0.23417673342283d0
        tbox_icoeffs(423) = 0.23530310923477d0
        tbox_icoeffs(424) = 0.23644033535269d0
        tbox_icoeffs(425) = 0.23758867682113d0
        tbox_icoeffs(426) = 0.23874840865386d0
        tbox_icoeffs(427) = 0.23991981634695d0
        tbox_icoeffs(428) = 0.24110319642551d0
        tbox_icoeffs(429) = 0.24229885702679d0
        tbox_icoeffs(430) = 0.24350711852265d0
        tbox_icoeffs(431) = 0.24472831418451d0
        tbox_icoeffs(432) = 0.24596279089449d0
        tbox_icoeffs(433) = 0.24721090990657d0
        tbox_icoeffs(434) = 0.24847304766201d0
        tbox_icoeffs(435) = 0.24974959666397d0
        tbox_icoeffs(436) = 0.25104096641631d0
        tbox_icoeffs(437) = 0.25234758443255d0
        tbox_icoeffs(438) = 0.25366989732133d0
        tbox_icoeffs(439) = 0.25500837195546d0
        tbox_icoeffs(440) = 0.25636349673238d0
        tbox_icoeffs(441) = 0.25773578293495d0
        tbox_icoeffs(442) = 0.25912576620222d0
        tbox_icoeffs(443) = 0.26053400812100d0
        tbox_icoeffs(444) = 0.26196109795067d0
        tbox_icoeffs(445) = 0.26340765449454d0
        tbox_icoeffs(446) = 0.26487432813332d0
        tbox_icoeffs(447) = 0.26636180303776d0
        tbox_icoeffs(448) = 0.26787079957991d0
        tbox_icoeffs(449) = 0.26940207696481d0
        tbox_icoeffs(450) = 0.27095643610742d0
        tbox_icoeffs(451) = 0.27253472278280d0
        tbox_icoeffs(452) = 0.27413783108146d0
        tbox_icoeffs(453) = 0.27576670720623d0
        tbox_icoeffs(454) = 0.27742235365205d0
        tbox_icoeffs(455) = 0.27910583381629d0
        tbox_icoeffs(456) = 0.28081827709420d0
        tbox_icoeffs(457) = 0.28256088452207d0
        tbox_icoeffs(458) = 0.28433493504105d0
        tbox_icoeffs(459) = 0.28614179246526d0
        tbox_icoeffs(460) = 0.28798291325225d0
        tbox_icoeffs(461) = 0.28985985518940d0
        tbox_icoeffs(462) = 0.29177428712984d0
        tbox_icoeffs(463) = 0.29372799993420d0
        tbox_icoeffs(464) = 0.29572291880292d0
        tbox_icoeffs(465) = 0.29776111721759d0
        tbox_icoeffs(466) = 0.29984483275103d0
        tbox_icoeffs(467) = 0.30197648505630d0
        tbox_icoeffs(468) = 0.30415869640678d0
        tbox_icoeffs(469) = 0.30639431523595d0
        tbox_icoeffs(470) = 0.30868644322070d0
        tbox_icoeffs(471) = 0.31103846657066d0
        tbox_icoeffs(472) = 0.31345409233622d0
        tbox_icoeffs(473) = 0.31593739073687d0
        tbox_icoeffs(474) = 0.31849284475431d0
        tbox_icoeffs(475) = 0.32112540854577d0
        tbox_icoeffs(476) = 0.32384057663740d0
        tbox_icoeffs(477) = 0.32664446638658d0
        tbox_icoeffs(478) = 0.32954391690101d0
        tbox_icoeffs(479) = 0.33254660853587d0
        tbox_icoeffs(480) = 0.33566120835051d0
        tbox_icoeffs(481) = 0.33889754862615d0
        tbox_icoeffs(482) = 0.34226684793211d0
        tbox_icoeffs(483) = 0.34578198755313d0
        tbox_icoeffs(484) = 0.34945786091584d0
        tbox_icoeffs(485) = 0.35331182031270d0
        tbox_icoeffs(486) = 0.35736425629338d0
        tbox_icoeffs(487) = 0.36163935709830d0
        tbox_icoeffs(488) = 0.36616613153999d0
        tbox_icoeffs(489) = 0.37097976970068d0
        tbox_icoeffs(490) = 0.37612364937436d0
        tbox_icoeffs(491) = 0.38165175247824d0
        tbox_icoeffs(492) = 0.38763384153846d0
        tbox_icoeffs(493) = 0.39415707842176d0
        tbox_icoeffs(494) = 0.40135193157999d0
        tbox_icoeffs(495) = 0.40934625195816d0
        tbox_icoeffs(496) = 0.41852139102809d0
        tbox_icoeffs(497) = 0.42870670111883d0
        tbox_icoeffs(498) = 0.44247426151116d0
        tbox_icoeffs(499) = 0.45345768744521d0
        tbox_icoeffs(500) = 0.49146550996677d0
        tbox_icoeffs(501) = 0.45882630283799d0
        tbox_icoeffs(502) = 0.76343194636976d0
        tbox_icoeffs(503) = 1.68678349423578d0

        ! derivative of cumulative normal distribution
        tbox_dernormal = 1.486719514734298d0 * 1E-6

        ! derivative of inverse cumulative normal distribution
        tbox_dernormalInv = 1d0 / tbox_dernormal

        ! set normal loaded to true
        tbox_normal_loaded = .true.

    end subroutine loadNormal















!##############################################################################
!##############################################################################
! MODULE polynomial
!##############################################################################
!##############################################################################


    !##############################################################################
    ! FUNCTION poly_interpol_1
    !
    ! Constructs interpolating polynomial given nodes xi and data yi.
    !##############################################################################
    function poly_interpol_1(x, xi, yi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate polynomial
        real*8, intent(in) :: x

        ! nodes of interpolation
        real*8, intent(in) :: xi(0:)

        ! data of interpolation
        real*8, intent(in) :: yi(0:)

        ! return value
        real*8 :: poly_interpol_1


        !##### OTHER VARIABLES ####################################################

        integer :: j, n, i
        real*8 :: lagrange_poly(0:size(xi, 1)-1)

        !##### ROUTINE CODE #######################################################

        ! get number of interpolation nodes
        n = assert_eq(size(xi, 1), size(yi, 1), 'poly_interpol') - 1

        ! initialize lagrange basis polynomials
        lagrange_poly(:) = 1d0

        ! span polynomials
        do j = 0, n
            do i = 0, n
                if(j /= i)then
                    lagrange_poly(j) = lagrange_poly(j) * (x-xi(i))/(xi(j)-xi(i))
                endif
            enddo
        enddo

        poly_interpol_1 = sum(yi*lagrange_poly, 1)

    end function poly_interpol_1


    !##############################################################################
    ! FUNCTION poly_interpol_m
    !
    ! Constructs interpolating polynomial given nodes xi and data yi and several
    !     points x.
    !##############################################################################
    function poly_interpol_m(x, xi, yi)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate polynomial
        real*8, intent(in) :: x(1:)

        ! nodes of interpolation
        real*8, intent(in) :: xi(0:)

        ! data of interpolation
        real*8, intent(in) :: yi(0:)

        ! return value
        real*8 :: poly_interpol_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        integer :: j, n, i
        real*8 :: lagrange_poly(size(x, 1), 0:size(xi, 1))


        !##### ROUTINE CODE #######################################################

        ! get number of interpolation nodes
        n = assert_eq(size(xi, 1), size(yi, 1), 'poly_interpol') - 1

        ! initialize lagrange basis polynomials
        lagrange_poly(:, :) = 1d0

        ! span polynomials
        do j = 0, n
            do i = 0, n
                if(j /= i)then
                    lagrange_poly(:, j) = lagrange_poly(:, j) * &
                        (x-xi(i))/(xi(j)-xi(i))
                endif
            enddo
        enddo

        do j = 1, size(x, 1)
            poly_interpol_m(j) = sum(yi*lagrange_poly(j, :), 1)
        enddo

    end function poly_interpol_m















!##############################################################################
!##############################################################################
! MODULE minimization
!
! Large parts of the procedures were taken from:
!     Press, Teukolsky, Vetterling and Flannery (1992): "Numerical Recipes in
!     FORTRAN: The Art of Scientific Computing", 2nd edition, Cambridge
!     Univeristy Press, Cambridge.
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE settol_min
    !
    ! For setting global tolerance level.
    !##############################################################################
    subroutine settol_min(tol)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! tolerance level
        real*8, intent(in) :: tol


        !##### ROUTINE CODE #######################################################

        ! check whether tolerance level is valid
        if(tol > 0d0)then
            tbox_gftol = tol
        else
            call warning('settol_min', 'tolerance level is not valid')
        endif

    end subroutine settol_min


    !##############################################################################
    ! SUBROUTINE setiter_min
    !
    ! For setting maximum number of iterations.
    !##############################################################################
    subroutine setiter_min(iter)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! tolerance level
        integer, intent(in) :: iter


        !##### ROUTINE CODE #######################################################

        ! check whether tolerance level is valid
        if(iter > 0)then
            tbox_itermax_min = iter
        else
            call warning('setiter_min', 'number of iterations is not valid')
        endif

    end subroutine setiter_min


    !##############################################################################
    ! SUBROUTINE brent
    !
    ! Minimizes a one dimensional function.
    !##############################################################################
    subroutine brent(xmin, fret, minimum, maximum, func)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! minimum value found
        real*8, intent(inout) :: xmin

        ! function value at minimum
        real*8, intent(out) :: fret

        ! left, middle and right interval points
        real*8, intent(in) :: minimum, maximum


        !##### OTHER VARIABLES ####################################################

        real*8 :: tol
        real*8, parameter :: cgold = 0.3819660d0
        real*8, parameter :: zeps = 1.0e-3*epsilon(xmin)
        real*8 :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, &
        u, v, w, x, xm, ax, bx, cx
        integer :: iter


        !##### INTERFACES #########################################################

        ! interface for the function
        interface
            function func(p)
                implicit none
                real*8, intent(in) :: p
                real*8 :: func
            end function func
        end interface


        !##### ROUTINE CODE #######################################################

        ! set tolerance level
        tol =  tbox_gftol

        ! set ax, bx and cx
        ax = minimum
        cx = maximum
        bx = xmin

        a = min(ax, cx)
        b = max(ax, cx)
        v = bx
        w = v
        x = v
        e = 0d0
        fx = func(x)
        fv = fx
        fw = fx

        do iter = 1,tbox_itermax_min
            xm = 0.5d0*(a+b)
            tol1 = tol*abs(x)+zeps
            tol2 = 2.0d0*tol1

            if(abs(x-xm) <= (tol2-0.5d0*(b-a)))then
                xmin = x
                fret = fx
                return
            endif

            if(abs(e) > tol1)then
                r = (x-w)*(fx-fv)
                q = (x-v)*(fx-fw)
                p = (x-v)*q-(x-w)*r
                q = 2.0d0*(q-r)
                if (q > 0.0d0) p = -p
                q = abs(q)
                etemp = e
                e = d
                if(abs(p) >= abs(0.5d0*q*etemp) .or. &
                        p <= q*(a-x) .or. p >= q*(b-x))then
                    e = merge(a-x, b-x, x >= xm )
                    d = CGOLD*e
                else
                    d = p/q
                    u = x+d
                    if(u-a < tol2 .or. b-u < tol2) d = sign(tol1, xm-x)
                endif

            else
                e = merge(a-x, b-x, x >= xm )
                d = CGOLD*e
            endif

            u = merge(x+d, x+sign(tol1, d), abs(d) >= tol1 )
            fu = func(u)
            if(fu <= fx)then
                if(u >= x)then
                    a = x
                else
                b = x
                endif
                call shft(v, w, x, u)
                call shft(fv, fw, fx, fu)
            else
                if(u < x)then
                    a = u
                else
                    b = u
                endif
                if(fu <= fw .or. abs(w-x)  <= 1d-100)then
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                elseif(fu <= fv .or. abs(v-x) <= 1d-100 .or. abs(v-w) <= 1d-100)then
                    v = u
                    fv = fu
                endif
            endif
        enddo

        call warning('fminsearch', 'maximum iterations exceeded')


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! SUBROUTINE shft
        !
        ! Shifts b to a, c to b and d to c.
        !##########################################################################
        subroutine shft(a, b, c, d)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            real*8, intent(out)   :: a
            real*8, intent(inout) :: b, c
            real*8, intent(in   ) :: d


            !##### ROUTINE CODE ###################################################
            a = b
            b = c
            c = d
        end subroutine shft

    end subroutine brent


    !##############################################################################
    ! SUBROUTINE powell
    !
    ! Powell is a multidimensional function minimizer.
    !##############################################################################
    subroutine powell(p, fret, minimum, maximum, func)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! starting and ending point
        real*8, intent(inout) :: p(:)

        ! value of function in minimum
        real*8, intent(out) :: fret

        ! minimum optimization interval point
        real*8, intent(in) :: minimum(:)

        ! maximum optimization interval point
        real*8, intent(in) :: maximum(:)


        !##### OTHER VARIABLES ####################################################

        real*8 :: xi(size(p, 1), size(p, 1))
        real*8 :: ftol
        real*8, parameter :: tiny = 1.0e-25
        integer :: i, ibig, n, iter
        real*8 :: del, fp, fptt, t
        real*8, dimension(size(p)) :: pt, ptt, xit
        real*8 :: xicom(size(p,1))


        !##### INTERFACES #########################################################

        ! interface for the function
        interface
            function func(p)
                implicit none
                real*8, intent(in) :: p(:)
                real*8 :: func
            end function func
        end interface


        !##### ROUTINE CODE #######################################################

        ! set tolerance level
        ftol = tbox_gftol

        ! Set number of points
        n = assert_eq(size(p), size(minimum, 1), size(maximum,1), 'fminsearch')

        ! initialize direction set
        xi = 0d0
        do i = 1, n
            xi(i, i) = 1d0
        enddo

        ! calculate function value
        fret = func(p)

        ! store old p
        pt(:) = p(:)

        ! start iteration
        iter = 0
        do
            ! step counter
            iter = iter+1

            ! save old function value
            fp = fret

            ! ibig will be direction of steepest decline
            ibig = 0
            del = 0.0d0

            ! iterate over all dimensions
            do i = 1, n

                ! copy direction i and store old function value
                xit(:) = xi(:,i)
                fptt = fret

                ! minimize along this direction
                call linmin(p, xit, n, fret, func)

                ! store i into i big if i is the direction of steepest decline
                if (fptt-fret > del) then
                    del=fptt-fret
                    ibig=i
                endif
            enddo

            ! termination criterion
            if (2d0*(fp - fret) <= ftol*(abs(fp) + abs(fret)) + tiny) return

            ! quit if maximum iterations reached
            if (iter == tbox_itermax_min)then
                call warning('fminsearch', 'maximum iterations exceeded')
            endif

            ! construct extrapolated point
            ptt(:) = 2d0*p(:) - pt(:)
            xit(:) = p(:) - pt(:)
            pt(:) = p(:)

            ! calculate function value at extrapolated point
            fptt = func(ptt)

            ! if function value greater than actual value
            ! -> no change of directions
            if (fptt >= fp) cycle

            ! calculate t
            t = 2d0*(fp - 2d0*fret + fptt) * (fp - fret - del)**2 &
                - del*(fp - fptt)**2

            ! if t > 0 -> no change of directions
            if (t >= 0d0) cycle

            ! else minimize along new direction and start new iteration
            call linmin(p, xit, n, fret, func)
            xi(:, ibig) = xi(:, n)
            xi(:,n) = xit(:)
        enddo


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! SUBROUTINE linmin
        !
        ! Minimizes multidimensional function along a given direction.
        !##########################################################################
        subroutine linmin(p, xi, n, fret, func)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! point where to start and minimum if minimization is done
            real*8, intent(inout) :: p(n)

            ! direction in which to optimize
            real*8, intent(inout) :: xi(n)

            ! number of dimensions
            integer, intent(in) :: n

            ! value of function at minimum
            real*8, intent(out)  :: fret


            !##### OTHER VARIABLES ################################################

            real*8 :: tol
            real*8  :: ax, bx, cx, xmin


            !##### INTERFACES #####################################################

                ! interface for the function
            interface
                function func(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: func
                end function func
            end interface


            !##### ROUTINE CODE ###################################################

            ! set tolerance level
            tol = tbox_gftol

            xicom(:) = xi(:)

            ! get optimization interval
            call getintervall(ax, bx, cx)

            ! minimize function using one dimensional optimizer
            fret = brent_pow(ax, bx, cx, tol, xmin, func)

            ! calculate new direction and endpoint
            xi(:) = xmin*xi(:)
            p(:) = p(:)+xi(:)

        end subroutine linmin


        !##########################################################################
        ! SUBROUTINE getinterval
        !
        ! Calculates optimization interval along a given direction.
        !##########################################################################
        subroutine getintervall(ax, bx, cx)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! left, middle and right interval point
            real*8, intent(out) :: ax, bx, cx


            !##### OTHER VARIABLES ################################################

            integer :: i
            real*8  :: w(0:1,n)


            !##### ROUTINE CODE ###################################################

            ! calculate right interval point
            cx = -1.e20
            do i = 1, n
                if(abs(xicom(i)) >= 1d-100)then
                    w(0,i) = (extr(i, 0, xicom(i))-p(i))/xicom(i)
                else
                    w(0,i) = -1.e20
                endif
                if(w(0,i) > cx)cx = w(0, i)
            enddo

            ! calculate left interval point
            ax = 1.e20
            do i=1, n
                if(abs(xicom(i)) >= 1d-100)then
                    w(1,i) = (extr(i, 1, xicom(i))-p(i))/xicom(i)
                else
                    w(1,i) = 1.e20
                endif
                if(w(1,i) < ax)ax = w(1,i)
            enddo

            ! calculate point in between [ax, cx]
            bx = 0d0

        end subroutine getintervall


        !##########################################################################
        ! FUNCTION extr
        !
        ! Calculates interval endpoint in a given dimension.
        !##########################################################################
        real*8 function extr(dim, maxxx, richt)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! dimension in which to search
            integer, intent(in) :: dim

            ! do you want the maximum or minimum endpoint
            integer, intent(in) :: maxxx

            ! what is the optimization direction
            real*8 :: richt


            !##### ROUTINE CODE ###################################################

            if(richt > 0d0 .and. maxxx == 1 .or. richt < 0d0 .and. maxxx == 0)then
                extr = maximum(dim)
            else
                extr = minimum(dim)
            endif

        end function extr


        !##########################################################################
        ! FUNCTION f1dim
        !
        ! Maps multidimensional function and point/direction combo into
        !     one-dimensional function.
        !##########################################################################
        function f1dim(x, func)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! point where to evaluate the multidimensional function
            real*8, intent(in)  :: x

            ! function value
            real*8 :: f1dim


            !##### OTHER VARIABLES ################################################

            real*8 :: xt(n)


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function func(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: func
                end function func
            end interface


            !##### ROUTINE CODE ###################################################

            ! create point where to evaluate func
            xt(:) = p(:)+x*xicom(:)

            ! evaluate func at this point
            f1dim = func(xt)

        end function f1dim


        !##########################################################################
        ! FUNCTION brent_pow
        !
        ! Minimizes a one dimensional function.
        !##########################################################################
        function brent_pow(ax, bx, cx, tol, xmin, func)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! left, middle and right interval points
            real*8, intent(in) :: ax, bx, cx

            ! level of tolerance
            real*8, intent(in) :: tol

            ! minimum value found
            real*8, intent(out) :: xmin

            ! function value at minimum
            real*8 :: brent_pow


            !##### OTHER VARIABLES ################################################

            real*8, parameter :: cgold = 0.3819660d0
            real*8, parameter :: zeps=1.0e-3*epsilon(ax)
            integer :: iter
            real*8 :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, &
                u, v, w, x, xm


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function func(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: func
                end function func
            end interface


            !##### ROUTINE CODE ###################################################

            a = min(ax, cx)
            b = max(ax, cx)
            v = bx
            w = v
            x = v
            e = 0d0
            fx = f1dim(x, func)
            fv = fx
            fw = fx

            do iter = 1,tbox_tbox_itermax_pow_b
                xm = 0.5d0*(a+b)
                tol1 = tol*abs(x)+zeps
                tol2 = 2.0d0*tol1

                if(abs(x-xm) <= (tol2-0.5d0*(b-a)))then
                    xmin = x
                    brent_pow = fx
                    return
                endif

                if(abs(e) > tol1)then
                    r = (x-w)*(fx-fv)
                    q = (x-v)*(fx-fw)
                    p = (x-v)*q-(x-w)*r
                    q = 2.0d0*(q-r)
                    if (q > 0.0d0) p = -p
                    q = abs(q)
                    etemp = e
                    e = d
                    if(abs(p) >= abs(0.5d0*q*etemp) .or. &
                        p <= q*(a-x) .or. p >= q*(b-x))then
                        e = merge(a-x, b-x, x >= xm )
                        d = CGOLD*e
                    else
                        d = p/q
                        u = x+d
                        if(u-a < tol2 .or. b-u < tol2) d = sign(tol1, xm-x)
                    endif

                else
                    e = merge(a-x, b-x, x >= xm )
                    d = CGOLD*e
                endif
                u = merge(x+d, x+sign(tol1, d), abs(d) >= tol1 )
                fu = f1dim(u, func)
                if(fu <= fx)then
                    if(u >= x)then
                        a = x
                    else
                        b = x
                    endif
                    call shft(v, w, x, u)
                    call shft(fv, fw, fx, fu)
                else
                    if(u < x)then
                        a = u
                    else
                        b = u
                    endif
                    if(fu <= fw .or. abs(w-x) <= 1d-100)then
                        v = w
                        fv = fw
                        w = u
                        fw = fu
                    elseif(fu <= fv .or. abs(v-x) <= 1d-100 .or. abs(v-w) <= 1d-100)then
                        v = u
                        fv = fu
                    endif
               endif
            enddo
            call warning('brent', 'maximum iterations exceeded')

        end function brent_pow


        !##########################################################################
        ! SUBROUTINE shft
        !
        ! Shifts b to a, c to b and d to c.
        !##########################################################################
        subroutine shft(a, b, c, d)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            real*8, intent(out)   :: a
            real*8, intent(inout) :: b, c
            real*8, intent(in   ) :: d


            !##### ROUTINE CODE ###################################################
            a = b
            b = c
            c = d
        end subroutine shft

    end subroutine powell















!##############################################################################
!##############################################################################
! MODULE simplex
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE solve_lin
    !
    ! For solving a linear program in normal form by means of the simplex
    !   algorithm.
    !##############################################################################
    subroutine solve_lin(x, c, A, b, numle, numge, numeq)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! solution of the linear program
        real*8, intent(inout) :: x(:)

        ! coefficients in the function to minimize
        real*8, intent(in) :: c(:)

        ! constraint matrix
        real*8, intent(in) :: A(:, :)

        ! target vectors of constraint
        real*8, intent(in) :: b(:)

        ! number of lower equal constraints
        integer, intent(in) :: numle

        ! number of greater equal constraints
        integer, intent(in) :: numge

        ! number of equality constraints
        integer, intent(in) :: numeq


        !##### OTHER VARIABLES ####################################################

        integer :: m, n, i1
        real*8 :: A_h(size(A, 1), size(A, 2)+numle+numge)
        real*8 :: c_h(size(c, 1)+numle+numge), b_h(size(b, 1))
        real*8 :: x_h(size(x, 1)+numle+numge)
        real*8 :: newA(size(A, 1), size(A, 2)+numle+numge)


        !##### ROUTINE CODE #######################################################

        ! check for sizes
        n = assert_eq(size(x, 1), size(c, 1), size(A, 2), 'solve_lin')
        m = assert_eq(size(A, 1), size(b, 1), 'solve_lin')

        ! check for correct inputs
        if(numle < 0)then
           call error('solve_lin', 'Number of lower equal constraints must '// &
               'not be negative')
        elseif(numge < 0)then
           call error('solve_lin', 'Number of greater equal constraints must '// &
               'not be negative')
        elseif(numeq < 0)then
           call error('solve_lin', 'Number of equality constraints must '// &
               'not be negative')
        elseif(numle+numge+numeq /= size(b,1))then
           call error('solve_lin', 'Number of equations does not match size of b')
        endif

        ! set up optimization problem
        A_h = 0d0
        A_h(1:m, 1:n) = A(:, :)
        do i1 = 1, numle
            A_h(i1, n+i1) = 1d0
        enddo
        do i1 = 1, numge
            A_h(numle+i1, n+numle+i1) = -1d0
        enddo

        ! check for negative bs
        b_h = b
        do i1 = 1, m
            if(b(i1) < 0d0)then
                A_h(i1, :) = -A_h(i1, :)
                b_h(i1) = -b_h(i1)
            endif
        enddo

        ! initialize c
        c_h = 0d0
        c_h(1:n) = c(:)

        call get_starting_value(x_h, A_h, b_h, newA)

        call solve_simplex(x_h, c_h, newA)

        x = x_h(1:n)

    contains


        !##############################################################################
        ! SUBROUTINE get_starting_value
        !
        ! Calculates a starting value for the linear program.
        !##############################################################################
        subroutine get_starting_value(x0, A, b, newA)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! starting value of the linear program
            real*8, intent(out) :: x0(:)

            ! constraint matrix
            real*8, intent(in) :: A(:, :)

            ! target vectors of constraint
            real*8, intent(in) :: b(:)

            ! new matrix for simplex
            real*8, intent(out) :: newA(:, :)


            !##### OTHER VARIABLES ####################################################

            integer :: m, n, j
            real*8 :: Astart(size(A,1), size(A,1)+size(A,2))
            real*8 :: cstart(size(A,1)+size(A,2)), xstart(size(A,1)+size(A,2))

            !##### ROUTINE CODE #######################################################

            ! get sizes
            n = size(A, 2)
            m = size(A, 1)

            ! set up help problem
            cstart(1:n) = 0d0
            cstart(n+1:n+m) = 1d0

            ! set up help matrix
            Astart(1:m, 1:n) = A
            Astart(1:m, n+1:n+m) = 0d0
            do j = 1, m
                Astart(j,n+j) = 1d0
            enddo

            ! get initial guess
            xstart(1:n) = 0d0
            xstart(n+1:n+m) = b

            ! solve linear program
            call solve_simplex(xstart, cstart, Astart, newA)

            ! set starting value
            x0 = xstart(1:n)

        end subroutine get_starting_value


        !##############################################################################
        ! SUBROUTINE solve_simplex
        !
        ! Solves a linear program in canonic form.
        !##############################################################################
        subroutine solve_simplex(x, c, A, newA)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #############################################

            ! starting value of the linear program and result
            real*8, intent(inout) :: x(:)

            ! coefficients of the program
            real*8, intent(in) :: c(:)

            ! constraint matrix
            real*8, intent(in) :: A(:, :)

            ! new tableau if starting value calculation
            real*8, intent(out), optional :: newA(:, :)


            !##### OTHER VARIABLES ####################################################

            integer :: k, m, n, j, ibas, inot, piv(2), ihelp, i1, i2, n1
            integer :: bas(size(A, 1)), nbas(size(A, 2)-size(A, 1))
            real*8 :: alpha(size(A, 1), size(A,2)-size(A, 1))
            real*8 :: gamma(size(A, 2)-size(A,1)), x0(size(A, 1))
            real*8 :: pivcheck(size(A,1)), phelp, bhelp, alpha_h(size(alpha, 2))

            !##### ROUTINE CODE #######################################################

            ! get sizes
            n = size(A, 2)
            m = size(A, 1)
            k = size(A, 2)-size(A,1)

            ! set up basis and non basis elements
            ibas = 1
            inot = 1
            do j = 1, n
                if(abs(x(j)) >= 1d-100)then
                    bas(ibas) = j
                    ibas = ibas + 1
                else
                    nbas(inot) = j
                    inot = inot + 1
                endif
            enddo

            ! set up x0 and c'x
            do ibas = 1, m
                x0(ibas) = x(bas(ibas))
            enddo

            ! set up alphas
            do inot = 1, k
                alpha(:, inot) = -A(:, nbas(inot))
            enddo

            ! set up gammas
            do inot = 1, k
                gamma(inot) = 0d0
                do ibas = 1, m
                    gamma(inot) = gamma(inot) + alpha(ibas, inot)*c(bas(ibas))
                enddo
                gamma(inot) = gamma(inot) + c(nbas(inot))
            enddo

            ! start algorithm
            do

                ! choose pivot column
                piv = 0
                do inot = 1, k
                    if(gamma(inot) < 0d0)then
                        piv(2) = inot
                        exit
                    endif
                enddo

                ! algorithm ends of no gamma < 0
                if(piv(2) == 0)exit

                ! else choose pivot row
                do ibas = 1, m
                    if(abs(alpha(ibas, piv(2))) >= 1d-100)then
                        pivcheck(ibas) = x0(ibas)/alpha(ibas, piv(2))
                    else
                        pivcheck(ibas) = -1d300
                    endif
                enddo
                phelp = -1d300
                do ibas = 1, m
                    if(alpha(ibas, piv(2)) < 0d0 .and. pivcheck(ibas) > phelp)then
                        phelp = pivcheck(ibas)
                        piv(1) = ibas
                    endif
                enddo

                ! no solution in piv(1) == 0
                if(piv(1) == 0)then
                    call error('solve_lin','Problem has no solution')
                endif

                ! Apply basis change
                Ihelp = nbas(piv(2))
                nbas(piv(2)) = bas(piv(1))
                bas(piv(1)) = Ihelp

                ! change pivot element
                alpha(piv(1), piv(2)) = 1d0/alpha(piv(1), piv(2))

                ! change pivot column
                do ibas = 1, m
                    if(ibas /= piv(1))then
                        alpha(ibas, piv(2)) = alpha(ibas, piv(2))* &
                            alpha(piv(1), piv(2))
                    endif
                enddo

                ! change pivot row
                do inot = 1, k
                    if(inot /= piv(2))then
                        alpha(piv(1), inot) = -alpha(piv(1), inot)* &
                            alpha(piv(1), piv(2))
                    endif
                enddo

                ! change other elements of alpha
                do ibas = 1, m
                    do inot = 1, k
                        if(ibas /= piv(1) .and. inot /= piv(2))then
                            alpha(ibas, inot) = alpha(ibas, inot) + &
                                alpha(ibas, piv(2))*alpha(piv(1), inot)/ &
                                alpha(piv(1), piv(2))
                        endif
                    enddo
                enddo

                ! change x0
                x0(piv(1)) = -x0(piv(1))*alpha(piv(1), piv(2))
                do ibas = 1, m
                    if(ibas /= piv(1))then
                        x0(ibas) = x0(ibas) + alpha(ibas, piv(2))*x0(piv(1))/ &
                            alpha(piv(1), piv(2))
                    endif
                enddo

                ! change gammas
                gamma(piv(2)) = gamma(piv(2))*alpha(piv(1), piv(2))
                do inot = 1, k
                    if(inot /= piv(2))then
                        gamma(inot) = gamma(inot) + alpha(piv(1), inot)* &
                            gamma(piv(2))/alpha(piv(1), piv(2))
                    endif
                enddo

            enddo

            ! get solution
            x = 0d0
            do ibas = 1, m
                x(bas(ibas)) = x0(ibas)
            enddo

            ! set up new tableau if needed
            if(present(newA))then

                ! check for existance of a solution
                do i1 = 1, n
                    if(c(i1) > 0d0)then
                        n1 = i1
                        exit
                    endif
                enddo
                if(any(x(n1:n) > 0d0))then
                    call error('solve_lin', &
                        'Linear Program does not have a solution')
                endif

                ! sort tableau in ascending order
                do i1 = m-1,1,-1
                    do i2 = 1, i1
                        if(bas(i2) > bas(i2+1))then
                            bhelp = bas(i2)
                            bas(i2) = bas(i2+1)
                            bas(i2+1) = bhelp

                            alpha_h = alpha(i2, :)
                            alpha(i2, :) = alpha(i2+1, :)
                            alpha(i2+1, :) = alpha_h
                        endif
                    enddo
                enddo

                ! get new matrix
                newA = 0d0
                do i1 = 1, k
                    if(nbas(i1) <= n1-1)then
                        newA(:, nbas(i1)) = -alpha(:, i1)
                    endif
                enddo
            endif

        end subroutine solve_simplex

    end subroutine solve_lin















!##############################################################################
!##############################################################################
! MODULE rootfinding
! Large parts of the procedures were taken from:
!     Press, Teukolsky, Vetterling and Flannery (1992): "Numerical Recipes in
!     FORTRAN: The Art of Scientific Computing", 2nd edition, Cambridge
!     Univeristy Press, Cambridge.
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE settol_root
    !
    ! For setting global tolerance level.
    !##############################################################################
    subroutine settol_root(tol)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! tolerance level
        real*8, intent(in) :: tol


        !##### ROUTINE CODE #######################################################

        ! check whether tolerance level is valid
        if(tol > 0d0)then
            tbox_gftol_root = tol
        else
            call warning('settol_root', 'tolerance level is not valid')
        endif

    end subroutine settol_root


    !##############################################################################
    ! SUBROUTINE setiter_root
    !
    ! For setting maximum number of iterations.
    !##############################################################################
    subroutine setiter_root(iter)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! number of iterations
        integer, intent(in) :: iter


        !##### ROUTINE CODE #######################################################

        ! check whether number of iterations is valid
        if(iter > 0)then
            itermax_root = iter
        else
            call warning('setiter_root', 'number of iterations is not valid')
        endif

    end subroutine setiter_root


    !##############################################################################
    ! SUBROUTINE newton_interpol
    !
    ! Find root of one-dimensional function by interpolatory newton method.
    !##############################################################################
    subroutine newton_interpol(x, funcv, check_return)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! initial guess and root of the function
        real*8, intent(inout) :: x

        ! check is true if broydn converged to local minimum or can make no
        !     further progress
        logical, intent(out), optional :: check_return

        !##### OTHER VARIABLES ####################################################

        real*8, parameter :: eps = epsilon(x)
        real*8 :: tolf, tolmin
        real*8, parameter :: tolx = eps
        real*8, parameter :: stpmx = 100d0
        real*8 :: x1, x2, f1, f2, xnew, fnew
        integer :: its


        !##### INTERFACES #########################################################

        ! interface for the function
        interface
            function funcv(p)
                implicit none
                real*8, intent(in) :: p
                real*8 :: funcv
            end function funcv
        end interface


        !##### ROUTINE CODE #######################################################

        ! set tolerance levels
        tolf = tbox_gftol_root
        tolmin = tbox_gftol_root

        ! initialize values
        x1 = x
        x2 = x + 0.001d0

        ! calculate function values at x1, x2
        f1 = funcv(x1)
        f2 = funcv(x2)

        ! check if already in zero
        if(abs(f1) < tolf)then
            x = x1
            if(present(check_return))check_return = .false.
            return
        endif

        if(abs(f2) < tolf)then
            x = x2
            if(present(check_return))check_return = .false.
            return
        endif

        ! start iteration
        do its = 1, itermax_root

            ! calculate new point xnew
            xnew = x2 - (x2-x1)/(f2-f1)*f2

            ! calculate new function value
            fnew = funcv(xnew)

            ! check wether function is small enough
            if(abs(fnew) < tolf)then
                x = xnew
                if(present(check_return))check_return = .false.
                return
            endif

            ! check whether you are in a minimum or cannot proceed further
            if(abs((f2-f1)/(x2-x1)) < tolmin .or. &
                2d0*abs(xnew-x2) < tolx*abs(xnew+x2))then
                x = x2
                if(present(check_return))check_return = .true.
                return
            endif

            ! else set new data and repeat step
            x1 = x2
            f1 = f2
            x2 = xnew
            f2 = fnew
        enddo

        ! throw warning if newton didn't converge
        call warning('fzero', 'fzero exceeds maximum iterations')

        x = xnew

    end subroutine newton_interpol


    !##############################################################################
    ! SUBROUTINE broydn
    !
    ! Find root of multidimensional function of.
    !##############################################################################
    subroutine broydn(x, funcv, check_return)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! initial guess and root of the function
        real*8, intent(inout) :: x(:)

        ! check is true if broydn converged to local minimum or can make no
        !     further progress
        logical, intent(out), optional :: check_return


        !##### OTHER VARIABLES ####################################################

        real*8, parameter :: eps = epsilon(x)
        real*8 :: tolf, tolmin
        real*8, parameter :: tolx = eps
        real*8, parameter :: stpmx = 100d0
        integer :: i, j, k, its, n
        real*8 :: f, fold, stpmax
        real*8, dimension(size(x)) :: c, d, fvcold, g, p, s, t, w, xold
        real*8, dimension(size(x),size(x)) :: qt, r
        logical :: restrt, sing, check
        real*8 :: fvec(size(x,1))


        !##### INTERFACES #########################################################

        ! interface for the function
        interface
            function funcv(p)
                implicit none
                real*8, intent(in) :: p(:)
                real*8 :: funcv(size(p, 1))
            end function funcv
        end interface


        !##### ROUTINE CODE #######################################################

        ! set tolerance levels
        tolf = tbox_gftol_root
        tolmin = tbox_gftol_root

        ! get size of x
        n = size(x)

        ! calculate function euklidean norm at starting point
        f = fmin(x, funcv)

        ! check if root has been found
        if (maxval(abs(fvec(:))) < 0.01d0*tolf) then
            if(present(check_return))check_return = .false.
                return
        endif

        stpmax = stpmx*max(sqrt(dot_product(x(:),x(:))),dble(n))
        restrt=.true.

        ! iterate broydn steps
        do its=1,itermax_root

            ! If restart then calculate jacobian of function
            if (restrt) then

                ! calculate jacobian of func at x
                call fdjac(x, fvec, r, funcv)

                ! make q-r-decomposition of jacobian
                call qrdcmp(r, c, d, sing)

                ! throw error if jacobian is singular
                if(sing)then
                    call warning('broydn', 'singular jacobian in broydn')
                    return
                endif

                ! create unity matrix
                qt(:,:) = 0d0
                    do j = 1, n
                            qt(j, j) = 1d0
                    enddo

                ! for Q^T explicitly
                do k = 1, n-1
                    if (abs(c(k)) >= 1d-100) then
                        qt(k:n,:) = qt(k:n, :)-outerprod(r(k:n, k), &
                            matmul(r(k:n, k), qt(k:n, :)))/c(k)
                    endif
                enddo
                where(lower_triangle(n,n))r(:, :) = 0d0

                ! puts diagonal elements of R matrix to r
                do j = 1, n
                            r(j, j) = d(j)
                    enddo

            ! else do Broydn update step
            else

                ! set up s as delta x
                s(:) = x(:)-xold(:)

                ! t = R*delta x
                do i = 1, n
                    t(i) = dot_product(r(i,i:n), s(i:n))
                enddo

                ! w = delta f - B*s = delta f - R*s*Q^T
                w(:) = fvec(:)-fvcold(:)-matmul(t(:), qt(:,:))

                ! if w entries are small enough, set them to zero
                where(abs(w(:)) < EPS*(abs(fvec(:))+abs(fvcold(:)))) w(:) = 0d0

                ! update for non-noisy components of w
                if(any(abs(w(:)) >= 1d-100))then

                    ! update t and s
                    t(:) = matmul(qt(:,:),w(:))
                    s(:)=s(:)/dot_product(s,s)

                    ! update R and Q^T
                    call qrupdt(r,qt,t,s)

                    ! get diagonal of matrix r
                    do j = 1, size(r,1)
                        d(j) = r(j,j)
                    enddo

                    ! if any diagonal value of r is 0, then jacobian is singular
                    if(any(abs(d(:)) <= 1d-100))then
                        call warning('broydn', 'singular jacobian in broydn')
                        return
                    endif
                endif
            endif

            ! perform the newton step by inverting jacobian
            p(:) = -matmul(qt(:,:), fvec(:))
            do i = 1, n
                g(i) = -dot_product(r(1:i,i), p(1:i))
            enddo

            ! store old x, function value and function norm
            xold(:) = x(:)
            fvcold(:) = fvec(:)
            fold = f

            ! solve linear equation with upper triangular matrix r
            call rsolv(r, d, p)

            ! searches along the new gradient direction for new x and f
            call lnsrch(xold, fold, g, p, x, f, stpmax, check, funcv)

            ! check whether root was found
            if(maxval(abs(fvec(:))) < tolf)then
                if(present(check_return))check_return = .false.
                return
            endif

                    ! if check is true
            if(check)then

                ! check if improvement can be made, if not, return
                if(restrt .or. maxval(abs(g(:))*max(abs(x(:)), &
                        1d0)/max(f, 0.5d0*n)) < tolmin)then
                    if(present(check_return))check_return = check
                    return
                endif

                ! else calculate new jacobian
                restrt=.true.

            ! if check is false
            else

                ! do broydn step
                restrt=.false.

                ! check for convergence
                if (maxval((abs(x(:)-xold(:)))/max(abs(x(:)), &
                    1.0d0)) < tolx)then
                    if(present(check_return))check_return = check
                    return
                endif
            endif
        enddo

        ! throw warning if broydn didn't converge
        call warning('fzero', 'fzero exceeds maximum iterations')


    !##### SUBROUTINES AND FUNCTIONS ##########################################

    contains


        !##########################################################################
        ! FUNCTION fdjac
        !
        ! Calculates finite difference jacobian.
        !##########################################################################
        subroutine fdjac(x, fvec, df, funcv)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! value where to calculate finite difference jacobian
            real*8, intent(inout) :: x(:)

            ! function value at x
            real*8, intent(in) :: fvec(:)

            ! resulting finite difference jacobian
            real*8, intent(out) :: df(:, :)


            !##### OTHER VARIABLES ################################################

            real*8, parameter :: eps = 1.0e-6
            integer :: j, n
            real*8, dimension(size(x)) :: xsav, xph, h


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function funcv(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: funcv(size(p, 1))
                end function funcv
            end interface


            !##### ROUTINE CODE ###################################################

            ! check equality of sizes
            n = assert_eq(size(x), size(fvec), size(df,1), size(df,2), 'fdjac')

            ! store old x
            xsav = x

            ! calculate difference
            h = eps*abs(xsav)
            where(abs(h) <= 1d-100)h = EPS

            ! calculate x + h
            xph = xsav + h
            h = xph - xsav

            ! itertate over dimensions and calculate difference
            do j = 1, n
                x(j) = xph(j)
                df(:,j) = (funcv(x)-fvec(:))/h(j)
                x(j) = xsav(j)
            enddo

        end subroutine fdjac


        !##########################################################################
        ! FUNCTION lnsrch
        !
        ! Finds point along a line, given function value and gradient, where
        !     function has decreased sufficiently (for one dimensional function).
        !##########################################################################
        subroutine lnsrch(xold, fold, g, p, x, f, stpmax, check, funcv)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! point where to start line search
            real*8, intent(in) :: xold(:)

            ! the old function value
            real*8, intent(in) :: fold

            ! gradient at this point
            real*8, intent(in) :: g(:)

            ! a line search direction
            real*8, intent(inout) :: p(:)

            ! new value along the search line
            real*8, intent(out) :: x(:)

            ! function value at new x
            real*8, intent(out) :: f

            ! maximum size of steps such that lnsrch does not search un undefined
            !     areas
            real*8, intent(in) :: stpmax

            ! is true if x is too close at xold
            logical, intent(out) :: check


            !##### OTHER VARIABLES ################################################

            real*8, parameter :: alf = 1.0e-4
            real*8, parameter :: tolx = epsilon(x)
            integer :: ndum
            real*8 :: a, alam, alam2, alamin, b, disc, f2, pabs, rhs1, rhs2, &
                slope, tmplam


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function funcv(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: funcv(size(p, 1))
                end function funcv
            end interface


            !##### ROUTINE CODE ###################################################

            ! assert sizes or arrays
            ndum = assert_eq(size(g), size(p), size(x), size(xold), 'lnsrch')
            ndum = ndum

            ! set check's default value
            check=.false.

            ! calculate norm of p
            pabs = sqrt(dot_product(p, p))

            ! restrict p to maximum stepsize
            if(pabs > stpmax)p(:) = p(:)*stpmax/pabs

            ! calculate slope
            slope = dot_product(g, p)

            ! throw error if you would go uphill
            if(slope >= 0d0)then
                call warning('lnsrch', 'roundoff problem, I cannot go uphill')
                return
            endif

            ! calculate newton stepsize
            alamin = tolx/maxval(abs(p(:))/max(abs(xold(:)),1d0))
            alam = 1d0

            ! start iteration
            do
                ! calculate calculate new x
                x(:) = xold(:)+alam*p(:)

                ! calculate new function value at x
                f = fmin(x, funcv)

                ! if new x is not away enough return with check=true
                if(alam < alamin)then
                    x(:) = xold(:)
                    check = .true.
                    return

                ! if optimal value found return with false
                elseif(f <= fold+alf*alam*slope)then
                    return

                ! else do backtracking
                else
                    if(abs(alam -1d0) <= 1d-100)then
                        tmplam = -slope/(2d0*(f-fold-slope))
                    else
                        rhs1 = f-fold-alam*slope
                        rhs2 = f2-fold-alam2*slope
                        a = (rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
                        b = (-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
                        if(abs(a) <= 1d-100)then
                            tmplam = -slope/(2d0*b)
                        else
                            disc = b*b-3d0*a*slope
                            if(disc < 0d0)then
                                tmplam = 0.5d0*alam
                            elseif(b <= 0d0)then
                                tmplam = (-b+sqrt(disc))/(3d0*a)
                            else
                                tmplam = -slope/(b+sqrt(disc))
                            endif
                        endif
                        if(tmplam > 0.5d0*alam)tmplam = 0.5d0*alam
                    endif
                endif
                alam2 = alam
                f2 = f
                alam = max(tmplam,0.1d0*alam)
            enddo

        end subroutine lnsrch


        !##########################################################################
        ! FUNCTION fmin
        !
        ! Calculates vector norm of multidimensional function.
        !##########################################################################
        function fmin(x, funcv)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! value where to evaluate function
            real*8, intent(in) :: x(:)

            ! euklidean square norm of function at x
            real*8 :: fmin


            !##### INTERFACES #####################################################

            ! interface for the function
            interface
                function funcv(p)
                    implicit none
                    real*8, intent(in) :: p(:)
                    real*8 :: funcv(size(p, 1))
                end function funcv
            end interface

            ! calculate function value
            fvec = funcv(x)

            ! calculate squared norm
            fmin = 0.5d0*dot_product(fvec, fvec)

        end function fmin


        !##########################################################################
        ! SUBROUTINE qrdcmp
        !
        ! Calculates QR decomposition of a matrix.
        !##########################################################################
        subroutine qrdcmp(a, c, d, sing)

            implicit none

            real*8, intent(inout) :: a(:, :)
            real*8, intent(out) :: c(:), d(:)
            logical, intent(out) :: sing
            integer :: k, n
            real*8 :: scale, sigma

            n = assert_eq(size(a,1), size(a,2), size(c), size(d), 'qrdcmp')
            sing = .false.
            do k = 1, n-1
                scale = maxval(abs(a(k:n, k)))
                if(abs(scale) <= 1d-100)then
                        sing = .true.
                        c(k) = 0d0
                        d(k) = 0d0
                else
                        a(k:n, k) = a(k:n, k)/scale
                        sigma = sign(sqrt(dot_product(a(k:n, k),a(k:n, k))),a(k, k))
                        a(k,k) = a(k, k)+sigma
                        c(k) = sigma*a(k, k)
                        d(k) = -scale*sigma
                        a(k:n, k+1:n) = a(k:n, k+1:n)-outerprod(a(k:n, k),&
                                matmul(a(k:n, k),a(k:n, k+1:n)))/c(k)
                endif
            enddo
            d(n) = a(n, n)
            if (abs(d(n)) <= 1d-100) sing = .true.

        end subroutine qrdcmp


        !##########################################################################
        ! SUBROUTINE qrupdt
        !
        ! Updates qr-matrices.
        !##########################################################################
        subroutine qrupdt(r,qt,u,v)

            implicit none

            real*8, intent(inout) :: r(:, :), qt(:, :)
            real*8, intent(inout) :: u(:)
            real*8, intent(in) :: v(:)
            integer :: i, k, n

            n = assert_eq((/ size(r,1), size(r,2), size(qt,1), size(qt,2), &
                size(u), size(v)/), 'qrupdt')
            k = n+1-ifirstloc(abs(u(n:1:-1)) >= 1d-100)
            if(k < 1)k=1
            do i = k-1, 1, -1
                call rotate(r,qt,i,u(i),-u(i+1))
                u(i) = pythag(u(i),u(i+1))
            enddo
            r(1,:) = r(1,:)+u(1)*v
            do i = 1,k-1
                call rotate(r,qt,i,r(i,i),-r(i+1,i))
            enddo
        end subroutine qrupdt

        !##########################################################################
        ! SUBROUTINE rsolv
        !
        ! Solves upper diagonal system.
        !##########################################################################
        subroutine rsolv(a, d, b)

            implicit none

            real*8, intent(in) :: a(:, :), d(:)
            real*8, intent(inout) :: b(:)
            integer :: i, n

            n = assert_eq(size(a,1), size(a,2), size(b), size(d), 'rsolv')
            b(n) = b(n)/d(n)
            do i = n-1, 1, -1
                    b(i) =( b(i)-dot_product(a(i, i+1:n),b(i+1:n)))/d(i)
            enddo

        end subroutine rsolv


        subroutine rotate(r, qt, i, a, b)

            implicit none

            real*8, intent(inout) :: r(:, :), qt(:, :)
            integer, intent(in) :: i
            real*8, intent(in) :: a, b
            integer :: n
            real*8 :: c, fact, s, temp(size(r,1))

            n = assert_eq(size(r,1), size(r,2), size(qt,1), size(qt,2), 'rotate')
            if(abs(a) <= 1d-100)then
                c = 0d0
                s = sign(1d0, b)
            elseif(abs(a) > abs(b))then
                fact = b/a
                c = sign(1d0/sqrt(1d0+fact**2), a)
                s = fact*c
            else
                fact = a/b
                s = sign(1d0/sqrt(1d0+fact**2), b)
                c=fact*s
            endif
            temp(i:n) = r(i, i:n)
            r(i, i:n) = c*temp(i:n)-s*r(i+1, i:n)
            r(i+1, i:n) = s*temp(i:n)+c*r(i+1, i:n)
            temp = qt(i, :)
            qt(i, :) = c*temp-s*qt(i+1, :)
            qt(i+1, :) = s*temp+c*qt(i+1, :)

        end subroutine rotate


        function pythag(a, b)

            implicit none

            real*8, intent(in) :: a, b
            real*8 :: pythag
            real*8 :: absa, absb

            absa = abs(a)
            absb = abs(b)
            if(absa > absb)then
                pythag = absa*sqrt(1d0+(absb/absa)**2)
            else
                if(abs(absb) <= 1d-100)then
                    pythag = 0d0
                else
                    pythag = absb*sqrt(1d0+(absa/absb)**2)
                endif
            endif

        end function pythag


        function ifirstloc(mask)

            logical, intent(in) :: mask(:)
            integer :: ifirstloc, loca(1)

            loca = maxloc(merge(1, 0, mask))
            ifirstloc = loca(1)
            if(.not. mask(ifirstloc))ifirstloc = size(mask)+1

        end function ifirstloc


        function lower_triangle(j, k, extra)

            integer, intent(in) :: j, k
            integer, intent(in), optional :: extra
            logical :: lower_triangle(j, k)
            integer :: n

            n = 0
            if(present(extra))n = extra

            lower_triangle = (outerdiff(arth_i(1, 1, j), arth_i(1, 1, k)) > -n)

        end function lower_triangle


        function outerdiff(a, b)

            integer, intent(in) :: a(:), b(:)
            integer :: outerdiff(size(a, 1),size(b, 1))

            outerdiff = spread(a, dim=2, ncopies=size(b, 1)) - &
                spread(b, dim=1, ncopies=size(a, 1))

        end function outerdiff


        function outerprod(a, b)

            real*8, intent(in) :: a(:), b(:)
            real*8 :: outerprod(size(a, 1),size(b, 1))

            outerprod = spread(a, dim=2, ncopies=size(b, 1)) * &
                spread(b, dim=1, ncopies=size(a, 1))

        end function outerprod


        function arth_i(first, increment, n)

            integer, intent(in) :: first, increment, n
            integer, parameter :: npar_arth = 16
            integer, parameter :: npar2_arth = 8
            integer :: arth_i(n)
            integer :: k, k2, temp

            if(n > 0)arth_i(1) = first
            if(n <= npar_arth) then
                do k = 2, n
                    arth_i(k) = arth_i(k-1) + increment
                enddo
            else
                do k = 2, npar2_arth
                    arth_i(k) = arth_i(k-1) + increment
                enddo
                temp = increment*npar2_arth
                k = npar2_arth
                do
                    if(k >= n)exit
                    k2 = k+k
                    arth_i(k+1:min(k2,n)) = temp+arth_i(1:min(k,n-k))
                    temp = temp + temp
                    k = k2
                enddo
            endif
        end function arth_i

    end subroutine broydn















!##############################################################################
!##############################################################################
! MODULE splines
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE spline_interp1
    !
    ! Subroutine for one-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp1(y, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real*8, intent(in) :: y(0:)

        ! coefficients for spline interpolation
        real*8, intent(out) :: c(1:)


        !##### OTHER VARIABLES ####################################################

        integer :: n, j
        real*8, allocatable :: r(:), d(:)


        !##### ROUTINE CODE #######################################################

        ! assert sizes for the two arrays do fit
        n = assert_eq(size(y,1)+2, size(c,1), 'spline_interp')

        ! deallocate help arrays
        if(allocated(r))deallocate(r)
        if(allocated(d))deallocate(d)

        ! allocate help arrays
        allocate(r(n))
        allocate(d(n))

        ! calculate real n
        n = n-3

        ! calculate numerical derivatives at end points
        r(1) = (2d0*y(0)-5d0*y(1)+4d0*y(2)-y(3))/6d0
        r(n+3) = (2d0*y(n)-5d0*y(n-1)+4d0*y(n-2)-y(n-3))/6d0

        ! set rest of right side of equation system
        r(2:n+2) = y(0:n)

        ! solve the spline interpolation equation system
        c(2) = (y(0)-r(1))/6d0
        c(n+2) = (y(n)-r(n+3))/6d0

        d(3) = 4d0
        r(3) = y(1)-c(2)
        r(n+1) = y(n-1)-c(n+2)

        do j = 4, n+1
            d(j) = 4d0-1d0/d(j-1)
            r(j) = r(j)-r(j-1)/d(j-1)
        enddo

        c(n+1) = r(n+1)/d(n+1)

        do j = n, 3, -1
            c(j) = (r(j)-c(j+1))/d(j)
        enddo

        c(1) = r(1)+2d0*c(2)-c(3)
        c(n+3) = r(n+3)+2d0*c(n+2)-c(n+1)

    end subroutine spline_interp1


    !##############################################################################
    ! SUBROUTINE spline_interp2
    !
    ! Subroutine for two-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp2(y, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real*8, intent(in) :: y(0:, 0:)

        ! coefficients for spline interpolation
        real*8, intent(out) :: c(1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer :: n(2), j
        real*8, allocatable :: tempc(:, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 2
            n(j) = assert_eq(size(y,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, 0:n(2)))

        ! calculate temporary coefficients
        do j = 0, n(2)
            call spline_interp1(y(:, j), tempc(:, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            call spline_interp1(tempc(j, :), c(j, :))
        enddo

    end subroutine spline_interp2


    !##############################################################################
    ! SUBROUTINE spline_interp3
    !
    ! Subroutine for three-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp3(y, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real*8, intent(in) :: y(0:, 0:, 0:)

        ! coefficients for spline interpolation
        real*8, intent(out) :: c(1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer :: n(3), j, j2
        real*8, allocatable :: tempc(:, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 3
            n(j) = assert_eq(size(y,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, 0:n(3)))

        ! calculate temporary coefficients
        do j = 0, n(3)
            call spline_interp2(y(:, :, j), tempc(:, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                call spline_interp1(tempc(j, j2, :), c(j, j2, :))
            enddo
        enddo

    end subroutine spline_interp3


    !##############################################################################
    ! SUBROUTINE spline_interp4
    !
    ! Subroutine for four-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp4(y, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real*8, intent(in) :: y(0:, 0:, 0:, 0:)

        ! coefficients for spline interpolation
        real*8, intent(out) :: c(1:, 1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer :: n(4), j, j2, j3
        real*8, allocatable :: tempc(:, :, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 4
            n(j) = assert_eq(size(y,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, n(3)+3, 0:n(4)))

        ! calculate temporary coefficients
        do j = 0, n(4)
            call spline_interp3(y(:, :, :, j), tempc(:, :, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                do j3 = 1, n(3)+3
                    call spline_interp1(tempc(j, j2, j3, :), c(j, j2, j3, :))
                enddo
            enddo
        enddo

    end subroutine spline_interp4


    !##############################################################################
    ! SUBROUTINE spline_interp5
    !
    ! Subroutine for five-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp5(y, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:)

        ! coefficients for spline interpolation
        real*8, intent(out) :: c(1:, 1:, 1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer :: n(5), j, j2, j3, j4
        real*8, allocatable :: tempc(:, :, :, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 5
            n(j) = assert_eq(size(y,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, n(3)+3, n(4)+3, 0:n(5)))

        ! calculate temporary coefficients
        do j = 0, n(5)
            call spline_interp4(y(:, :, :, :, j), tempc(:, :, :, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                do j3 = 1, n(3)+3
                    do j4 = 1, n(4)+3
                        call spline_interp1(tempc(j, j2, j3, j4, :), &
                            c(j, j2, j3, j4, :))
                    enddo
                enddo
            enddo
        enddo

    end subroutine spline_interp5


    !##############################################################################
    ! SUBROUTINE spline_interp6
    !
    ! Subroutine for six-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp6(y, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:, 0:)

        ! coefficients for spline interpolation
        real*8, intent(out) :: c(1:, 1:, 1:, 1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer :: n(6), j, j2, j3, j4, j5
        real*8, allocatable :: tempc(:, :, :, :, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 6
            n(j) = assert_eq(size(y,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, n(3)+3, n(4)+3, n(5)+3, 0:n(6)))

        ! calculate temporary coefficients
        do j = 0, n(6)
            call spline_interp5(y(:, :, :, :, :, j), tempc(:, :, :, :, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                do j3 = 1, n(3)+3
                    do j4 = 1, n(4)+3
                        do j5 = 1, n(5)+3
                            call spline_interp1(tempc(j, j2, j3, j4, j5, :), &
                                c(j, j2, j3, j4, j5, :))
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine spline_interp6


    !##############################################################################
    ! SUBROUTINE spline_interp7
    !
    ! Subroutine for seven-dimensional spline interpolation.
    !##############################################################################
    subroutine spline_interp7(y, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! interpolation data
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:, 0:, 0:)

        ! coefficients for spline interpolation
        real*8, intent(out) :: c(1:, 1:, 1:, 1:, 1:, 1:, 1:)


        !##### OTHER VARIABLES ####################################################

        integer :: n(7), j, j2, j3, j4, j5, j6
        real*8, allocatable :: tempc(:, :, :, :, :, :, :)


        !##### ROUTINE CODE #######################################################

        ! calculate array sizes
        do j = 1, 7
            n(j) = assert_eq(size(y,j)+2, size(c,j), 'spline_interp')
        enddo

        ! calculate real n
        n = n-3

        ! deallocate tempc
        if(allocated(tempc))deallocate(tempc)

        ! allocate tempc
        allocate(tempc(n(1)+3, n(2)+3, n(3)+3, n(4)+3, n(5)+3, n(6)+3, 0:n(7)))

        ! calculate temporary coefficients
        do j = 0, n(7)
            call spline_interp6(y(:, :, :, :, :, :, j), tempc(:, :, :, :, :, :, j))
        enddo

        ! calculate actual coefficients
        do j = 1, n(1)+3
            do j2 = 1, n(2)+3
                do j3 = 1, n(3)+3
                    do j4 = 1, n(4)+3
                        do j5 = 1, n(5)+3
                            do j6 = 1, n(6)+3
                                call spline_interp1(tempc(j, j2, j3, &
                                    j4, j5, j6, :),c(j, j2, j3, j4, j5, j6, :))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine spline_interp7


    !##############################################################################
    ! FUNCTION spline1
    !
    ! Function for evaluation of one-dimensional spline.
    !##############################################################################
    function spline1(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:)

        ! value of spline function
        real*8 :: spline1


        !##### OTHER VARIABLES ####################################################

        integer :: n1, j1, p1, q1
        real*8 :: phi1, xtemp1


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n1 = size(c, 1)

        ! calculate left and right summation end point
        p1 = max(floor(x)+1, 1)
        q1 = min(p1+3, n1)

        spline1 = 0d0

        do j1 = p1, q1

            ! calculate value where to evaluate basis function
            xtemp1 = abs(x-j1+2)

            ! calculate basis function
            if(xtemp1 <= 1d0)then
                phi1 = 4d0+xtemp1**2*(3d0*xtemp1-6d0)
            elseif(xtemp1 <= 2d0)then
                phi1 = (2d0-xtemp1)**3
            else
                phi1 = 0d0
            endif

            ! calculate spline value
            spline1 = spline1+c(j1)*phi1
        enddo

    end function spline1


    !##############################################################################
    ! FUNCTION spline1_grid
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline1_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:)

        ! left interval endpoint
        real*8, intent(in) :: left

        ! right interval endpoint
        real*8, intent(in) :: right

        ! growth rate of grid
        real*8, intent(in), optional :: growth

        ! value of spline function
        real*8 :: spline1_grid


        !##### OTHER VARIABLES ####################################################

        integer :: n
        real*8 :: xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        n = size(c, 1)-3

        ! invert grid
        if(present(growth))then
            xtemp = grid_Inv_Grow(x, left, right, growth, n)
        else
            xtemp = grid_Inv_Equi(x, left, right, n)
        endif

        ! calculate spline value
        spline1_grid = spline1(xtemp, c)

    end function spline1_grid


    !##############################################################################
    ! FUNCTION spline1_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline1_complete(x, y, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x

        ! data for spline interpolation
        real*8, intent(in) :: y(0:)

        ! left interval endpoint
        real*8, intent(in) :: left

        ! right interval endpoint
        real*8, intent(in) :: right

        ! growth rate of grid
        real*8, intent(in), optional :: growth

        ! value of spline function
        real*8 :: spline1_complete


        !##### OTHER VARIABLES ####################################################

        real*8 :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline1_complete_m((/x/), y, left, right, growth)
        else
            spline_temp = spline1_complete_m((/x/), y, left, right)
        endif

        ! paste data
        spline1_complete = spline_temp(1)

    end function spline1_complete


    !##############################################################################
    ! FUNCTION spline1_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline1_complete_m(x, y, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(1:)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:)

        ! left interval endpoint
        real*8, intent(in) :: left

        ! right interval endpoint
        real*8, intent(in) :: right

        ! growth rate of grid
        real*8, intent(in), optional :: growth

        ! value of spline function
        real*8 :: spline1_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: c(1:size(y, 1)+2)
        real*8 :: xtemp(1:size(x, 1))
        integer :: n, m, j


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        n = size(y, 1)-1

        ! calculate number of evaluation points
        m = size(x, 1)

        ! invert grid for every evaluation point
        if(present(growth))then
            xtemp(:) = grid_Inv_Grow(x(:), left, right, growth, n)
        else
            xtemp(:) = grid_Inv_Equi(x(:), left, right, n)
        endif

        ! interpolate data
        call spline_interp1(y, c)

        ! calculate spline values at point
        do j = 1, m
            spline1_complete_m(j) = spline1(xtemp(j), c)
        enddo

    end function spline1_complete_m


    !##############################################################################
    ! FUNCTION spline2
    !
    ! Function for evaluation of two-dimensional spline.
    !##############################################################################
    function spline2(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(2)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:)

        ! value of spline function
        real*8 :: spline2


        !##### OTHER VARIABLES ####################################################

        integer :: n(2), p(2), q(2)
        integer :: j1, j2
        real*8 :: phi1, xtemp1, phi2, xtemp2
        real*8 :: s2


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n(1) = size(c, 1)
        n(2) = size(c, 2)

        ! calculate left and right summation end point
        p = max(floor(x)+1, 1)
        q = min(p+3, n)

        spline2 = 0d0

        do j1 = p(1), q(1)

            ! calculate value where to evaluate basis function
            xtemp1 = abs(x(1)-j1+2)

            ! calculate basis function
            if(xtemp1 <= 1d0)then
                phi1 = 4d0+xtemp1**2*(3d0*xtemp1-6d0)
            elseif(xtemp1 <= 2d0)then
                phi1 = (2d0-xtemp1)**3
            else
                phi1 = 0d0
            endif


            !#### calculate spline for second dimension ###########################

            s2 = 0d0

            do j2 = p(2), q(2)

                ! calculate value where to evaluate basis function
                xtemp2 = abs(x(2)-j2+2)

                ! calculate basis function
                if(xtemp2 <= 1d0)then
                    phi2 = 4d0+xtemp2**2*(3d0*xtemp2-6d0)
                elseif(xtemp2 <= 2d0)then
                    phi2 = (2d0-xtemp2)**3
                else
                    phi2 = 0d0
                endif

                ! calculate spline value
                s2 = s2+c(j1, j2)*phi2
            enddo

            ! calculate spline value
            spline2 = spline2+s2*phi1
        enddo

    end function spline2


    !##############################################################################
    ! FUNCTION spline2_grid
    !
    ! Function for evaluation of two-dimensional spline, includ inverting grid.
    !##############################################################################
    function spline2_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(2)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:)

        ! left interval endpoint
        real*8, intent(in) :: left(2)

        ! right interval endpoint
        real*8, intent(in) :: right(2)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(2)

        ! value of spline function
        real*8 :: spline2_grid


        !##### OTHER VARIABLES ####################################################

        integer :: n, j
        real*8 :: xtemp(2)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline2_grid = spline2(xtemp, c)

    end function spline2_grid


    !##############################################################################
    ! FUNCTION spline2_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline2_complete(x, y, left, right, growth)


        integer, parameter :: dim = 2

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(dim)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline2_complete


        !##### OTHER VARIABLES ####################################################

        real*8 :: xtemp(1, dim)
        real*8 :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline2_complete_m(xtemp, y, left, right, growth)
        else
            spline_temp = spline2_complete_m(xtemp, y, left, right)
        endif

        ! paste data
        spline2_complete = spline_temp(1)

    end function spline2_complete


    !##############################################################################
    ! FUNCTION spline2_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline2_complete_m(x, y, left, right, growth)

        integer, parameter :: dim = 2

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline2_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: c(1:size(y, 1)+2, 1:size(y, 2)+2)
        real*8 :: xtemp(1:size(x, 1), dim)
        integer :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp2(y, c)

        ! calculate spline values at point
        do j = 1, m
            spline2_complete_m(j) = spline2(xtemp(j, :), c)
        enddo

    end function spline2_complete_m


    !##############################################################################
    ! FUNCTION spline3
    !
    ! Function for evaluation of three-dimensional spline.
    !##############################################################################
    function spline3(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(3)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:)

        ! value of spline function
        real*8 :: spline3


        !##### OTHER VARIABLES ####################################################

        integer :: n(3), p(3), q(3)
        integer :: j1, j2, j3
        real*8 :: phi1, xtemp1, phi2, xtemp2, phi3, xtemp3
        real*8 :: s2, s3


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n(1) = size(c, 1)
        n(2) = size(c, 2)
        n(3) = size(c, 3)

        ! calculate left and right summation end point
        p = max(floor(x)+1, 1)
        q = min(p+3, n)

        spline3 = 0d0

        do j1 = p(1), q(1)

            ! calculate value where to evaluate basis function
            xtemp1 = abs(x(1)-j1+2)

            ! calculate basis function
            if(xtemp1 <= 1d0)then
                phi1 = 4d0+xtemp1**2*(3d0*xtemp1-6d0)
            elseif(xtemp1 <= 2d0)then
                phi1 = (2d0-xtemp1)**3
            else
                phi1 = 0d0
            endif


            !#### calculate spline for second dimension ###########################

            s2 = 0d0

            do j2 = p(2), q(2)

                ! calculate value where to evaluate basis function
                xtemp2 = abs(x(2)-j2+2)

                ! calculate basis function
                if(xtemp2 <= 1d0)then
                    phi2 = 4d0+xtemp2**2*(3d0*xtemp2-6d0)
                elseif(xtemp2 <= 2d0)then
                    phi2 = (2d0-xtemp2)**3
                else
                    phi2 = 0d0
                endif


                !#### calculate spline for second dimension #######################

                s3 = 0d0

                do j3 = p(3), q(3)

                    ! calculate value where to evaluate basis function
                    xtemp3 = abs(x(3)-j3+2)

                    ! calculate basis function
                    if(xtemp3 <= 1d0)then
                        phi3 = 4d0+xtemp3**2*(3d0*xtemp3-6d0)
                    elseif(xtemp3 <= 2d0)then
                        phi3 = (2d0-xtemp3)**3
                    else
                        phi3 = 0d0
                    endif

                    ! calculate spline value
                    s3 = s3+c(j1, j2, j3)*phi3
                enddo

                ! calculate spline value
                s2 = s2+s3*phi2
            enddo

            ! calculate spline value
            spline3 = spline3+s2*phi1
        enddo

    end function spline3


    !##############################################################################
    ! FUNCTION spline3_grid
    !
    ! Function for evaluation of three-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline3_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(3)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:)

        ! left interval endpoint
        real*8, intent(in) :: left(3)

        ! right interval endpoint
        real*8, intent(in) :: right(3)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(3)

        ! value of spline function
        real*8 :: spline3_grid


        !##### OTHER VARIABLES ####################################################

        integer :: n, j
        real*8 :: xtemp(3)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline3_grid = spline3(xtemp, c)

    end function spline3_grid


    !##############################################################################
    ! FUNCTION spline3_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline3_complete(x, y, left, right, growth)


        integer, parameter :: dim = 3

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(dim)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline3_complete


        !##### OTHER VARIABLES ####################################################

        real*8 :: xtemp(1, dim)
        real*8 :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline3_complete_m(xtemp, y, left, right, growth)
        else
            spline_temp = spline3_complete_m(xtemp, y, left, right)
        endif

        ! paste data
        spline3_complete = spline_temp(1)

    end function spline3_complete


    !##############################################################################
    ! FUNCTION spline3_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline3_complete_m(x, y, left, right, growth)

        integer, parameter :: dim = 3

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline3_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: c(1:size(y, 1)+2, 1:size(y, 2)+2, 1:size(y, 3)+2)
        real*8 :: xtemp(1:size(x, 1), dim)
        integer :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp3(y, c)

        ! calculate spline values at point
        do j = 1, m
            spline3_complete_m(j) = spline3(xtemp(j, :), c)
        enddo

    end function spline3_complete_m


    !##############################################################################
    ! FUNCTION spline4
    !
    ! Function for evaluation of four-dimensional spline.
    !##############################################################################
    function spline4(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(4)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:, 1:)

        ! value of spline function
        real*8 :: spline4


        !##### OTHER VARIABLES ####################################################

        integer :: n, j, p, q
        real*8 :: phi, xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n = size(c, 1)

        ! calculate left and right summation end point
        p = max(floor(x(1))+1, 1)
        q = min(p+3, n)

        spline4 = 0d0

        do j = p, q

            ! calculate value where to evaluate basis function
            xtemp = abs(x(1)-j+2)

            ! calculate basis function
            if(xtemp <= 1d0)then
                phi = 4d0+xtemp**2*(3d0*xtemp-6d0)
            elseif(xtemp <= 2d0)then
                phi = (2d0-xtemp)**3
            else
                phi = 0d0
            endif

            ! calculate spline value
            spline4 = spline4+spline3(x(2:4), c(j, :, :, :))*phi
        enddo

    end function spline4


    !##############################################################################
    ! FUNCTION spline4_grid
    !
    ! Function for evaluation of four-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline4_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(4)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:, 1:)

        ! left interval endpoint
        real*8, intent(in) :: left(4)

        ! right interval endpoint
        real*8, intent(in) :: right(4)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(4)

        ! value of spline function
        real*8 :: spline4_grid


        !##### OTHER VARIABLES ####################################################

        integer :: n, j
        real*8 :: xtemp(4)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline4_grid = spline4(xtemp, c)

    end function spline4_grid


    !##############################################################################
    ! FUNCTION spline4_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline4_complete(x, y, left, right, growth)


        integer, parameter :: dim = 4

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(dim)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline4_complete


        !##### OTHER VARIABLES ####################################################

        real*8 :: xtemp(1, dim)
        real*8 :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline4_complete_m(xtemp, y, left, right, growth)
        else
            spline_temp = spline4_complete_m(xtemp, y, left, right)
        endif

        ! paste data
        spline4_complete = spline_temp(1)

    end function spline4_complete


    !##############################################################################
    ! FUNCTION spline4_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline4_complete_m(x, y, left, right, growth)

        integer, parameter :: dim = 4

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline4_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: c(1:size(y, 1)+2, 1:size(y, 2)+2, 1:size(y, 3)+2, &
            1:size(y, 4)+2)
        real*8 :: xtemp(1:size(x, 1), dim)
        integer :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp4(y, c)

        ! calculate spline values at point
        do j = 1, m
            spline4_complete_m(j) = spline4(xtemp(j, :), c)
        enddo

    end function spline4_complete_m


    !##############################################################################
    ! FUNCTION spline5
    !
    ! Function for evaluation of five-dimensional spline.
    !##############################################################################
    function spline5(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(5)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:, 1:, 1:)

        ! value of spline function
        real*8 :: spline5


        !##### OTHER VARIABLES ####################################################

        integer :: n, j, p, q
        real*8 :: phi, xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n = size(c, 1)

        ! calculate left and right summation end point
        p = max(floor(x(1))+1, 1)
        q = min(p+3, n)

        spline5 = 0d0

        do j = p, q

            ! calculate value where to evaluate basis function
            xtemp = abs(x(1)-j+2)

            ! calculate basis function
            if(xtemp <= 1d0)then
                phi = 4d0+xtemp**2*(3d0*xtemp-6d0)
            elseif(xtemp <= 2d0)then
                phi = (2d0-xtemp)**3
            else
                phi = 0d0
            endif

            ! calculate spline value
            spline5 = spline5+spline4(x(2:5), c(j, :, :, :, :))*phi
        enddo

    end function spline5


    !##############################################################################
    ! FUNCTION spline5_grid
    !
    ! Function for evaluation of five-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline5_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(5)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:, 1:, 1:)

        ! left interval endpoint
        real*8, intent(in) :: left(5)

        ! right interval endpoint
        real*8, intent(in) :: right(5)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(5)

        ! value of spline function
        real*8 :: spline5_grid


        !##### OTHER VARIABLES ####################################################

        integer :: n, j
        real*8 :: xtemp(5)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline5_grid = spline5(xtemp, c)

    end function spline5_grid


    !##############################################################################
    ! FUNCTION spline5_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline5_complete(x, y, left, right, growth)


        integer, parameter :: dim = 5

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(dim)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline5_complete


        !##### OTHER VARIABLES ####################################################

        real*8 :: xtemp(1, dim)
        real*8 :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline5_complete_m(xtemp, y, left, right, growth)
        else
            spline_temp = spline5_complete_m(xtemp, y, left, right)
        endif

        ! paste data
        spline5_complete = spline_temp(1)

    end function spline5_complete


    !##############################################################################
    ! FUNCTION spline5_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline5_complete_m(x, y, left, right, growth)

        integer, parameter :: dim = 5

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline5_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: c(1:size(y, 1)+2, 1:size(y, 2)+2, 1:size(y, 3)+2, &
            1:size(y, 4)+2, 1:size(y, 5)+2)
        real*8 :: xtemp(1:size(x, 1), dim)
        integer :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp5(y, c)

        ! calculate spline values at point
        do j = 1, m
            spline5_complete_m(j) = spline5(xtemp(j, :), c)
        enddo

    end function spline5_complete_m


    !##############################################################################
    ! FUNCTION spline6
    !
    ! Function for evaluation of six-dimensional spline.
    !##############################################################################
    function spline6(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(6)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:, 1:, 1:, 1:)

        ! value of spline function
        real*8 :: spline6


        !##### OTHER VARIABLES ####################################################

        integer :: n, j, p, q
        real*8 :: phi, xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n = size(c, 1)

        ! calculate left and right summation end point
        p = max(floor(x(1))+1, 1)
        q = min(p+3, n)

        spline6 = 0d0

        do j = p, q

            ! calculate value where to evaluate basis function
            xtemp = abs(x(1)-j+2)

            ! calculate basis function
            if(xtemp <= 1d0)then
                phi = 4d0+xtemp**2*(3d0*xtemp-6d0)
            elseif(xtemp <= 2d0)then
                phi = (2d0-xtemp)**3
            else
                phi = 0d0
            endif

            ! calculate spline value
            spline6 = spline6+spline5(x(2:6), c(j, :, :, :, :, :))*phi
        enddo

    end function spline6


    !##############################################################################
    ! FUNCTION spline6_grid
    !
    ! Function for evaluation of six-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline6_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(6)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:, 1:, 1:, 1:)

        ! left interval endpoint
        real*8, intent(in) :: left(6)

        ! right interval endpoint
        real*8, intent(in) :: right(6)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(6)

        ! value of spline function
        real*8 :: spline6_grid


        !##### OTHER VARIABLES ####################################################

        integer :: n, j
        real*8 :: xtemp(6)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline6_grid = spline6(xtemp, c)

    end function spline6_grid


    !##############################################################################
    ! FUNCTION spline6_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline6_complete(x, y, left, right, growth)


        integer, parameter :: dim = 6

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(dim)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline6_complete


        !##### OTHER VARIABLES ####################################################

        real*8 :: xtemp(1, dim)
        real*8 :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline6_complete_m(xtemp, y, left, right, growth)
        else
            spline_temp = spline6_complete_m(xtemp, y, left, right)
        endif

        ! paste data
        spline6_complete = spline_temp(1)

    end function spline6_complete


    !##############################################################################
    ! FUNCTION spline6_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline6_complete_m(x, y, left, right, growth)

        integer, parameter :: dim = 6

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline6_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: c(1:size(y, 1)+2, 1:size(y, 2)+2, 1:size(y, 3)+2, &
            1:size(y, 4)+2, 1:size(y, 5)+2, 1:size(y, 6)+2)
        real*8 :: xtemp(1:size(x, 1), dim)
        integer :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp6(y, c)

        ! calculate spline values at point
        do j = 1, m
            spline6_complete_m(j) = spline6(xtemp(j, :), c)
        enddo

    end function spline6_complete_m


    !##############################################################################
    ! FUNCTION spline7
    !
    ! Function for evaluation of seven-dimensional spline.
    !##############################################################################
    function spline7(x, c)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(7)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:, 1:, 1:, 1:, 1:)

        ! value of spline function
        real*8 :: spline7


        !##### OTHER VARIABLES ####################################################

        integer :: n, j, p, q
        real*8 :: phi, xtemp


        !##### ROUTINE CODE #######################################################

        ! calculate number of points used
        n = size(c, 1)

        ! calculate left and right summation end point
        p = max(floor(x(1))+1, 1)
        q = min(p+3, n)

        spline7 = 0d0

        do j = p, q

            ! calculate value where to evaluate basis function
            xtemp = abs(x(1)-j+2)

            ! calculate basis function
            if(xtemp <= 1d0)then
                phi = 4d0+xtemp**2*(3d0*xtemp-6d0)
            elseif(xtemp <= 2d0)then
                phi = (2d0-xtemp)**3
            else
                phi = 0d0
            endif

            ! calculate spline value
            spline7 = spline7+spline6(x(2:7), c(j, :, :, :, :, :, :))*phi
        enddo

    end function spline7


    !##############################################################################
    ! FUNCTION spline7_grid
    !
    ! Function for evaluation of seven-dimensional spline, includes inverting grid.
    !##############################################################################
    function spline7_grid(x, c, left, right, growth)


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(7)

        ! coefficients for spline interpolation
        real*8, intent(in) :: c(1:, 1:, 1:, 1:, 1:, 1:, 1:)

        ! left interval endpoint
        real*8, intent(in) :: left(7)

        ! right interval endpoint
        real*8, intent(in) :: right(7)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(7)

        ! value of spline function
        real*8 :: spline7_grid


        !##### OTHER VARIABLES ####################################################

        integer :: n, j
        real*8 :: xtemp(7)


        !##### ROUTINE CODE #######################################################

        ! calculate number of grid-points
        do j = 1, size(x, 1)
            n = size(c, j)-3

            ! invert grid
            if(present(growth))then
                xtemp(j) = grid_Inv_Grow(x(j), left(j), right(j), growth(j), n)
            else
                xtemp(j) = grid_Inv_Equi(x(j), left(j), right(j), n)
            endif
        enddo

        ! calculate spline value
        spline7_grid = spline7(xtemp, c)

    end function spline7_grid


    !##############################################################################
    ! FUNCTION spline7_complete
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for a single point.
    !##############################################################################
    function spline7_complete(x, y, left, right, growth)


        integer, parameter :: dim = 7

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(dim)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline7_complete


        !##### OTHER VARIABLES ####################################################

        real*8 :: xtemp(1, dim)
        real*8 :: spline_temp(1)


        !##### ROUTINE CODE #######################################################

        ! set xtemp
        xtemp(1, :) = x

        ! invert grid for every evaluation point
        if(present(growth))then
            spline_temp = spline7_complete_m(xtemp, y, left, right, growth)
        else
            spline_temp = spline7_complete_m(xtemp, y, left, right)
        endif

        ! paste data
        spline7_complete = spline_temp(1)

    end function spline7_complete


    !##############################################################################
    ! FUNCTION spline7_complete_m
    !
    ! Function for evaluation of one-dimensional spline, includes inverting grid
    !     and interpolation method for many points.
    !##############################################################################
    function spline7_complete_m(x, y, left, right, growth)

        integer, parameter :: dim = 7

        !##### INPUT/OUTPUT VARIABLES #############################################

        ! value where to evaluate spline
        real*8, intent(in) :: x(1:, 1:)

        ! data for spline interpolation
        real*8, intent(in) :: y(0:, 0:, 0:, 0:, 0:, 0:, 0:)

        ! left interval endpoint
        real*8, intent(in) :: left(dim)

        ! right interval endpoint
        real*8, intent(in) :: right(dim)

        ! growth rate of grid
        real*8, intent(in), optional :: growth(dim)

        ! value of spline function
        real*8 :: spline7_complete_m(size(x, 1))


        !##### OTHER VARIABLES ####################################################

        real*8 :: c(1:size(y, 1)+2, 1:size(y, 2)+2, 1:size(y, 3)+2, &
            1:size(y, 4)+2, 1:size(y, 5)+2, 1:size(y, 6)+2, 1:size(y, 7)+2)
        real*8 :: xtemp(1:size(x, 1), dim)
        integer :: n, m, j, k


        !##### ROUTINE CODE #######################################################

        ! calculate number of evaluation points
        m = size(x, 1)

        ! check whether x has the right dimension
        n = assert_eq(size(x, 2), dim, 'spline')

        ! invert grid for every evaluation point
        if(present(growth))then
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Grow(x(:, k), &
                        left(k), right(k), growth(k), n)
            enddo
        else
            do k = 1, dim
                ! calculate number of grid-points
                n = size(y, k)-1

                xtemp(:, k) = grid_Inv_Equi(x(:, k), left(k), right(k), n)
            enddo
        endif

        ! interpolate data
        call spline_interp7(y, c)

        ! calculate spline values at point
        do j = 1, m
            spline7_complete_m(j) = spline7(xtemp(j, :), c)
        enddo

    end function spline7_complete_m















!##############################################################################
!##############################################################################
! MODULE gaussian_int
! Large parts of the procedures were taken from:
!     Press, Teukolsky, Vetterling and Flannery (1992): "Numerical Recipes in
!     FORTRAN: The Art of Scientific Computing", 2nd edition, Cambridge
!     Univeristy Press, Cambridge.
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE legendre
    !
    ! Calculates Gauss-Legendre abscissas and weights on [x1, x2].
    !##############################################################################
    subroutine legendre(x1, x2, x, w)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! left interval point on which integral should be calculated
        real*8, intent(in) :: x1

        ! left interval point on which integral should be calculated
        real*8, intent(in) :: x2

        ! abscissas of gaussian integration formula
        real*8, intent(out) :: x(:)

        ! weights of gaussian integration formula
        real*8, intent(out) :: w(:)


        !##### OTHER VARIABLES ####################################################

        real*8, parameter :: eps = 3.0e-14
        integer, parameter :: maxits = 10
        real*8, parameter :: pi = 3.14159265358979d0
        integer :: its, j, m, n
        real*8 :: xl, xm
        real*8, dimension((size(x)+1)/2) :: p1, p2, p3, pp, z, z1
        logical, dimension((size(x)+1)/2) :: unfinished


        !##### ROUTINE CODE #######################################################

        ! assert size equality
        n = assert_eq(size(x), size(w), 'legendre')

        ! calculate only up to (n+1)/2 due to symmetry
        m = (n+1)/2

        ! calculate interval midpoint
        xm = 0.5d0*(x2+x1)

        ! calculate half of interval length
        xl = 0.5d0*(x2-x1)

        ! set initial guess for the roots
        z = cos(pi*(arth(1,1,m)-0.25d0)/(n+0.5d0))

        ! initialized unfinished
        unfinished = .true.

        ! iterate Newton steps up to maximum iterations
        do its = 1, maxits

            ! calculate Legendre polynomial at z where root has not yet been found

            ! initialize p1 and p2
            where (unfinished)
                p1 = 1d0
                p2 = 0d0
            endwhere

            ! calculate polynomial value at z by recursive formula
            do j = 1, n

                ! only where root has not yet been found
                where (unfinished)

                    ! the polynomial of order n - 2
                    p3 = p2

                    ! the polynomial of order n - 1
                    p2 = p1

                    ! the legendre polynomial
                    p1 = ((2d0*j-1d0)*z*p2-(j-1d0)*p3)/j
                endwhere
            enddo

            ! calculate derivative of polynomial p1 at z
            where (unfinished)

                ! derivative
                pp = n*(z*p1-p2)/(z*z-1d0)

                ! store old z
                z1 = z

                ! perform the newton step
                z = z1-p1/pp

                ! check for difference between old and new guess being small enough
                unfinished=(abs(z-z1) > EPS)
            endwhere

            ! if all values have sufficiently converged, stop iteration
            if (.not. any(unfinished)) exit
        end do

        ! throw error message if not sufficiently convergerd
        if(its == maxits+1)call error('legendre', 'too many iterations')

        ! else calculate abscissas
        x(1:m) = xm-xl*z

        ! symmetry for abscissas
        x(n:n-m+1:-1) = xm+xl*z

        ! calculate weights
        w(1:m) = 2d0*xl/((1d0-z**2)*pp**2)

        ! symmetry for weights
        w(n:n-m+1:-1) = w(1:m)

    end subroutine legendre


    !##############################################################################
    ! FUNCTION arth
    !
    ! Calculates incremented array from first with n entries.
    !##############################################################################
    function arth(first, increment, n)

        integer, intent(in) :: first, increment, n
        integer, parameter :: npar_arth = 16
        integer, parameter :: npar2_arth = 8
        integer :: arth(n)
        integer :: k, k2, temp

        ! initialize first element
        if(n > 0)arth(1) = first

        ! calculate by hand if n <= 16
        if(n <= npar_arth) then
            do k = 2, n
                arth(k) = arth(k-1) + increment
            enddo

        ! else set entries stepwise by 8 steps
        else
            do k = 2, npar2_arth
                arth(k) = arth(k-1) + increment
            enddo
            temp = increment*npar2_arth
            k = npar2_arth
            do
                if(k >= n)exit
                k2 = k+k
                arth(k+1:min(k2,n)) = temp+arth(1:min(k,n-k))
                temp = temp + temp
                k = k2
            enddo
        endif

    end function arth















!##############################################################################
!##############################################################################
! MODULE AR_discrete
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE discretize_AR
    !
    ! Discretizes an AR(1) process of the form z_j = \rho*z_{j-1} + eps.
    !##############################################################################
    subroutine discretize_AR(rho, mu, sigma_eps, z, pi, w)

       implicit none


       !##### INPUT/OUTPUT VARIABLES #############################################

       ! autoregression parameter
       real*8, intent(in) :: rho

       ! unconditional mean of the process
       real*8, intent(in) :: mu

       ! variance of the shock
       real*8, intent(in) :: sigma_eps

       ! discrete shock values
       real*8, intent(out) :: z(:)

       ! transition matrix
       real*8, intent(out) :: pi(:, :)

       ! the stationary distribution
       real*8, intent(out), optional :: w(:)


       !##### OTHER VARIABLES ####################################################

       integer :: n, in
       real*8 :: psi, sigma_eta


       !##### ROUTINE CODE #######################################################

       ! assert size equality and get approximation points
       n = assert_eq(size(z), size(pi,1), size(pi,2), 'discretize_AR')

       ! calculate variance of the overall process
       sigma_eta = sigma_eps/(1d0-rho**2)

       ! determine the transition matrix
       call rouwenhorst_matrix(rho, pi)

       ! determine the nodes
       psi = sqrt(dble(n-1))*sqrt(sigma_eta)
       do in = 1, n
           z(in) = -psi + 2d0*psi*dble(in-1)/dble(n-1)
       enddo
       z = z + mu

       if(present(w))then
           w = 1d0/dble(n)
           do in = 1, 10000
               w = matmul(transpose(pi), w)
           enddo
       endif

       !##########################################################################
       ! Subroutines and functions
       !##########################################################################

       contains


       !##########################################################################
       ! subroutine rouwenhorst_matrix
       !
       ! Calculates value of function that should be integrated for pis.
       !##########################################################################
       recursive subroutine rouwenhorst_matrix(rho, pi_new)

           implicit none
           real*8, intent(in) :: rho
           real*8, intent(out) :: pi_new(:, :)
           integer :: n
           real*8 :: p, pi_old(size(pi_new,1)-1, size(pi_new,1)-1)

           n = size(pi_new, 1)
           p = (1d0 + rho)/2d0

           if(n == 2)then
               pi_new(1, :) = (/p, 1d0-p/)
               pi_new(2, :) = (/1d0-p, p/)
           else
               call rouwenhorst_matrix(rho, pi_old)
               pi_new = 0d0

               pi_new(1:n-1, 1:n-1) = pi_new(1:n-1, 1:n-1) + p*pi_old
               pi_new(1:n-1, 2:n  ) = pi_new(1:n-1, 2:n  ) + (1d0-p)*pi_old
               pi_new(2:n  , 1:n-1) = pi_new(2:n  , 1:n-1) + (1d0-p)*pi_old
               pi_new(2:n  , 2:n  ) = pi_new(2:n  , 2:n  ) + p*pi_old

               pi_new(2:n-1, :) = pi_new(2:n-1, :)/2d0
           endif
       end subroutine

    end subroutine discretize_AR


    !##############################################################################
    ! SUBROUTINE discretize_log_AR
    !
    ! Discretizes a log-AR(1) process
    !##############################################################################
    subroutine discretize_log_AR(rho, mu, sigma_eps, z, pi, w)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! autoregression parameter
        real*8, intent(in) :: rho

        ! unconditional mean of the process
        real*8, intent(in) :: mu

        ! variance of the shock
        real*8, intent(in) :: sigma_eps

        ! discrete shock values
        real*8, intent(out) :: z(:)

        ! transition matrix
        real*8, intent(out) :: pi(:, :)

        ! the stationary distribution
        real*8, intent(out), optional :: w(:)


        !##### OTHER VARIABLES ####################################################

        real*8 :: sigma_eta, mu_c, sigma_c


        !##### ROUTINE CODE #######################################################

        ! calculate variance of the overall process
        sigma_eta = sigma_eps/(1d0-rho**2)

        ! get the transformed variance and expectation
        sigma_c = log(1d0+sigma_eta/mu**2)
        mu_c  = log(mu)-0.5d0*sigma_c
        sigma_c = sigma_c*(1d0-rho**2)

        ! discretize the log distribution
        if(present(w))then
            call discretize_AR(rho, mu_c, sigma_c, z, pi, w)
        else
            call discretize_AR(rho, mu_c, sigma_c, z, pi)
        endif

        ! take exponentials
        z = exp(z)

    end subroutine discretize_log_AR


    !##############################################################################
    ! SUBROUTINE simulate_AR
    !
    ! Simulates a discrete AR(1) process.
    !##############################################################################
    subroutine simulate_AR(pi, shocks)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! transition matrix
        real*8, intent(in) :: pi(:, :)

        ! simulated schocks
        integer, intent(out) :: shocks(:)


        !##### OTHER VARIABLES ####################################################

        integer :: T, n, j


        !##### ROUTINE CODE #######################################################

        ! assert size equality and get number of simulated schocks
        n = assert_eq(size(pi,1), size(pi,2), 'tauchen')
        T = size(shocks)

        if(tbox_seed)then
            call init_random_seed()
            tbox_seed = .false.
        endif

        ! get first entry
        shocks(1) = n/2+1

        ! now calculate other shocks
        do j = 2, T
            shocks(j) = get_tomorrow(pi(shocks(j-1), :))
        enddo


    !##########################################################################
    ! Subroutines and functions
    !##########################################################################

    contains


        !##########################################################################
        ! FUNCTION get_tomorrow
        !
        ! Calculates value of function that should be integrated for pis.
        !##########################################################################
        function get_tomorrow(pi)

            implicit none


            !##### INPUT/OUTPUT VARIABLES #########################################

            ! transition probabilities
            real*8, intent(in) :: pi(:)

            ! tomorrows shock
            integer :: get_tomorrow


            !##### OTHER VARIABLES ################################################

            real*8 :: rand
            integer :: i1


            !##### ROUTINE CODE ###################################################

            ! get random number
            call random_number(rand)

            ! get tomorrows value
            do i1 = 1, size(pi, 1)-1

                if(rand <= sum(pi(1:i1), 1))then
                    get_tomorrow = i1
                    return
                endif
            enddo

            ! else choose last value
            get_tomorrow = i1
            return

        end function

    end subroutine simulate_AR












!##############################################################################
!##############################################################################
! MODULE gnuplot
!##############################################################################
!##############################################################################


    !##############################################################################
    ! SUBROUTINE execplot
    !
    ! Actually creates the plot files.
    !##############################################################################
    subroutine execplot(xlim, xticks, xlabel, ylim, yticks, ylabel, title, &
            legend, filename, filetype, output)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        ! x axis definitial
        real*8, optional :: xlim(2)

        ! x axis tick definitions
        real*8, optional :: xticks

        ! y axis definitial
        real*8, optional :: ylim(2)

        ! y axis tick definitions
        real*8, optional :: yticks

        ! output file name
        character(LEN=*), optional :: title

        ! output file name
        character(LEN=*), optional :: xlabel

        ! output file name
        character(LEN=*), optional :: ylabel

        ! legend position
        character(LEN=2), optional :: legend

        ! output file name
        character(LEN=*), optional :: filename

        ! file type
        character(LEN=*), optional :: filetype

        ! output file name
        character(LEN=*), optional :: output

        !##### OTHER VARIABLES ####################################################

        integer :: i1, i2
        character(LEN=3) :: ft
        character(LEN=150) :: cfile, dfile


        !##### ROUTINE CODE #######################################################

        if(present(output))then
            cfile = output//'_c.dat'
            dfile = output//'_d.dat'
        else
            cfile = 'command13545431.dat'
            dfile = 'plotdata13545431.dat'
        endif

        ! write the output data file
        open(213659,file=trim(dfile))
        do i1 = 1, gnu_nmax
            write(213659, '(2000e21.10e4)')(gnu_x(i1, i2), gnu_y(i1, i2), i2=1, gnu_mmax)
        enddo
        close(213659)

        ! write the command file
        open(213659,file=trim(cfile))

        ! set terminal and grid
        write(213659, '(a)')'set terminal wxt title  "Gnuplot"'
        write(213659, '(a)')'set grid'
        if(gnu_histogram)then
            write(213659, '(a)')'set style data histograms'
            write(213659, '(a)')'set style fill solid border -1'
        endif

        ! set x axis
        if(present(xlim))write(213659, '(a,e13.5,a,e13.5,a)')'set xrange [',minval(xlim),':',maxval(xlim),']'
        if(present(xticks))write(213659, '(a,e13.5)')'set xtics ',xticks
        if(present(xlabel))write(213659, '(a)')'set xlabel "'//xlabel//'"font ",12"'

        ! set y axis
        if(present(ylim))write(213659, '(a,e13.5,a,e13.5,a)')'set yrange [',minval(ylim),':',maxval(ylim),']'
        if(present(yticks))write(213659, '(a,e13.5)')'set ytics ',yticks
        if(present(ylabel))write(213659, '(a)')'set ylabel "'//ylabel//'"font ",12"'

        ! set title
        if(present(title))write(213659, '(a)')'set title "'//title//'" font ",16"'

        ! legend statement
        if(gnu_dolegend)then
            if(present(legend))then
                select case (legend(1:2))
                    case("ln")
                        write(213659, '(a)')'set key left top'
                    case("ls")
                        write(213659, '(a)')'set key left bottom'
                    case("lo")
                        write(213659, '(a)')'set key left outside'
                    case("lb")
                        write(213659, '(a)')'set key left below'
                    case("rn")
                        write(213659, '(a)')'set key right top'
                    case("rs")
                        write(213659, '(a)')'set key right bottom'
                    case("ro")
                        write(213659, '(a)')'set key right outside'
                    case("rb")
                        write(213659, '(a)')'set key right below'
                    case("cn")
                        write(213659, '(a)')'set key center top'
                    case("cs")
                        write(213659, '(a)')'set key center bottom'
                    case("co")
                        write(213659, '(a)')'set key center outside'
                    case("cb")
                        write(213659, '(a)')'set key center below'
                    case default
                        write(213659, '(a)')'set key center below'
                end select
            else
                write(213659, '(a)')'set key center below'
            endif
        else
            write(213659, '(a)')'unset key'
        endif

        ! write plot lines
        if(gnu_mmax == 1)then
            write(213659, '(a)')'plot "'//trim(dfile)//'" using 1:2 '//trim(gnu_definitions(1))
        else
            write(213659, '(a)')'plot "'//trim(dfile)//'" using 1:2 '//trim(gnu_definitions(1))//',\'
            do i1 = 2, gnu_mmax-1
                write(213659, '(a,i4,a,i4,a)')'     "'//trim(dfile)//'" using '&
                    ,2*(i1-1)+1,':',2*(i1-1)+2,' '//trim(gnu_definitions(i1))//',\'
            enddo
            write(213659, '(a,i4,a,i4,a)')'     "'//trim(dfile)//'" using ',  &
                2*(gnu_mmax-1)+1,':',2*(gnu_mmax-1)+2,' '//trim(gnu_definitions(i1))
        endif

        write(213659, '(a)')'pause -1 "Press RETURN to continue..."'

        ! write graph to file
        if(present(filename))then
            ft = "eps"
            if(present(filetype))then
                if(filetype(1:3) == "png")ft= "png"
            endif
            write(213659, *)
            write(213659, '(a)')'set terminal '//ft
            write(213659, '(a)')'set output "'//filename//'.'//ft//'"'
            write(213659, '(a)')'replot'
        endif

        write(213659, '(a)')'q'
        close(213659)

        call system('gnuplot "'//trim(cfile)//'"')
        if(.not.present(output))then
            open(213659,file=trim(dfile))
            close(213659, status='delete')
            open(213659,file=trim(cfile))
            close(213659, status='delete')
        endif

        gnu_addtoplot = .false.
        gnu_dolegend = .false.
        gnu_histogram = .false.

    end subroutine


    !##############################################################################
    ! SUBROUTINE plot
    !
    ! Plots a x-y-data column to the output file.
    !##############################################################################
    subroutine plot(xin, yin, color, linewidth, marker, markersize, noline, legend)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        real*8, intent(in) :: xin(:), yin(:)
        character(LEN=*), optional :: color
        real*8, optional :: linewidth
        integer, optional :: marker
        real*8, optional :: markersize
        logical, optional :: noline
        character(LEN=*), optional :: legend


        !##### OTHER VARIABLES ####################################################

        integer :: n, i1
        logical :: lines, points


        !##### ROUTINE CODE #######################################################

        n = assert_eq(size(xin, 1), size(yin, 1), 'plot')

        ! generate the respective plot data
        if(.not.gnu_addtoplot)then

            ! allocate new x array
            if(allocated(gnu_x))deallocate(gnu_x)
            allocate(gnu_x(n, 1))
            gnu_x(:, 1) = xin

            ! allocate new y array
            if(allocated(gnu_y))deallocate(gnu_y)
            allocate(gnu_y(n, 1))
            gnu_y(:, 1) = yin

            ! allocate new x_temp array
            if(allocated(gnu_x_temp))deallocate(gnu_x_temp)
            allocate(gnu_x_temp(n, 1))
            gnu_x_temp(:, 1) = xin

            ! allocate new y_temp array
            if(allocated(gnu_y_temp))deallocate(gnu_y_temp)
            allocate(gnu_y_temp(n, 1))
            gnu_y_temp(:, 1) = yin

            gnu_addtoplot = .true.
            gnu_nmax = n
            gnu_mmax = 1

        else

            ! get new number of lines
            gnu_mmax = gnu_mmax+1
            if(gnu_mmax > 1000)then
                write(*,'(/a/)')'SORRY: I CANNOT PLOT MORE THAN 1000 LINES'
                return
            endif

            ! deallocate arrays
            deallocate(gnu_x)
            deallocate(gnu_y)

            ! if the new array is the longer one
            if(n > gnu_nmax)then
                allocate(gnu_x(n, gnu_mmax))
                allocate(gnu_y(n, gnu_mmax))

                gnu_x(1:gnu_nmax, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(1:gnu_nmax, 1:gnu_mmax-1) = gnu_y_temp

                ! fill up with the same values
                do i1 = gnu_nmax+1, n
                    gnu_x(i1, 1:gnu_mmax-1) = gnu_x(gnu_nmax, 1:gnu_mmax-1)
                    gnu_y(i1, 1:gnu_mmax-1) = gnu_y(gnu_nmax, 1:gnu_mmax-1)
                enddo

                gnu_x(:, gnu_mmax) = xin
                gnu_y(:, gnu_mmax) = yin

                ! set new nmax
                gnu_nmax = n

            ! if the old array is the longer one
            else

                allocate(gnu_x(gnu_nmax, gnu_mmax))
                allocate(gnu_y(gnu_nmax, gnu_mmax))

                gnu_x(:, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(:, 1:gnu_mmax-1) = gnu_y_temp

                gnu_x(1:n, gnu_mmax) = xin
                gnu_y(1:n, gnu_mmax) = yin

                ! fill up with same values
                do i1 = n+1, gnu_nmax
                    gnu_x(i1, gnu_mmax) = gnu_x(n, gnu_mmax)
                    gnu_y(i1, gnu_mmax) = gnu_y(n, gnu_mmax)
                enddo
            endif

            deallocate(gnu_x_temp)
            allocate(gnu_x_temp(size(gnu_x,1), size(gnu_x,2)))
            gnu_x_temp = gnu_x

            deallocate(gnu_y_temp)
            allocate(gnu_y_temp(size(gnu_y,1), size(gnu_y,2)))
            gnu_y_temp = gnu_y

        endif

        ! check for lines and points
        lines = .true.
        if(present(noline))then
            if(noline)lines = .false.
        endif

        if(.not.lines)then
            points = .true.
        else
            points = .false.
        endif
        if(present(marker))points = .true.

        ! set up definitions
        gnu_definitions(gnu_mmax) = 'with'

        ! get lines and points
        if(lines .and. points)then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' linespoints'
        elseif(lines)then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lines'
        else
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' points'
        endif

        ! get the line color
        if(present(color))then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lc rgb "'//adjustl(trim(color))//'"'
        else
            if(gnu_mmax == 1)then
                gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lc rgb "blue"'
            else
                write(gnu_definitions(gnu_mmax), '(a,i4)')trim(gnu_definitions(gnu_mmax))//' lc ', gnu_mmax-1
            endif
        endif

        ! get line width
        if(lines)then
            if(present(linewidth)) then
                write(gnu_definitions(gnu_mmax), '(a,f8.2)')trim(gnu_definitions(gnu_mmax))//' lw ', max(linewidth, 0d0)
            else
                write(gnu_definitions(gnu_mmax), '(a,f8.2)')trim(gnu_definitions(gnu_mmax))//' lw ', 2d0
            endif
        endif

        ! get marker definition
        if(points)then

            if(present(marker)) then
                write(gnu_definitions(gnu_mmax), '(a,i2)')trim(gnu_definitions(gnu_mmax))//' pt ', min(marker, 13)
            else
                write(gnu_definitions(gnu_mmax), '(a,i2)')trim(gnu_definitions(gnu_mmax))//' pt ', 1
            endif

            if(present(markersize))then
                write(gnu_definitions(gnu_mmax), '(a,f8.2)')trim(gnu_definitions(gnu_mmax))//' ps ', max(markersize, 0d0)
            else
                write(gnu_definitions(gnu_mmax), '(a,f8.2)')trim(gnu_definitions(gnu_mmax))//' ps ', 1d0
            endif
        endif

        ! set the legend
        if(present(legend))then
            gnu_dolegend = .true.
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' title "'//adjustl(trim(legend))//'"'
        else
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' notitle '
        endif

    end subroutine


    !##############################################################################
    ! SUBROUTINE plot_hist
    !
    ! Plots a histogram.
    !##############################################################################
    subroutine plot_hist(xin, yin, color, legend)

        implicit none


        !##### INPUT/OUTPUT VARIABLES #############################################

        real*8, intent(in) :: xin(:), yin(:)
        character(LEN=*), optional :: color
        character(LEN=*), optional :: legend


        !##### OTHER VARIABLES ####################################################

        integer :: n, i1


        !##### ROUTINE CODE #######################################################

        n = assert_eq(size(xin, 1), size(yin, 1), 'plot')

        ! generate the respective plot data
        if(.not.gnu_addtoplot)then

            ! allocate new x array
            if(allocated(gnu_x))deallocate(gnu_x)
            allocate(gnu_x(n, 1))
            gnu_x(:, 1) = xin

            ! allocate new y array
            if(allocated(gnu_y))deallocate(gnu_y)
            allocate(gnu_y(n, 1))
            gnu_y(:, 1) = yin

            ! allocate new x_temp array
            if(allocated(gnu_x_temp))deallocate(gnu_x_temp)
            allocate(gnu_x_temp(n, 1))
            gnu_x_temp(:, 1) = xin

            ! allocate new y_temp array
            if(allocated(gnu_y_temp))deallocate(gnu_y_temp)
            allocate(gnu_y_temp(n, 1))
            gnu_y_temp(:, 1) = yin

            gnu_addtoplot = .true.
            gnu_nmax = n
            gnu_mmax = 1

        else

            ! get new number of lines
            gnu_mmax = gnu_mmax+1
            if(gnu_mmax > 1000)then
                write(*,'(/a/)')'SORRY: I CANNOT PLOT MORE THAN 1000 LINES'
                return
            endif

            ! deallocate arrays
            deallocate(gnu_x)
            deallocate(gnu_y)

            ! if the new array is the longer one
            if(n > gnu_nmax)then
                allocate(gnu_x(n, gnu_mmax))
                allocate(gnu_y(n, gnu_mmax))

                gnu_x(1:gnu_nmax, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(1:gnu_nmax, 1:gnu_mmax-1) = gnu_y_temp

                ! fill up with the same values
                do i1 = gnu_nmax+1, n
                    gnu_x(i1, 1:gnu_mmax-1) = gnu_x(gnu_nmax, 1:gnu_mmax-1)
                    gnu_y(i1, 1:gnu_mmax-1) = gnu_y(gnu_nmax, 1:gnu_mmax-1)
                enddo

                gnu_x(:, gnu_mmax) = xin
                gnu_y(:, gnu_mmax) = yin

                ! set new nmax
                gnu_nmax = n

            ! if the old array is the longer one
            else

                allocate(gnu_x(gnu_nmax, gnu_mmax))
                allocate(gnu_y(gnu_nmax, gnu_mmax))

                gnu_x(:, 1:gnu_mmax-1) = gnu_x_temp
                gnu_y(:, 1:gnu_mmax-1) = gnu_y_temp

                gnu_x(1:n, gnu_mmax) = xin
                gnu_y(1:n, gnu_mmax) = yin

                ! fill up with same values
                do i1 = n+1, gnu_nmax
                    gnu_x(i1, gnu_mmax) = gnu_x(n, gnu_mmax)
                    gnu_y(i1, gnu_mmax) = gnu_y(n, gnu_mmax)
                enddo
            endif

            deallocate(gnu_x_temp)
            allocate(gnu_x_temp(size(gnu_x,1), size(gnu_x,2)))
            gnu_x_temp = gnu_x

            deallocate(gnu_y_temp)
            allocate(gnu_y_temp(size(gnu_y,1), size(gnu_y,2)))
            gnu_y_temp = gnu_y

        endif

        ! set up definitions
        gnu_definitions(gnu_mmax) = 'with boxes'

        ! get the line color
        if(present(color))then
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' lc rgb "'//adjustl(trim(color))//'"'
        else
            write(gnu_definitions(gnu_mmax), '(a,i4)')trim(gnu_definitions(gnu_mmax))//' lc ', gnu_mmax
        endif

        ! set the legend
        if(present(legend))then
            gnu_dolegend = .true.
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' title "'//adjustl(trim(legend))//'"'
        else
            gnu_definitions(gnu_mmax) = trim(gnu_definitions(gnu_mmax))//' notitle '
        endif

        gnu_histogram = .true.

    end subroutine



end module
