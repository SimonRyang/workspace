program main

  ! modules
  use omp_lib

  implicit none

  integer, parameter :: numthreads = 4
  integer, parameter :: M = 10000
  integer, parameter :: L = 2000
  integer, parameter :: N = 4000
	integer :: j, i, k
  real*8 :: A(N,M), B(M, L), C(N, L)
	real*8 :: seconds

  seconds = omp_get_wtime ( )

	A = 3d0
	B = 9d0
	C = 2d0

	!$omp parallel do num_threads(numthreads)
	
	do j = 1, L
		do i = 1, N
			do k = 1, M
				C(i,j) = C(i,j) + A(i, k)*B(k, j)
			enddo
		enddo
	enddo
	!$omp end parallel do

	write(*,*)'  Done!'

  seconds = omp_get_wtime ( ) - seconds;
	write(*,*)seconds

	call sleep(10)

end program
