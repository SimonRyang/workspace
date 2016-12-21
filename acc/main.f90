program main

  ! modules
  !use omp_lib

  implicit none

  integer, parameter :: numthreads = 4
	integer, parameter :: L = 1200
	integer, parameter :: N = 4000
	integer, parameter :: M = 300
	real*8 :: A(N,M), B(M,L), C(N,L)
	integer :: j, i, k
	real*8 :: seconds

	seconds = omp_get_wtime()

	A = 3d0
	B = 2d0
	C = 0d0

	do j = 1, L
		do i = 1, N
			do k = 1, M
				C(i,j) = C(i,j) + A(i,k)*B(k,j)
			enddo
		enddo
	enddo

	write(*,*)'  Done!'
	seconds = seconds - omp_get_wtime()
	write(*,*) seconds	

end program
	

