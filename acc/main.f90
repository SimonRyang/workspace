program main

  ! modules
  !use omp_lib

  implicit none

  integer, parameter :: numthreads = 4
	integer, parameter :: L = 1200
	integer, parameter :: N = 400
	integer, parameter :: M = 300
	real*8 :: A(N,M), B(M,L), C(N,L)
	integer :: j, i, k

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

end program
	

