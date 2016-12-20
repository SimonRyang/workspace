program main

  ! modules
  use omp_lib


  implicit none

  integer, parameter :: numthreads = 4
  integer, parameter :: M = 100
  integer, parameter :: L = 200
  integer, parameter :: N = 400
	integer :: j, i, k
  real*8 :: A(N,M), B(M, L), C(N, L)

	call tic()

	A = 3d0
	B = 9d0
	C = 2d0

	do j = 1, L
		do i = 1, N
			do k = 1, M
				C(i,j) = C(i,j) + A(i, k)*B(k, j)
			enddo
		enddo
	enddo

	write(*,*)'Done!'

	call toc()

	call sleep(10)

end program
