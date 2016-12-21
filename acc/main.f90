program main

  ! modules
  use omp_lib
	use accel_lib

  implicit none

  integer, parameter :: numthreads = 14
	integer, parameter :: L = 600
	integer, parameter :: N = 500
	integer, parameter :: M = 700
	real*8 :: A(N,M), B(M,L), C(L,N,M)
	integer :: j, i, k
	real*8 :: seconds

	seconds = omp_get_wtime()

	A = 3d0
	B = 2d0
	C = 0d0

	!!$omp parallel do num_threads(numthreads)
	!$acc parallel copyin(A,B) copy(C)
	!$acc loop
	do j = 1, L
		do i = 1, N
			do k = 1, M
				C(j,i,k) = sqrt(A(i,k)*B(k,j)/2d0**0.3d0)
			enddo
		enddo
	enddo
	!$acc end parallel
	!!$omp end parallel do

	write(*,*)'  Done!'
	seconds = omp_get_wtime() - seconds
	write(*,*) seconds	

	write(*,'(a, i10)')'sum = ', int(sum(C))

end program
	

