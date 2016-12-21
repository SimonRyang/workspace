program main

  ! modules
  use omp_lib

  implicit none

  integer, parameter :: numthreads = 28
	integer, parameter :: L = 1200
	integer, parameter :: N = 800
	integer, parameter :: M = 600
	integer :: A(N,M), B(M,L), C(N,L)
	integer :: j, i, k, sum
	real*8 :: seconds

	seconds = omp_get_wtime()

	A = 3
	B = 2
	C = 0

	!!$acc data copyin(A,B) copy(C)
	!!$acc kernels loop
	!$omp parallel do reduction(+:C) num_threads(numthreads)
	do j = 1, L
		do i = 1, N
			do k = 1, M
				C(i,j) = C(i,j) + A(i,k)*B(k,j)
			enddo
		enddo
	enddo
	!$omp end parallel do
	!!$acc end data

	write(*,*)'  Done!'
	seconds = omp_get_wtime() - seconds
	write(*,*) seconds	

	sum = 0
	do j = 1, L
		sum = sum + C(j, j)
	enddo

	write(*,'(i10)')'sum = ', sum

end program
	

