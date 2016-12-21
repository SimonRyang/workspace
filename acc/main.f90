program main

  ! modules
  use omp_lib

  implicit none

  integer, parameter :: numthreads = 14
	integer, parameter :: L = 1200
	integer, parameter :: N = 800
	integer, parameter :: M = 800
	real*8 :: A(N,M), B(M,L), C(N,L)
	integer :: j, i, k
	real*8 :: seconds, summ

	seconds = omp_get_wtime()

	A = 3d0
	B = 2d0
	C = 0d0

	!$omp parallel do num_threads(numthreads)
	!!$acc data copyin(A,B) copy(C)
	!!$acc kernels loop
	do j = 1, L
		do i = 1, N
			do k = 1, M
				C(i,j) = C(i,j) + A(i,k)*B(k,j)
			enddo
		enddo
	enddo
	!!$acc end data
	!$omp end parallel do

	write(*,*)'  Done!'
	seconds = omp_get_wtime() - seconds
	write(*,*) seconds	

	summ = 0
	do j = 1, L
		summ = summ + C(j, j)
	enddo

	write(*,'(a, i10)')'sum = ', int(summ)

end program
	

