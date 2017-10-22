module linear_solvers

contains

   subroutine tridiagonal(upper, diag, lower, x, b, n)
      !*********************************************************************************
      !**                                                                             **
      !**  This function solves a tridiagonal linear system where:                    **
      !**                                                                             **
      !**   upper - upper diagonal entries of matrix.                                 **
      !**   diag  -       diagonal entries of matrix.                                 **
      !**   lower - lower diagonal entries of matrix.                                 **
      !**   x     - solution vector.                                                  **
      !**   b     - right hand side of linear system.                                 **
      !**   n     - size of linear system.                                            **
      !**                                                                             **
      !*********************************************************************************
      implicit none
      !---------------------------------------------------------------------------------
      integer, intent(IN) :: n
      double precision, dimension(1:n) :: upper, diag, lower, b, X
      !---------------------------------------------------------------------------------
      double precision, dimension(1:n) :: gamma
      double precision :: beta
      integer :: i
      !---------------------------------------------------------------------------------
      gamma = 0; beta = 0

      if (dabs(diag(1)) < 1d-10) stop 'Rewrite equations: trisolve'

      gamma(1) = upper(1)/diag(1)
      x(1) = b(1)/diag(1)

      ! Forward substitution

      do i = 2, n
         beta = diag(i) - gamma(i - 1)*lower(i)

         if (dabs(beta) < 1d-10) stop 'trisolve failed'

         gamma(i) = upper(i)/beta
         x(i) = (b(i) - x(i - 1)*lower(i))/beta
      end do

      ! Backward substitution
      do i = n - 1, 1, -1
         x(i) = x(i) - gamma(i)*x(i + 1)
      end do

      return

   end subroutine tridiagonal

   subroutine gaussian_elimination(a, x, b, n)
      !*********************************************************************************
      !**                                                                             **
      !**  This function solves a linear system of the form ax-b by Gaussian          **
      !**  elimination. Here:                                                         **
      !**   a     - lower diagonal entries of matrix.                                 **
      !**   x     - solution vector.                                                  **
      !**   b     - right hand side of linear system.                                 **
      !**   n     - size of linear system.                                            **
      !**                                                                             **
      !*********************************************************************************
      implicit none
      !---------------------------------------------------------------------------------
      integer, intent(IN) :: n
      double precision, intent(IN), dimension(1:n, 1:n) :: a
      double precision, intent(INOUT), dimension(1:n) :: x
      double precision, intent(IN), dimension(1:n) :: b
      !---------------------------------------------------------------------------------
      double precision, dimension(1:n, 1:n + 1) :: aug
      double precision :: pmax, aux
      integer :: i, j, k, n1, ll
      !---------------------------------------------------------------------------------

      n1 = n + 1

      aug = 0.0d0

      do i = 1, n
         do j = 1, n
            aug(i, j) = a(i, j)
         end do
      end do

      aug(:, n1) = b(:)

      do i = 1, n - 1
         pmax = 0.0

         do j = i, n
            if (abs(aug(j, i)) .gt. pmax) then
               pmax = abs(aug(j, i))
               ll = j
            endif
         enddo

         if (ll .ne. i) then
            do j = i, n1
               aux = aug(i, j)
               aug(i, j) = aug(ll, j)
               aug(ll, j) = aux
            enddo
         endif

         aux = aug(i, i)

         do j = i + 1, n1
            aug(i, j) = aug(i, j)/aux
         end do

         do k = i + 1, n
            aux = aug(k, i)

            do j = i + 1, n1
               aug(k, j) = aug(k, j) - aux*aug(i, j)
            end do
         end do

      end do

      aug(n, n1) = aug(n, n1)/aug(n, n)

      do i = n - 1, 1, -1
         do j = i + 1, n
            aug(i, n1) = aug(i, n1) - aug(i, j)*aug(j, n1)
         enddo
      enddo

      x = aug(:, n1)

      return

   end subroutine gaussian_elimination

end module linear_solvers
