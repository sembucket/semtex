c12345678901234567890123456789012345678901234567890123456789012345678901
c     $Id$
c
c     Matrix-matrix, matrix-vector multiply routines,
c     designed to be called from C.
c     E.g. where C = A * B; the FORTRAN equivalent is C' = B' * A'.
c
c
      subroutine dmxm (A, nra, B, nca, C, ncb)
c
c     C = A * B.
c
      implicit none
c
      integer          nra, nca, ncb, i, j, k
      double precision A(nca, nra), B(ncb, nca), C(ncb, nra)
c
      do j = 1, ncb
         do i = 1, nra
            c(j, i) = 0.0D0
            do k = 1, nca
               c(j, i) = c(j, i) + b(j, k) * a(k, i)
            enddo
         enddo
      enddo
      return
      end
c
c
      subroutine smxm (A, nra, B, nca, C, ncb)
c
c     C = A * B.
c
      implicit none
c
      integer  nra, nca, ncb, i, j, k
      real     A(nca, nra), B(ncb, nca), C(ncb, nra)
c
      do j = 1, ncb
         do i = 1, nra
            c(j, i) = 0.0
            do k = 1, nca
               c(j, i) = c(j, i) + b(j, k) * a(k, i)
            enddo
         enddo
      enddo
      return
      end
c
c
      subroutine dmxma (A, nra, B, nca, C, ncb)
c
c     C += A * B.
c
      implicit none
c
      integer          nra, nca, ncb, i, j, k
      double precision A(nca, nra), B(ncb, nca), C(ncb, nra)
c
      do j = 1, ncb
         do i = 1, nra
            do k = 1, nca
               c(j, i) = c(j, i) + b(j, k) * a(k, i)
            enddo
         enddo
      enddo
      return
      end
c
c
      subroutine smxma (A, nra, B, nca, C, ncb)
c
c     C += A * B.
c
      implicit none
c
      integer  nra, nca, ncb, i, j, k
      real     A(nca, nra), B(ncb, nca), C(ncb, nra)
c
      do j = 1, ncb
         do i = 1, nra
            do k = 1, nca
               c(j, i) = c(j, i) + b(j, k) * a(k, i)
            enddo
         enddo
      enddo
      return
      end

c
c
c
c
      subroutine dmxms (A, nra, B, nca, C, ncb)
c
c     C -= A * B.
c
      implicit none
c
      integer          nra, nca, ncb, i, j, k
      double precision A(nca, nra), B(ncb, nca), C(ncb, nra)
c
      do j = 1, ncb
         do i = 1, nra
            do k = 1, nca
               c(j, i) = c(j, i) - b(j, k) * a(k, i)
            enddo
         enddo
      enddo
      return
      end
c
c
      subroutine smxms (A, nra, B, nca, C, ncb)
c
c     C -= A * B.
c
      implicit none
c
      integer  nra, nca, ncb, i, j, k
      real     A(nca, nra), B(ncb, nca), C(ncb, nra)
c
      do j = 1, ncb
         do i = 1, nra
            do k = 1, nca
               c(j, i) = c(j, i) - b(j, k) * a(k, i)
            enddo
         enddo
      enddo
      return
      end

c
c
      subroutine dmxv (A, nra, x, nca, y)
c
c     y = A * x.
c
      implicit none
c
      integer          nra, nca, i, j
      double precision A(nca, nra), x(nra), y(nca)
c     
c     -- Alternative ordering for testing:
c
c      do i = 1, nca
c         y(i) = 0.0D0
c      enddo
c      
c      do j = 1, nra
c         do i = 1, nca
c            y(i) = y(i) + a(j, i) * x(j)
c         enddo
c      enddo
      do i = 1, nca
         y(i) = 0.0D0
         do j = 1, nra
            y(i) = y(i) + a(j, i) * x(j)
         enddo
      enddo
      return
      end
c
c
      subroutine smxv (A, nra, x, nca, y)
c
c     y = A * x.
c
      implicit none
c
      integer nra, nca, i, j
      real    A(nca, nra), x(nra), y(nca)
c
      do i = 1, nca
         y(i) = 0.0
         do j = 1, nra
            y(i) = y(i) + a(j, i) * x(j)
         enddo
      enddo
      return
      end
