c12345678901234567890123456789012345678901234567890123456789012345678901
c     $Id$
c
c     Matrix-matrix multiply routines, designed to be called from C
c     where C = A * B.  The FORTRAN equivalent is C' = B' * A'.
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
