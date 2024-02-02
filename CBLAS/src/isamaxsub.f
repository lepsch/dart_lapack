      subroutine isamaxsub(n,x,incx,iamax)
c
      external isamax
      integer  isamax,iamax
      integer n,incx
      real x(*)
c
      iamax=isamax(n,x,incx)
      return
      end