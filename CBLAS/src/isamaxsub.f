      subroutine isamaxsub(n,x,incx,iamax)
c
      external isamax
      int      isamax,iamax
      int     n,incx
      real x(*)
c
      iamax=isamax(n,x,incx)
      return
      end
