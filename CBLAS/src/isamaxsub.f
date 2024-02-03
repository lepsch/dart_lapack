      subroutine isamaxsub(n,x,incx,iamax)

      // external isamax
      int      isamax,iamax;
      int     n,incx;
      real x(*)

      iamax=isamax(n,x,incx)
      return
      end
