      subroutine izamaxsub(n,x,incx,iamax)
c
      external izamax
      int      izamax,iamax;
      int     n,incx;
      double complex x(*)
c
      iamax=izamax(n,x,incx)
      return
      end
