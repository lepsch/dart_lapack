      subroutine icamaxsub(n,x,incx,iamax)
c
      external icamax
      int      icamax,iamax;
      int     n,incx;
      complex x(*)
c
      iamax=icamax(n,x,incx)
      return
      end
