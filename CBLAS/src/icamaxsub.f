      subroutine icamaxsub(n,x,incx,iamax)
c
      external icamax
      integer  icamax,iamax
      integer n,incx
      complex x(*)
c
      iamax=icamax(n,x,incx)
      return
      end