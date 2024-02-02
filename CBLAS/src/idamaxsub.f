      subroutine idamaxsub(n,x,incx,iamax)
c
      external idamax
      integer  idamax,iamax
      integer n,incx
      double precision x(*)
c
      iamax=idamax(n,x,incx)
      return
      end