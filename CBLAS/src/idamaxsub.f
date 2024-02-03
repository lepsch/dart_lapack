      subroutine idamaxsub(n,x,incx,iamax)
c
      external idamax
      int      idamax,iamax
      int     n,incx
      double           x(*);
c
      iamax=idamax(n,x,incx)
      return
      end
