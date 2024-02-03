      subroutine icamaxsub(n,x,incx,iamax);

      // external icamax
      int      icamax,iamax;
      int     n,incx;
      complex x(*);

      iamax=icamax(n,x,incx);
      return;
      end;
