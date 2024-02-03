      subroutine izamaxsub(n,x,incx,iamax);

      // external izamax
      int      izamax,iamax;
      int     n,incx;
      double complex x(*);

      iamax=izamax(n,x,incx);
      return;
      end;
