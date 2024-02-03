      subroutine scasumsub(n,x,incx,asum);

      // external scasum
      real scasum,asum;
      int     n,incx;
      complex x(*);

      asum=scasum(n,x,incx);
      return;
      end;
