      subroutine dasumsub(n,x,incx,asum)

      // external dasum
      double           dasum,asum;
      int     n,incx;
      double           x(*);

      asum=dasum(n,x,incx)
      return
      end
