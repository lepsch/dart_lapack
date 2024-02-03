      subroutine dasumsub(n,x,incx,asum)
c
      external dasum
      double           dasum,asum;
      int     n,incx;
      double           x(*);
c
      asum=dasum(n,x,incx)
      return
      end
