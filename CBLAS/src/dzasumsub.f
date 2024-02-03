      subroutine dzasumsub(n,x,incx,asum)

      // external dzasum
      double           dzasum,asum;
      int     n,incx;
      double complex x(*)

      asum=dzasum(n,x,incx)
      return
      end
