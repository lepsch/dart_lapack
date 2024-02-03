      subroutine dzasumsub(n,x,incx,asum)
c
      external dzasum
      double           dzasum,asum;
      int     n,incx
      double complex x(*)
c
      asum=dzasum(n,x,incx)
      return
      end
