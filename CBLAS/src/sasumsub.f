      subroutine sasumsub(n,x,incx,asum)
c
      external sasum
      real sasum,asum
      int     n,incx;
      real x(*)
c
      asum=sasum(n,x,incx)
      return
      end
