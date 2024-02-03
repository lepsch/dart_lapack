      subroutine dasumsub(n,x,incx,asum)
c
      external dasum
      double precision dasum,asum
      int     n,incx
      double precision x(*)
c
      asum=dasum(n,x,incx)
      return
      end
