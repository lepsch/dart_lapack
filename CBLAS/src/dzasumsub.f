      subroutine dzasumsub(n,x,incx,asum)
c
      external dzasum
      double precision dzasum,asum
      integer n,incx
      double complex x(*)
c
      asum=dzasum(n,x,incx)
      return
      end