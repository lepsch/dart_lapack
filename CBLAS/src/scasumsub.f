      subroutine scasumsub(n,x,incx,asum)
c
      external scasum
      real scasum,asum
      integer n,incx
      complex x(*)
c
      asum=scasum(n,x,incx)
      return
      end