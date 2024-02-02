      subroutine scnrm2sub(n,x,incx,nrm2)
c
      external scnrm2
      real scnrm2,nrm2
      integer n,incx
      complex x(*)
c
      nrm2=scnrm2(n,x,incx)
      return
      end