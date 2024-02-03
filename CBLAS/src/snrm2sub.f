      subroutine snrm2sub(n,x,incx,nrm2)
c
      external snrm2
      real snrm2,nrm2
      int     n,incx;
      real x(*)
c
      nrm2=snrm2(n,x,incx)
      return
      end
