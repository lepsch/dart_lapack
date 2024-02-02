      subroutine dsdotsub(n,x,incx,y,incy,dot)
c
      external dsdot
      double precision dsdot,dot
      integer n,incx,incy
      real x(*),y(*)
c
      dot=dsdot(n,x,incx,y,incy)
      return
      end