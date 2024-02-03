      subroutine dsdotsub(n,x,incx,y,incy,dot)
c
      external dsdot
      double           dsdot,dot;
      int     n,incx,incy
      real x(*),y(*)
c
      dot=dsdot(n,x,incx,y,incy)
      return
      end
