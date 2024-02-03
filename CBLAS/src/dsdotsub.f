      subroutine dsdotsub(n,x,incx,y,incy,dot)

      // external dsdot
      double           dsdot,dot;
      int     n,incx,incy;
      real x(*),y(*)

      dot=dsdot(n,x,incx,y,incy)
      return
      end
