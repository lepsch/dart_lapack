      subroutine ddotsub(n,x,incx,y,incy,dot)
c
      external ddot
      double           ddot;
      int     n,incx,incy
      double           x(*),y(*),dot;
c
      dot=ddot(n,x,incx,y,incy)
      return
      end
