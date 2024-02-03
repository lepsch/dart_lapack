      subroutine ddotsub(n,x,incx,y,incy,dot);

      // external ddot
      double           ddot;
      int     n,incx,incy;
      double           x(*),y(*),dot;

      dot=ddot(n,x,incx,y,incy);
      return;
      end;
