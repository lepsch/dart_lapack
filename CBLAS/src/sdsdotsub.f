      subroutine sdsdotsub(n,sb,x,incx,y,incy,dot);

      // external sdsdot
      real sb,sdsdot,dot;
      int     n,incx,incy;
      real x(*),y(*);

      dot=sdsdot(n,sb,x,incx,y,incy);
      return;
      end;
