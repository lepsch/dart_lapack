      subroutine sdsdotsub(n,sb,x,incx,y,incy,dot)
c
      external sdsdot
      real sb,sdsdot,dot
      integer n,incx,incy
      real x(*),y(*)
c
      dot=sdsdot(n,sb,x,incx,y,incy)
      return
      end