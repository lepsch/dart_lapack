      subroutine sdotsub(n,x,incx,y,incy,dot)
c
      external sdot
      real sdot
      integer n,incx,incy
      real x(*),y(*),dot
c
      dot=sdot(n,x,incx,y,incy)
      return
      end