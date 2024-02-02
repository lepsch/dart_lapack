      subroutine ddotsub(n,x,incx,y,incy,dot)
c
      external ddot
      double precision ddot
      integer n,incx,incy
      double precision x(*),y(*),dot
c
      dot=ddot(n,x,incx,y,incy)
      return
      end