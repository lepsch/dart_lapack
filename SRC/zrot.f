      SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, INCY, N;
      double             C;
      COMPLEX*16         S
      // ..
      // .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
      // ..

* =====================================================================

      // .. Local Scalars ..
      int                I, IX, IY;
      COMPLEX*16         STEMP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG
      // ..
      // .. Executable Statements ..

      IF( N.LE.0 ) RETURN       IF( INCX.EQ.1 .AND. INCY.EQ.1 ) GO TO 20

      // Code for unequal increments or equal increments not equal to 1

      IX = 1
      IY = 1
      IF( INCX.LT.0 ) IX = ( -N+1 )*INCX + 1       IF( INCY.LT.0 ) IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         STEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - DCONJG( S )*CX( IX )
         CX( IX ) = STEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN

      // Code for both increments equal to 1

   20 CONTINUE
      DO 30 I = 1, N
         STEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - DCONJG( S )*CX( I )
         CX( I ) = STEMP
   30 CONTINUE
      RETURN
      }
