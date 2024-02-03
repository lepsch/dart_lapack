      SUBROUTINE ZLACRT( N, CX, INCX, CY, INCY, C, S )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, INCY, N;
      COMPLEX*16         C, S
      // ..
      // .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
      // ..

* =====================================================================

      // .. Local Scalars ..
      int                I, IX, IY;
      COMPLEX*16         CTEMP
      // ..
      // .. Executable Statements ..

      IF( N.LE.0 ) RETURN       IF( INCX.EQ.1 .AND. INCY.EQ.1 ) GO TO 20

      // Code for unequal increments or equal increments not equal to 1

      IX = 1
      IY = 1
      IF( INCX.LT.0 ) IX = ( -N+1 )*INCX + 1       IF( INCY.LT.0 ) IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         CTEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - S*CX( IX )
         CX( IX ) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN

      // Code for both increments equal to 1

   20 CONTINUE
      DO 30 I = 1, N
         CTEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - S*CX( I )
         CX( I ) = CTEMP
   30 CONTINUE
      RETURN
      }
