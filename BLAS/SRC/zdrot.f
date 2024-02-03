      SUBROUTINE ZDROT( N, ZX, INCX, ZY, INCY, C, S )

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, INCY, N;
      double             C, S;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         ZX( * ), ZY( * )
      // ..

* =====================================================================

      // .. Local Scalars ..
      int                I, IX, IY;
      COMPLEX*16         CTEMP
      // ..
      // .. Executable Statements ..

      IF( N.LE.0 ) RETURN
      if ( INCX.EQ.1 .AND. INCY.EQ.1 ) {

         // code for both increments equal to 1

         DO I = 1, N
            CTEMP = C*ZX( I ) + S*ZY( I )
            ZY( I ) = C*ZY( I ) - S*ZX( I )
            ZX( I ) = CTEMP
         END DO
      } else {

         // code for unequal increments or equal increments not equal
          t // o 1

         IX = 1
         IY = 1
         IF( INCX.LT.0 ) IX = ( -N+1 )*INCX + 1          IF( INCY.LT.0 ) IY = ( -N+1 )*INCY + 1
         DO I = 1, N
            CTEMP = C*ZX( IX ) + S*ZY( IY )
            ZY( IY ) = C*ZY( IY ) - S*ZX( IX )
            ZX( IX ) = CTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      }
      RETURN

      // End of ZDROT

      }
