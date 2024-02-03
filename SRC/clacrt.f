      SUBROUTINE CLACRT( N, CX, INCX, CY, INCY, C, S )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, INCY, N;
      COMPLEX            C, S
      // ..
      // .. Array Arguments ..
      COMPLEX            CX( * ), CY( * )
      // ..

* =====================================================================

      // .. Local Scalars ..
      int                I, IX, IY;
      COMPLEX            CTEMP
      // ..
      // .. Executable Statements ..

      if (N.LE.0) RETURN       IF( INCX == 1 && INCY == 1 ) GO TO 20;

      // Code for unequal increments or equal increments not equal to 1

      IX = 1
      IY = 1
      if (INCX.LT.0) IX = ( -N+1 )*INCX + 1       IF( INCY.LT.0 ) IY = ( -N+1 )*INCY + 1;
      for (I = 1; I <= N; I++) { // 10
         CTEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - S*CX( IX )
         CX( IX ) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
      } // 10
      RETURN

      // Code for both increments equal to 1

      } // 20
      for (I = 1; I <= N; I++) { // 30
         CTEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - S*CX( I )
         CX( I ) = CTEMP
      } // 30
      RETURN
      }
