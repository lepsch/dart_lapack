      void zrot(final int N, final int CX, final int INCX, final int CY, final int INCY, final int C, final int S) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, INCY, N;
      double             C;
      Complex         S;
      Complex         CX( * ), CY( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IX, IY;
      Complex         STEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG

      if (N <= 0) return;
      IF( INCX == 1 && INCY == 1 ) GO TO 20;

      // Code for unequal increments or equal increments not equal to 1

      IX = 1;
      IY = 1;
      if (INCX < 0) IX = ( -N+1 )*INCX + 1;
      IF( INCY < 0 ) IY = ( -N+1 )*INCY + 1;
      for (I = 1; I <= N; I++) { // 10
         STEMP = C*CX( IX ) + S*CY( IY );
         CY[IY] = C*CY( IY ) - DCONJG( S )*CX( IX );
         CX[IX] = STEMP;
         IX = IX + INCX;
         IY = IY + INCY;
      } // 10
      return;

      // Code for both increments equal to 1

      } // 20
      for (I = 1; I <= N; I++) { // 30
         STEMP = C*CX( I ) + S*CY( I );
         CY[I] = C*CY( I ) - DCONJG( S )*CX( I );
         CX[I] = STEMP;
      } // 30
      }
