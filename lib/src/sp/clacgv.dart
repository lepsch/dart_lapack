      void clacgv(N, X, INCX ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, N;
      Complex            X( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IOFF;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG

      if ( INCX == 1 ) {
         for (I = 1; I <= N; I++) { // 10
            X[I] = CONJG( X( I ) );
         } // 10
      } else {
         IOFF = 1;
         if (INCX < 0) IOFF = 1 - ( N-1 )*INCX;
         for (I = 1; I <= N; I++) { // 20
            X[IOFF] = CONJG( X( IOFF ) );
            IOFF = IOFF + INCX;
         } // 20
      }
      }
