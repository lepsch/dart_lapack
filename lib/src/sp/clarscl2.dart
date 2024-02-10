      void clarscl2(final int M, final int N, final int D, final int X, final int LDX) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                M, N, LDX;
      Complex            X( LDX, * );
      double               D( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;

      for (J = 1; J <= N; J++) {
         for (I = 1; I <= M; I++) {
            X[I][J] = X( I, J ) / D( I );
         }
      }

      }
