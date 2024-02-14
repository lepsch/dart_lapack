      void clag2z(final int M, final int N, final Matrix<double> SA_, final int LDSA, final Matrix<double> A_, final int LDA, final Box<int> INFO,) {
  final SA = SA_.dim();
  final A = A_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDSA, M, N;
      Complex            SA( LDSA, * );
      Complex         A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;

      INFO = 0;
      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            A[I][J] = SA( I, J );
         } // 10
      } // 20
      }
