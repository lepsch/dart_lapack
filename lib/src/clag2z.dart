      void clag2z(M, N, SA, LDSA, A, LDA, Box<int> INFO ) {

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
