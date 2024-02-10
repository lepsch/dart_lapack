      void zlarzt(DIRECT, STOREV, N, K, final Matrix<double> V, final int LDV, TAU, T, final int LDT) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIRECT, STOREV;
      int                K, LDT, LDV, N;
      Complex         T( LDT, * ), TAU( * ), V( LDV, * );
      // ..

      Complex         ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      int                I, INFO, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV, ZLACGV, ZTRMV
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame

      // Check for currently supported options

      INFO = 0;
      if ( !lsame( DIRECT, 'B' ) ) {
         INFO = -1;
      } else if ( !lsame( STOREV, 'R' ) ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('ZLARZT', -INFO );
         return;
      }

      for (I = K; I >= 1; I--) { // 20
         if ( TAU( I ) == ZERO ) {

            // H(i)  =  I

            for (J = I; J <= K; J++) { // 10
               T[J][I] = ZERO;
            } // 10
         } else {

            // general case

            if ( I < K ) {

               // T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**H

               zlacgv(N, V( I, 1 ), LDV );
               zgemv('No transpose', K-I, N, -TAU( I ), V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, T( I+1, I ), 1 );
               zlacgv(N, V( I, 1 ), LDV );

               // T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)

               ztrmv('Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 );
            }
            T[I][I] = TAU( I );
         }
      } // 20
      }
