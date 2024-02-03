      void clarzt(DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIRECT, STOREV;
      int                K, LDT, LDV, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            T( LDT, * ), TAU( * ), V( LDV, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLACGV, CTRMV, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      // Check for currently supported options

      INFO = 0;
      if ( !LSAME( DIRECT, 'B' ) ) {
         INFO = -1;
      } else if ( !LSAME( STOREV, 'R' ) ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('CLARZT', -INFO );
         return;
      }

      DO 20 I = K, 1, -1;
         if ( TAU( I ) == ZERO ) {

            // H(i)  =  I

            for (J = I; J <= K; J++) { // 10
               T( J, I ) = ZERO;
            } // 10
         } else {

            // general case

            if ( I < K ) {

               // T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**H

               clacgv(N, V( I, 1 ), LDV );
               cgemv('No transpose', K-I, N, -TAU( I ), V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, T( I+1, I ), 1 );
               clacgv(N, V( I, 1 ), LDV );

               // T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)

               ctrmv('Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 );
            }
            T( I, I ) = TAU( I );
         }
      } // 20
      return;

      // End of CLARZT

      }
