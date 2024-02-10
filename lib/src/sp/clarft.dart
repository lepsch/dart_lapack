      void clarft(final int DIRECT, final int STOREV, final int N, final int K, final Matrix<double> V, final int LDV, final int TAU, final int T, final int LDT) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIRECT, STOREV;
      int                K, LDT, LDV, N;
      Complex            T( LDT, * ), TAU( * ), V( LDV, * );
      // ..

      Complex            ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      int                I, J, PREVLASTV, LASTV;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGEMV, CTRMV
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame

      // Quick return if possible

      if (N == 0) return;

      if ( lsame( DIRECT, 'F' ) ) {
         PREVLASTV = N;
         for (I = 1; I <= K; I++) {
            PREVLASTV = max( PREVLASTV, I );
            if ( TAU( I ) == ZERO ) {

               // H(i)  =  I

               for (J = 1; J <= I; J++) {
                  T[J][I] = ZERO;
               }
            } else {

               // general case

               if ( lsame( STOREV, 'C' ) ) {
                  // Skip any trailing zeros.
                  for (LASTV = N; LASTV >= I+1; LASTV--) {
                     if( V( LASTV, I ) != ZERO ) break;
                  }
                  for (J = 1; J <= I-1; J++) {
                     T[J][I] = -TAU( I ) * CONJG( V( I , J ) );
                  }
                  J = min( LASTV, PREVLASTV );

                  // T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)

                  cgemv('Conjugate transpose', J-I, I-1, -TAU( I ), V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE, T( 1, I ), 1 );
               } else {
                  // Skip any trailing zeros.
                  for (LASTV = N; LASTV >= I+1; LASTV--) {
                     if( V( I, LASTV ) != ZERO ) break;
                  }
                  for (J = 1; J <= I-1; J++) {
                     T[J][I] = -TAU( I ) * V( J , I );
                  }
                  J = min( LASTV, PREVLASTV );

                  // T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H

                  cgemm('N', 'C', I-1, 1, J-I, -TAU( I ), V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE, T( 1, I ), LDT );
               }

               // T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)

               ctrmv('Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1, I ), 1 );
               T[I][I] = TAU( I );
               if ( I > 1 ) {
                  PREVLASTV = max( PREVLASTV, LASTV );
               } else {
                  PREVLASTV = LASTV;
               }
            }
         }
      } else {
         PREVLASTV = 1;
         for (I = K; I >= 1; I--) {
            if ( TAU( I ) == ZERO ) {

               // H(i)  =  I

               for (J = I; J <= K; J++) {
                  T[J][I] = ZERO;
               }
            } else {

               // general case

               if ( I < K ) {
                  if ( lsame( STOREV, 'C' ) ) {
                     // Skip any leading zeros.
                     for (LASTV = 1; LASTV <= I-1; LASTV++) {
                        if( V( LASTV, I ) != ZERO ) break;
                     }
                     for (J = I+1; J <= K; J++) {
                        T[J][I] = -TAU( I ) * CONJG( V( N-K+I , J ) );
                     }
                     J = max( LASTV, PREVLASTV );

                     // T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)

                     cgemv('Conjugate transpose', N-K+I-J, K-I, -TAU( I ), V( J, I+1 ), LDV, V( J, I ), 1, ONE, T( I+1, I ), 1 );
                  } else {
                     // Skip any leading zeros.
                     for (LASTV = 1; LASTV <= I-1; LASTV++) {
                        if( V( I, LASTV ) != ZERO ) break;
                     }
                     for (J = I+1; J <= K; J++) {
                        T[J][I] = -TAU( I ) * V( J, N-K+I );
                     }
                     J = max( LASTV, PREVLASTV );

                     // T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H

                     cgemm('N', 'C', K-I, 1, N-K+I-J, -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, ONE, T( I+1, I ), LDT );
                  }

                  // T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)

                  ctrmv('Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 );
                  if ( I > 1 ) {
                     PREVLASTV = min( PREVLASTV, LASTV );
                  } else {
                     PREVLASTV = LASTV;
                  }
               }
               T[I][I] = TAU( I );
            }
         }
      }
      }
