      void dlarft(DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIRECT, STOREV;
      int                K, LDT, LDV, N;
      // ..
      // .. Array Arguments ..
      double             T( LDT, * ), TAU( * ), V( LDV, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, PREVLASTV, LASTV;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DTRMV
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N == 0) return;

      if ( LSAME( DIRECT, 'F' ) ) {
         PREVLASTV = N;
         for (I = 1; I <= K; I++) {
            PREVLASTV = max( I, PREVLASTV );
            if ( TAU( I ) == ZERO ) {

               // H(i)  =  I

               for (J = 1; J <= I; J++) {
                  T( J, I ) = ZERO;
               }
            } else {

               // general case

               if ( LSAME( STOREV, 'C' ) ) {
                  // Skip any trailing zeros.
                  DO LASTV = N, I+1, -1;
                     if( V( LASTV, I ) != ZERO ) EXIT;
                  }
                  for (J = 1; J <= I-1; J++) {
                     T( J, I ) = -TAU( I ) * V( I , J );
                  }
                  J = min( LASTV, PREVLASTV );

                  // T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)

                  dgemv('Transpose', J-I, I-1, -TAU( I ), V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE, T( 1, I ), 1 );
               } else {
                  // Skip any trailing zeros.
                  DO LASTV = N, I+1, -1;
                     if( V( I, LASTV ) != ZERO ) EXIT;
                  }
                  for (J = 1; J <= I-1; J++) {
                     T( J, I ) = -TAU( I ) * V( J , I );
                  }
                  J = min( LASTV, PREVLASTV );

                  // T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T

                  dgemv('No transpose', I-1, J-I, -TAU( I ), V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE, T( 1, I ), 1 );
               }

               // T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)

               dtrmv('Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1, I ), 1 );
               T( I, I ) = TAU( I );
               if ( I > 1 ) {
                  PREVLASTV = max( PREVLASTV, LASTV );
               } else {
                  PREVLASTV = LASTV;
               }
            }
         }
      } else {
         PREVLASTV = 1;
         DO I = K, 1, -1;
            if ( TAU( I ) == ZERO ) {

               // H(i)  =  I

               for (J = I; J <= K; J++) {
                  T( J, I ) = ZERO;
               }
            } else {

               // general case

               if ( I < K ) {
                  if ( LSAME( STOREV, 'C' ) ) {
                     // Skip any leading zeros.
                     for (LASTV = 1; LASTV <= I-1; LASTV++) {
                        if( V( LASTV, I ) != ZERO ) EXIT;
                     }
                     for (J = I+1; J <= K; J++) {
                        T( J, I ) = -TAU( I ) * V( N-K+I , J );
                     }
                     J = max( LASTV, PREVLASTV );

                     // T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)

                     dgemv('Transpose', N-K+I-J, K-I, -TAU( I ), V( J, I+1 ), LDV, V( J, I ), 1, ONE, T( I+1, I ), 1 );
                  } else {
                     // Skip any leading zeros.
                     for (LASTV = 1; LASTV <= I-1; LASTV++) {
                        if( V( I, LASTV ) != ZERO ) EXIT;
                     }
                     for (J = I+1; J <= K; J++) {
                        T( J, I ) = -TAU( I ) * V( J, N-K+I );
                     }
                     J = max( LASTV, PREVLASTV );

                     // T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T

                     dgemv('No transpose', K-I, N-K+I-J, -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, ONE, T( I+1, I ), 1 );
                  }

                  // T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)

                  dtrmv('Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 );
                  if ( I > 1 ) {
                     PREVLASTV = min( PREVLASTV, LASTV );
                  } else {
                     PREVLASTV = LASTV;
                  }
               }
               T( I, I ) = TAU( I );
            }
         }
      }
      return;

      // End of DLARFT

      }
