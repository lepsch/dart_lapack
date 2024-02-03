      SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIRECT, STOREV;
      int                K, LDT, LDV, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO
      const              ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, PREVLASTV, LASTV;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMV, ZTRMV, ZGEMM
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( LSAME( DIRECT, 'F' ) ) {
         PREVLASTV = N
         for (I = 1; I <= K; I++) {
            PREVLASTV = MAX( PREVLASTV, I )
            if ( TAU( I ).EQ.ZERO ) {

               // H(i)  =  I

               for (J = 1; J <= I; J++) {
                  T( J, I ) = ZERO
               END DO
            } else {

               // general case

               if ( LSAME( STOREV, 'C' ) ) {
                  // Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  for (J = 1; J <= I-1; J++) {
                     T( J, I ) = -TAU( I ) * CONJG( V( I , J ) )
                  END DO
                  J = MIN( LASTV, PREVLASTV )

                  // T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)

                  zgemv('Conjugate transpose', J-I, I-1, -TAU( I ), V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE, T( 1, I ), 1 );
               } else {
                  // Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  for (J = 1; J <= I-1; J++) {
                     T( J, I ) = -TAU( I ) * V( J , I )
                  END DO
                  J = MIN( LASTV, PREVLASTV )

                  // T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H

                  zgemm('N', 'C', I-1, 1, J-I, -TAU( I ), V( 1, I+1 ), LDV, V( I, I+1 ), LDV, ONE, T( 1, I ), LDT );
               }

               // T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)

               ztrmv('Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1, I ), 1 );
               T( I, I ) = TAU( I )
               if ( I.GT.1 ) {
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               } else {
                  PREVLASTV = LASTV
               }
             }
         END DO
      } else {
         PREVLASTV = 1
         DO I = K, 1, -1
            if ( TAU( I ).EQ.ZERO ) {

               // H(i)  =  I

               for (J = I; J <= K; J++) {
                  T( J, I ) = ZERO
               END DO
            } else {

               // general case

               if ( I.LT.K ) {
                  if ( LSAME( STOREV, 'C' ) ) {
                     // Skip any leading zeros.
                     for (LASTV = 1; LASTV <= I-1; LASTV++) {
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     for (J = I+1; J <= K; J++) {
                        T( J, I ) = -TAU( I ) * CONJG( V( N-K+I , J ) )
                     END DO
                     J = MAX( LASTV, PREVLASTV )

                     // T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)

                     zgemv('Conjugate transpose', N-K+I-J, K-I, -TAU( I ), V( J, I+1 ), LDV, V( J, I ), 1, ONE, T( I+1, I ), 1 );
                  } else {
                     // Skip any leading zeros.
                     for (LASTV = 1; LASTV <= I-1; LASTV++) {
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     for (J = I+1; J <= K; J++) {
                        T( J, I ) = -TAU( I ) * V( J, N-K+I )
                     END DO
                     J = MAX( LASTV, PREVLASTV )

                     // T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H

                     zgemm('N', 'C', K-I, 1, N-K+I-J, -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, ONE, T( I+1, I ), LDT );
                  }

                  // T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)

                  ztrmv('Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 );
                  if ( I.GT.1 ) {
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  } else {
                     PREVLASTV = LASTV
                  }
               }
               T( I, I ) = TAU( I )
            }
         END DO
      }
      RETURN

      // End of ZLARFT

      }
