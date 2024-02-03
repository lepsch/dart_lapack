      SUBROUTINE ZLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )

*  -- LAPACK computational routine --
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
      COMPLEX*16         ZERO
      const              ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV, ZLACGV, ZTRMV
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      // Check for currently supported options

      INFO = 0
      if ( .NOT.LSAME( DIRECT, 'B' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( STOREV, 'R' ) ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLARZT', -INFO );
         RETURN
      }

      DO 20 I = K, 1, -1
         if ( TAU( I ).EQ.ZERO ) {

            // H(i)  =  I

            DO 10 J = I, K
               T( J, I ) = ZERO
   10       CONTINUE
         } else {

            // general case

            if ( I.LT.K ) {

               // T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**H

               zlacgv(N, V( I, 1 ), LDV );
               zgemv('No transpose', K-I, N, -TAU( I ), V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, T( I+1, I ), 1 );
               zlacgv(N, V( I, 1 ), LDV );

               // T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)

               ztrmv('Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 );
            }
            T( I, I ) = TAU( I )
         }
   20 CONTINUE
      RETURN

      // End of ZLARZT

      }
