      SUBROUTINE SLARZT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIRECT, STOREV;
      int                K, LDT, LDV, N;
      // ..
      // .. Array Arguments ..
      REAL               T( LDT, * ), TAU( * ), V( LDV, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, STRMV, XERBLA
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
         xerbla('SLARZT', -INFO );
         RETURN
      }

      DO 20 I = K, 1, -1
         if ( TAU( I ).EQ.ZERO ) {

            // H(i)  =  I

            for (J = I; J <= K; J++) { // 10
               T( J, I ) = ZERO
   10       CONTINUE
         } else {

            // general case

            if ( I.LT.K ) {

               // T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**T

               sgemv('No transpose', K-I, N, -TAU( I ), V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, T( I+1, I ), 1 );

               // T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)

               strmv('Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 );
            }
            T( I, I ) = TAU( I )
         }
   20 CONTINUE
      RETURN

      // End of SLARZT

      }
