      SUBROUTINE SLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMBDA, W, S, LDS, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, KSTART, KSTOP, LDQ, LDS, N;
      REAL               RHO
      // ..
      // .. Array Arguments ..
      REAL               D( * ), DLAMBDA( * ), Q( LDQ, * ), S( LDS, * ), W( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      REAL               TEMP
      // ..
      // .. External Functions ..
      REAL               SNRM2
      // EXTERNAL SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLAED4, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( K.LT.0 ) {
         INFO = -1
      } else if ( KSTART.LT.1 || KSTART.GT.MAX( 1, K ) ) {
         INFO = -2
      } else if ( MAX( 1, KSTOP ).LT.KSTART || KSTOP.GT.MAX( 1, K ) ) {
         INFO = -3
      } else if ( N.LT.K ) {
         INFO = -4
      } else if ( LDQ.LT.MAX( 1, K ) ) {
         INFO = -7
      } else if ( LDS.LT.MAX( 1, K ) ) {
         INFO = -12
      }
      if ( INFO != 0 ) {
         xerbla('SLAED9', -INFO );
         RETURN
      }

      // Quick return if possible

      if (K == 0) RETURN;

      for (J = KSTART; J <= KSTOP; J++) { // 20
         slaed4(K, J, DLAMBDA, W, Q( 1, J ), RHO, D( J ), INFO );

         // If the zero finder fails, the computation is terminated.

         if (INFO != 0) GO TO 120;
      } // 20

      if ( K == 1 || K == 2 ) {
         for (I = 1; I <= K; I++) { // 40
            for (J = 1; J <= K; J++) { // 30
               S( J, I ) = Q( J, I )
            } // 30
         } // 40
         GO TO 120
      }

      // Compute updated W.

      scopy(K, W, 1, S, 1 );

      // Initialize W(I) = Q(I,I)

      scopy(K, Q, LDQ+1, W, 1 );
      for (J = 1; J <= K; J++) { // 70
         for (I = 1; I <= J - 1; I++) { // 50
            W( I ) = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) )
         } // 50
         for (I = J + 1; I <= K; I++) { // 60
            W( I ) = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) )
         } // 60
      } // 70
      for (I = 1; I <= K; I++) { // 80
         W( I ) = SIGN( SQRT( -W( I ) ), S( I, 1 ) )
      } // 80

      // Compute eigenvectors of the modified rank-1 modification.

      for (J = 1; J <= K; J++) { // 110
         for (I = 1; I <= K; I++) { // 90
            Q( I, J ) = W( I ) / Q( I, J )
         } // 90
         TEMP = SNRM2( K, Q( 1, J ), 1 )
         for (I = 1; I <= K; I++) { // 100
            S( I, J ) = Q( I, J ) / TEMP
         } // 100
      } // 110

      } // 120
      RETURN

      // End of SLAED9

      }
