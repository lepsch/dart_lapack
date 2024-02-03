      SUBROUTINE DLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMBDA, Q2, INDX, CTOT, W, S, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, LDQ, N, N1;
      double             RHO;
      // ..
      // .. Array Arguments ..
      int                CTOT( * ), INDX( * );
      double             D( * ), DLAMBDA( * ), Q( LDQ, * ), Q2( * ), S( * ), W( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D0, ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, II, IQ2, J, N12, N2, N23;
      double             TEMP;
      // ..
      // .. External Functions ..
      double             DNRM2;
      // EXTERNAL DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DLACPY, DLAED4, DLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( K.LT.0 ) {
         INFO = -1
      } else if ( N.LT.K ) {
         INFO = -2
      } else if ( LDQ.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('DLAED3', -INFO );
         RETURN
      }

      // Quick return if possible

      if (K == 0) RETURN;


      for (J = 1; J <= K; J++) { // 20
         dlaed4(K, J, DLAMBDA, W, Q( 1, J ), RHO, D( J ), INFO );

         // If the zero finder fails, the computation is terminated.

         if (INFO.NE.0) GO TO 120;
      } // 20

      if (K == 1) GO TO 110;
      if ( K == 2 ) {
         for (J = 1; J <= K; J++) { // 30
            W( 1 ) = Q( 1, J )
            W( 2 ) = Q( 2, J )
            II = INDX( 1 )
            Q( 1, J ) = W( II )
            II = INDX( 2 )
            Q( 2, J ) = W( II )
         } // 30
         GO TO 110
      }

      // Compute updated W.

      dcopy(K, W, 1, S, 1 );

      // Initialize W(I) = Q(I,I)

      dcopy(K, Q, LDQ+1, W, 1 );
      for (J = 1; J <= K; J++) { // 60
         for (I = 1; I <= J - 1; I++) { // 40
            W( I ) = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) )
         } // 40
         for (I = J + 1; I <= K; I++) { // 50
            W( I ) = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) )
         } // 50
      } // 60
      for (I = 1; I <= K; I++) { // 70
         W( I ) = SIGN( SQRT( -W( I ) ), S( I ) )
      } // 70

      // Compute eigenvectors of the modified rank-1 modification.

      for (J = 1; J <= K; J++) { // 100
         for (I = 1; I <= K; I++) { // 80
            S( I ) = W( I ) / Q( I, J )
         } // 80
         TEMP = DNRM2( K, S, 1 )
         for (I = 1; I <= K; I++) { // 90
            II = INDX( I )
            Q( I, J ) = S( II ) / TEMP
         } // 90
      } // 100

      // Compute the updated eigenvectors.

      } // 110

      N2 = N - N1
      N12 = CTOT( 1 ) + CTOT( 2 )
      N23 = CTOT( 2 ) + CTOT( 3 )

      dlacpy('A', N23, K, Q( CTOT( 1 )+1, 1 ), LDQ, S, N23 );
      IQ2 = N1*N12 + 1
      if ( N23.NE.0 ) {
         dgemm('N', 'N', N2, K, N23, ONE, Q2( IQ2 ), N2, S, N23, ZERO, Q( N1+1, 1 ), LDQ );
      } else {
         dlaset('A', N2, K, ZERO, ZERO, Q( N1+1, 1 ), LDQ );
      }

      dlacpy('A', N12, K, Q, LDQ, S, N12 );
      if ( N12.NE.0 ) {
         dgemm('N', 'N', N1, K, N12, ONE, Q2, N1, S, N12, ZERO, Q, LDQ );
      } else {
         dlaset('A', N1, K, ZERO, ZERO, Q( 1, 1 ), LDQ );
      }


      } // 120
      RETURN

      // End of DLAED3

      }
