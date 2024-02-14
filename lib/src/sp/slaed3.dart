      void slaed3(final int K, final int N, final int N1, final int D, final Matrix<double> Q_, final int LDQ, final int RHO, final int DLAMBDA, final int Q2, final int INDX, final int CTOT, final int W, final int S, final Box<int> INFO,) {
  final Q = Q_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDQ, N, N1;
      double               RHO;
      int                CTOT( * ), INDX( * );
      double               D( * ), DLAMBDA( * ), Q( LDQ, * ), Q2( * ), S( * ), W( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, II, IQ2, J, N12, N2, N23;
      double               TEMP;
      // ..
      // .. External Functions ..
      //- REAL               SNRM2;
      // EXTERNAL SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMM, SLACPY, SLAED4, SLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT

      // Test the input parameters.

      INFO = 0;

      if ( K < 0 ) {
         INFO = -1;
      } else if ( N < K ) {
         INFO = -2;
      } else if ( LDQ < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('SLAED3', -INFO );
         return;
      }

      // Quick return if possible

      if (K == 0) return;

      for (J = 1; J <= K; J++) { // 20
         slaed4(K, J, DLAMBDA, W, Q( 1, J ), RHO, D( J ), INFO );

         // If the zero finder fails, the computation is terminated.

         if (INFO != 0) GO TO 120;
      } // 20

      if (K == 1) GO TO 110;
      if ( K == 2 ) {
         for (J = 1; J <= K; J++) { // 30
            W[1] = Q( 1, J );
            W[2] = Q( 2, J );
            II = INDX( 1 );
            Q[1][J] = W( II );
            II = INDX( 2 );
            Q[2][J] = W( II );
         } // 30
         GO TO 110;
      }

      // Compute updated W.

      scopy(K, W, 1, S, 1 );

      // Initialize W(I) = Q(I,I)

      scopy(K, Q, LDQ+1, W, 1 );
      for (J = 1; J <= K; J++) { // 60
         for (I = 1; I <= J - 1; I++) { // 40
            W[I] = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) );
         } // 40
         for (I = J + 1; I <= K; I++) { // 50
            W[I] = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) );
         } // 50
      } // 60
      for (I = 1; I <= K; I++) { // 70
         W[I] = sign( sqrt( -W( I ) ), S( I ) );
      } // 70

      // Compute eigenvectors of the modified rank-1 modification.

      for (J = 1; J <= K; J++) { // 100
         for (I = 1; I <= K; I++) { // 80
            S[I] = W( I ) / Q( I, J );
         } // 80
         TEMP = SNRM2( K, S, 1 );
         for (I = 1; I <= K; I++) { // 90
            II = INDX( I );
            Q[I][J] = S( II ) / TEMP;
         } // 90
      } // 100

      // Compute the updated eigenvectors.

      } // 110

      N2 = N - N1;
      N12 = CTOT( 1 ) + CTOT( 2 );
      N23 = CTOT( 2 ) + CTOT( 3 );

      slacpy('A', N23, K, Q( CTOT( 1 )+1, 1 ), LDQ, S, N23 );
      IQ2 = N1*N12 + 1;
      if ( N23 != 0 ) {
         sgemm('N', 'N', N2, K, N23, ONE, Q2( IQ2 ), N2, S, N23, ZERO, Q( N1+1, 1 ), LDQ );
      } else {
         slaset('A', N2, K, ZERO, ZERO, Q( N1+1, 1 ), LDQ );
      }

      slacpy('A', N12, K, Q, LDQ, S, N12 );
      if ( N12 != 0 ) {
         sgemm('N', 'N', N1, K, N12, ONE, Q2, N1, S, N12, ZERO, Q, LDQ );
      } else {
         slaset('A', N1, K, ZERO, ZERO, Q( 1, 1 ), LDQ );
      }


      } // 120
      }
