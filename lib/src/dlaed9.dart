      void dlaed9(K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMBDA, W, S, LDS, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, K, KSTART, KSTOP, LDQ, LDS, N;
      double             RHO;
      // ..
      // .. Array Arguments ..
      double             D( * ), DLAMBDA( * ), Q( LDQ, * ), S( LDS, * ), W( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      double             TEMP;
      // ..
      // .. External Functions ..
      //- double             DNRM2;
      // EXTERNAL DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAED4, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;

      if ( K < 0 ) {
         INFO = -1;
      } else if ( KSTART < 1 || KSTART > max( 1, K ) ) {
         INFO = -2;
      } else if ( max( 1, KSTOP ) < KSTART || KSTOP > max( 1, K ) ) {
         INFO = -3;
      } else if ( N < K ) {
         INFO = -4;
      } else if ( LDQ < max( 1, K ) ) {
         INFO = -7;
      } else if ( LDS < max( 1, K ) ) {
         INFO = -12;
      }
      if ( INFO != 0 ) {
         xerbla('DLAED9', -INFO );
         return;
      }

      // Quick return if possible

      if (K == 0) return;

      for (J = KSTART; J <= KSTOP; J++) { // 20
         dlaed4(K, J, DLAMBDA, W, Q( 1, J ), RHO, D( J ), INFO );

         // If the zero finder fails, the computation is terminated.

         if (INFO != 0) GO TO 120;
      } // 20

      if ( K == 1 || K == 2 ) {
         for (I = 1; I <= K; I++) { // 40
            for (J = 1; J <= K; J++) { // 30
               S[J, I] = Q( J, I );
            } // 30
         } // 40
         GO TO 120;
      }

      // Compute updated W.

      dcopy(K, W, 1, S, 1 );

      // Initialize W(I) = Q(I,I)

      dcopy(K, Q, LDQ+1, W, 1 );
      for (J = 1; J <= K; J++) { // 70
         for (I = 1; I <= J - 1; I++) { // 50
            W[I] = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) );
         } // 50
         for (I = J + 1; I <= K; I++) { // 60
            W[I] = W( I )*( Q( I, J )/( DLAMBDA( I )-DLAMBDA( J ) ) );
         } // 60
      } // 70
      for (I = 1; I <= K; I++) { // 80
         W[I] = sign( sqrt( -W( I ) ), S( I, 1 ) );
      } // 80

      // Compute eigenvectors of the modified rank-1 modification.

      for (J = 1; J <= K; J++) { // 110
         for (I = 1; I <= K; I++) { // 90
            Q[I, J] = W( I ) / Q( I, J );
         } // 90
         TEMP = dnrm2( K, Q( 1, J ), 1 );
         for (I = 1; I <= K; I++) { // 100
            S[I, J] = Q( I, J ) / TEMP;
         } // 100
      } // 110

      } // 120
      return;
      }
