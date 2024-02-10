      void dgeqpf(M, N, A, LDA, JPVT, TAU, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      int                JPVT( * );
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, ITEMP, J, MA, MN, PVT;
      double             AII, TEMP, TEMP2, TOL3Z;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQR2, DLARF, DLARFG, DORM2R, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      //- int                idamax;
      //- double             DLAMCH, DNRM2;
      // EXTERNAL idamax, DLAMCH, DNRM2

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DGEQPF', -INFO );
         return;
      }

      MN = min( M, N );
      TOL3Z = sqrt(dlamch('Epsilon'));

      // Move initial columns up front

      ITEMP = 1;
      for (I = 1; I <= N; I++) { // 10
         if ( JPVT( I ) != 0 ) {
            if ( I != ITEMP ) {
               dswap(M, A( 1, I ), 1, A( 1, ITEMP ), 1 );
               JPVT[I] = JPVT( ITEMP );
               JPVT[ITEMP] = I;
            } else {
               JPVT[I] = I;
            }
            ITEMP = ITEMP + 1;
         } else {
            JPVT[I] = I;
         }
      } // 10
      ITEMP = ITEMP - 1;

      // Compute the QR factorization and update remaining columns

      if ( ITEMP > 0 ) {
         MA = min( ITEMP, M );
         dgeqr2(M, MA, A, LDA, TAU, WORK, INFO );
         if ( MA < N ) {
            dorm2r('Left', 'Transpose', M, N-MA, MA, A, LDA, TAU, A( 1, MA+1 ), LDA, WORK, INFO );
         }
      }

      if ( ITEMP < MN ) {

         // Initialize partial column norms. The first n elements of
         // work store the exact column norms.

         for (I = ITEMP + 1; I <= N; I++) { // 20
            WORK[I] = dnrm2( M-ITEMP, A( ITEMP+1, I ), 1 );
            WORK[N+I] = WORK( I );
         } // 20

         // Compute factorization

         for (I = ITEMP + 1; I <= MN; I++) { // 40

            // Determine ith pivot column and swap if necessary

            PVT = ( I-1 ) + idamax( N-I+1, WORK( I ), 1 );

            if ( PVT != I ) {
               dswap(M, A( 1, PVT ), 1, A( 1, I ), 1 );
               ITEMP = JPVT( PVT );
               JPVT[PVT] = JPVT( I );
               JPVT[I] = ITEMP;
               WORK[PVT] = WORK( I );
               WORK[N+PVT] = WORK( N+I );
            }

            // Generate elementary reflector H(i)

            if ( I < M ) {
               dlarfg(M-I+1, A( I, I ), A( I+1, I ), 1, TAU( I ) );
            } else {
               dlarfg(1, A( M, M ), A( M, M ), 1, TAU( M ) );
            }

            if ( I < N ) {

               // Apply H(i) to A(i:m,i+1:n) from the left

               AII = A( I, I );
               A[I][I] = ONE;
               dlarf('LEFT', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK( 2*N+1 ) );
               A[I][I] = AII;
            }

            // Update partial column norms

            for (J = I + 1; J <= N; J++) { // 30
               if ( WORK( J ) != ZERO ) {

                  // NOTE: The following 4 lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ( A( I, J ) ).abs() / WORK( J );
                  TEMP = max( ZERO, ( ONE+TEMP )*( ONE-TEMP ) );
                  TEMP2 = TEMP*( WORK( J ) / WORK( N+J ) )**2;
                  if ( TEMP2 <= TOL3Z ) {
                     if ( M-I > 0 ) {
                        WORK[J] = dnrm2( M-I, A( I+1, J ), 1 );
                        WORK[N+J] = WORK( J );
                     } else {
                        WORK[J] = ZERO;
                        WORK[N+J] = ZERO;
                     }
                  } else {
                     WORK[J] = WORK( J )*sqrt( TEMP );
                  }
               }
            } // 30

         } // 40
      }
      }
