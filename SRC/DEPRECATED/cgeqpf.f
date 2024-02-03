      SUBROUTINE CGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               RWORK( * );
      COMPLEX            A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, ITEMP, J, MA, MN, PVT;
      REAL               TEMP, TEMP2, TOL3Z;
      COMPLEX            AII;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQR2, CLARF, CLARFG, CSWAP, CUNM2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, CONJG, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SCNRM2, SLAMCH;
      // EXTERNAL ISAMAX, SCNRM2, SLAMCH
      // ..
      // .. Executable Statements ..

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
         xerbla('CGEQPF', -INFO );
         return;
      }

      MN = min( M, N );
      TOL3Z = sqrt(SLAMCH('Epsilon'));

      // Move initial columns up front

      ITEMP = 1;
      for (I = 1; I <= N; I++) { // 10
         if ( JPVT( I ) != 0 ) {
            if ( I != ITEMP ) {
               cswap(M, A( 1, I ), 1, A( 1, ITEMP ), 1 );
               JPVT( I ) = JPVT( ITEMP );
               JPVT( ITEMP ) = I;
            } else {
               JPVT( I ) = I;
            }
            ITEMP = ITEMP + 1;
         } else {
            JPVT( I ) = I;
         }
      } // 10
      ITEMP = ITEMP - 1;

      // Compute the QR factorization and update remaining columns

      if ( ITEMP > 0 ) {
         MA = min( ITEMP, M );
         cgeqr2(M, MA, A, LDA, TAU, WORK, INFO );
         if ( MA < N ) {
            cunm2r('Left', 'Conjugate transpose', M, N-MA, MA, A, LDA, TAU, A( 1, MA+1 ), LDA, WORK, INFO );
         }
      }

      if ( ITEMP < MN ) {

         // Initialize partial column norms. The first n elements of
         // work store the exact column norms.

         for (I = ITEMP + 1; I <= N; I++) { // 20
            RWORK( I ) = SCNRM2( M-ITEMP, A( ITEMP+1, I ), 1 );
            RWORK( N+I ) = RWORK( I );
         } // 20

         // Compute factorization

         for (I = ITEMP + 1; I <= MN; I++) { // 40

            // Determine ith pivot column and swap if necessary

            PVT = ( I-1 ) + ISAMAX( N-I+1, RWORK( I ), 1 );

            if ( PVT != I ) {
               cswap(M, A( 1, PVT ), 1, A( 1, I ), 1 );
               ITEMP = JPVT( PVT );
               JPVT( PVT ) = JPVT( I );
               JPVT( I ) = ITEMP;
               RWORK( PVT ) = RWORK( I );
               RWORK( N+PVT ) = RWORK( N+I );
            }

            // Generate elementary reflector H(i)

            AII = A( I, I );
            clarfg(M-I+1, AII, A( min( I+1, M ), I ), 1, TAU( I ) );
            A( I, I ) = AII;

            if ( I < N ) {

               // Apply H(i) to A(i:m,i+1:n) from the left

               AII = A( I, I );
               A( I, I ) = CMPLX( ONE );
               clarf('Left', M-I+1, N-I, A( I, I ), 1, CONJG( TAU( I ) ), A( I, I+1 ), LDA, WORK );
               A( I, I ) = AII;
            }

            // Update partial column norms

            for (J = I + 1; J <= N; J++) { // 30
               if ( RWORK( J ) != ZERO ) {

                  // NOTE: The following 4 lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ABS( A( I, J ) ) / RWORK( J );
                  TEMP = max( ZERO, ( ONE+TEMP )*( ONE-TEMP ) );
                  TEMP2 = TEMP*( RWORK( J ) / RWORK( N+J ) )**2;
                  if ( TEMP2 <= TOL3Z ) {
                     if ( M-I > 0 ) {
                        RWORK( J ) = SCNRM2( M-I, A( I+1, J ), 1 );
                        RWORK( N+J ) = RWORK( J );
                     } else {
                        RWORK( J ) = ZERO;
                        RWORK( N+J ) = ZERO;
                     }
                  } else {
                     RWORK( J ) = RWORK( J )*sqrt( TEMP );
                  }
               }
            } // 30

         } // 40
      }
      return;

      // End of CGEQPF

      }
