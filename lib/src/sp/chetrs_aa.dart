      void chetrs_aa(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      String             UPLO;
      int                N, NRHS, LDA, LDB, LWORK, INFO;
      int                IPIV( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      Complex            ONE;
      const              ONE = 1.0 ;
      bool               LQUERY, UPPER;
      int                K, KP, LWKMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACPY, CLACGV, CGTSV, CSWAP, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );
      if ( min( N, NRHS ) == 0 ) {
         LWKMIN = 1;
      } else {
         LWKMIN = 3*N-2;
      }

      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('CHETRS_AA', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK[1] = SROUNDUP_LWORK( LWKMIN );
         return;
      }

      // Quick return if possible

      if( min( N, NRHS ) == 0 ) return;

      if ( UPPER ) {

         // Solve A*X = B, where A = U**H*T*U.

         // 1) Forward substitution with U**H

         if ( N > 1 ) {

            // Pivot, P**T * B -> B

            K = 1;
            while (K <= N) {
               KP = IPIV( K );
               if (KP != K) cswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               K = K + 1;
            }

            // Compute U**H \ B -> B    [ (U**H \P**T * B) ]

            ctrsm('L', 'U', 'C', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B( 2, 1 ), LDB);
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (U**H \P**T * B) ]

         clacpy('F', 1, N, A(1, 1), LDA+1, WORK(N), 1);
         if ( N > 1 ) {
             clacpy('F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 2*N ), 1);
             clacpy('F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 1 ), 1);
             clacgv(N-1, WORK( 1 ), 1 );
         }
         cgtsv(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO);

         // 3) Backward substitution with U

         if ( N > 1 ) {

            // Compute U \ B -> B   [ U \ (T \ (U**H \P**T * B) ) ]

            ctrsm('L', 'U', 'N', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B(2, 1), LDB);

            // Pivot, P * B  -> B [ P * (U \ (T \ (U**H \P**T * B) )) ]

            K = N;
            while (K >= 1) {
               KP = IPIV( K );
               if (KP != K) cswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               K = K - 1;
            }
         }

      } else {

         // Solve A*X = B, where A = L*T*L**H.

         // 1) Forward substitution with L

         if ( N > 1 ) {

            // Pivot, P**T * B -> B

            K = 1;
            while (K <= N) {
               KP = IPIV( K );
               if (KP != K) cswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               K = K + 1;
            }

            // Compute L \ B -> B    [ (L \P**T * B) ]

            ctrsm('L', 'L', 'N', 'U', N-1, NRHS, ONE, A( 2, 1), LDA, B(2, 1), LDB );
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (L \P**T * B) ]

         clacpy('F', 1, N, A(1, 1), LDA+1, WORK(N), 1);
         if ( N > 1 ) {
             clacpy('F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 1 ), 1 );
             clacpy('F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 2*N ), 1);
             clacgv(N-1, WORK( 2*N ), 1 );
         }
         cgtsv(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO);

         // 3) Backward substitution with L**H

         if ( N > 1 ) {

            // Compute (L**H \ B) -> B   [ L**H \ (T \ (L \P**T * B) ) ]

            ctrsm('L', 'L', 'C', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B( 2, 1 ), LDB );

            // Pivot, P * B -> B  [ P * (L**H \ (T \ (L \P**T * B) )) ]

            K = N;
            while (K >= 1) {
               KP = IPIV( K );
               if (KP != K) cswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
               K = K - 1;
            }
         }

      }

      }
