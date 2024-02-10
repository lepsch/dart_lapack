      void zhetrs_aa(final int UPLO, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Array<int> IPIV, final Matrix<double> B, final int LDB, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      String             UPLO;
      int                N, NRHS, LDA, LDB, LWORK, INFO;
      int                IPIV( * );
      Complex         A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

// =====================================================================

      Complex         ONE;
      const              ONE = 1.0 ;
      bool               LQUERY, UPPER;
      int                K, KP, LWKMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGTSV, ZSWAP, ZTRSM, ZLACGV, ZLACPY, XERBLA
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
         xerbla('ZHETRS_AA', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK[1] = LWKMIN;
         return;
      }

      // Quick return if possible

      if( min( N, NRHS ) == 0 ) return;

      if ( UPPER ) {

         // Solve A*X = B, where A = U**H*T*U.

         // 1) Forward substitution with U**H

         if ( N > 1 ) {

            // Pivot, P**T * B -> B

            for (K = 1; K <= N; K++) {
               KP = IPIV( K );
               if (KP != K) zswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }

            // Compute U**H \ B -> B    [ (U**H \P**T * B) ]

            ztrsm('L', 'U', 'C', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B( 2, 1 ), LDB );
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (U**H \P**T * B) ]

         zlacpy('F', 1, N, A(1, 1), LDA+1, WORK(N), 1 );
         if ( N > 1 ) {
             zlacpy('F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 2*N ), 1);
             zlacpy('F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 1 ), 1 );
             zlacgv(N-1, WORK( 1 ), 1 );
         }
         zgtsv(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO );

         // 3) Backward substitution with U

         if ( N > 1 ) {

            // Compute U \ B -> B   [ U \ (T \ (U**H \P**T * B) ) ]

            ztrsm('L', 'U', 'N', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B(2, 1), LDB);

            // Pivot, P * B  [ P * (U**H \ (T \ (U \P**T * B) )) ]

            for (K = N; K >= 1; K--) {
               KP = IPIV( K );
               if (KP != K) zswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

      } else {

         // Solve A*X = B, where A = L*T*L**H.

         // 1) Forward substitution with L

         if ( N > 1 ) {

            // Pivot, P**T * B -> B

            for (K = 1; K <= N; K++) {
               KP = IPIV( K );
               if (KP != K) zswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }

            // Compute L \ B -> B    [ (L \P**T * B) ]

            ztrsm('L', 'L', 'N', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B(2, 1), LDB);
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (L \P**T * B) ]

         zlacpy('F', 1, N, A(1, 1), LDA+1, WORK(N), 1);
         if ( N > 1 ) {
             zlacpy('F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 1 ), 1);
             zlacpy('F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 2*N ), 1);
             zlacgv(N-1, WORK( 2*N ), 1 );
         }
         zgtsv(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO);

         // 3) Backward substitution with L**H

         if ( N > 1 ) {

            // Compute L**H \ B -> B   [ L**H \ (T \ (L \P**T * B) ) ]

            ztrsm('L', 'L', 'C', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B( 2, 1 ), LDB);

            // Pivot, P * B  [ P * (L**H \ (T \ (L \P**T * B) )) ]

            for (K = N; K >= 1; K--) {
               KP = IPIV( K );
               if (KP != K) zswap( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB );
            }
         }

      }

      }
