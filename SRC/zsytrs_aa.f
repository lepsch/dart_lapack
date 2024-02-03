      SUBROUTINE ZSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LDB, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      COMPLEX*16         ONE
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                K, KP, LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGTSV, ZSWAP, ZLACPY, ZTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LWORK.LT.MAX( 1, 3*N-2 ) .AND. .NOT.LQUERY ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         xerbla('ZSYTRS_AA', -INFO );
         RETURN
      } else if ( LQUERY ) {
         LWKOPT = (3*N-2)
         WORK( 1 ) = LWKOPT
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      if ( UPPER ) {

         // Solve A*X = B, where A = U**T*T*U.

         // 1) Forward substitution with U**T

         if ( N.GT.1 ) {

            // Pivot, P**T * B -> B

            for (K = 1; K <= N; K++) {
               KP = IPIV( K )
               IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            }

            // Compute U**T \ B -> B    [ (U**T \P**T * B) ]

            ztrsm('L', 'U', 'T', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B( 2, 1 ), LDB);
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (U**T \P**T * B) ]

         zlacpy('F', 1, N, A( 1, 1 ), LDA+1, WORK( N ), 1);
         if ( N.GT.1 ) {
            zlacpy('F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 1 ), 1 );
            zlacpy('F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 2*N ), 1 );
         }
         zgtsv(N, NRHS, WORK( 1 ), WORK( N ), WORK( 2*N ), B, LDB, INFO );

         // 3) Backward substitution with U

         if ( N.GT.1 ) {

            // Compute U \ B -> B   [ U \ (T \ (U**T \P**T * B) ) ]

            ztrsm('L', 'U', 'N', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B( 2, 1 ), LDB);

            // Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]

            DO K = N, 1, -1
               KP = IPIV( K )
               IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            }
         }

      } else {

         // Solve A*X = B, where A = L*T*L**T.

         // 1) Forward substitution with L

         if ( N.GT.1 ) {

            // Pivot, P**T * B -> B

            for (K = 1; K <= N; K++) {
               KP = IPIV( K )
               IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            }

            // Compute L \ B -> B    [ (L \P**T * B) ]

            ztrsm('L', 'L', 'N', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B( 2, 1 ), LDB);
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (L \P**T * B) ]

         zlacpy('F', 1, N, A(1, 1), LDA+1, WORK(N), 1);
         if ( N.GT.1 ) {
            zlacpy('F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 1 ), 1 );
            zlacpy('F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 2*N ), 1 );
         }
         zgtsv(N, NRHS, WORK( 1 ), WORK(N), WORK( 2*N ), B, LDB, INFO);

         // 3) Backward substitution with L**T

         if ( N.GT.1 ) {

            // Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]

            ztrsm('L', 'L', 'T', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B( 2, 1 ), LDB);

            // Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]

            DO K = N, 1, -1
               KP = IPIV( K )
               IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            }
         }

      }

      RETURN

      // End of ZSYTRS_AA

      }
