      SUBROUTINE ZHETRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

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
      int                K, KP, LWKMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGTSV, ZSWAP, ZTRSM, ZLACGV, ZLACPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      if ( MIN( N, NRHS ).EQ.0 ) {
         LWKMIN = 1
      } else {
         LWKMIN = 3*N-2
      }

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
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZHETRS_AA', -INFO )
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = LWKMIN
         RETURN
      }

      // Quick return if possible

      IF( MIN( N, NRHS ).EQ.0 ) RETURN

      if ( UPPER ) {

         // Solve A*X = B, where A = U**H*T*U.

         // 1) Forward substitution with U**H

         if ( N.GT.1 ) {

            // Pivot, P**T * B -> B

            DO K = 1, N
               KP = IPIV( K )
               IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END DO

            // Compute U**H \ B -> B    [ (U**H \P**T * B) ]

            CALL ZTRSM( 'L', 'U', 'C', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B( 2, 1 ), LDB )
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (U**H \P**T * B) ]

         CALL ZLACPY( 'F', 1, N, A(1, 1), LDA+1, WORK(N), 1 )
         if ( N.GT.1 ) {
             CALL ZLACPY( 'F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 2*N ), 1)
             CALL ZLACPY( 'F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 1 ), 1 )
             CALL ZLACGV( N-1, WORK( 1 ), 1 )
         }
         CALL ZGTSV( N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO )

         // 3) Backward substitution with U

         if ( N.GT.1 ) {

            // Compute U \ B -> B   [ U \ (T \ (U**H \P**T * B) ) ]

            CALL ZTRSM( 'L', 'U', 'N', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B(2, 1), LDB)

            // Pivot, P * B  [ P * (U**H \ (T \ (U \P**T * B) )) ]

            DO K = N, 1, -1
               KP = IPIV( K )
               IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END DO
         }

      } else {

         // Solve A*X = B, where A = L*T*L**H.

         // 1) Forward substitution with L

         if ( N.GT.1 ) {

            // Pivot, P**T * B -> B

            DO K = 1, N
               KP = IPIV( K )
               IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END DO

            // Compute L \ B -> B    [ (L \P**T * B) ]

            CALL ZTRSM( 'L', 'L', 'N', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B(2, 1), LDB)
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (L \P**T * B) ]

         CALL ZLACPY( 'F', 1, N, A(1, 1), LDA+1, WORK(N), 1)
         if ( N.GT.1 ) {
             CALL ZLACPY( 'F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 1 ), 1)
             CALL ZLACPY( 'F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 2*N ), 1)
             CALL ZLACGV( N-1, WORK( 2*N ), 1 )
         }
         CALL ZGTSV(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO)

         // 3) Backward substitution with L**H

         if ( N.GT.1 ) {

            // Compute L**H \ B -> B   [ L**H \ (T \ (L \P**T * B) ) ]

            CALL ZTRSM( 'L', 'L', 'C', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B( 2, 1 ), LDB)

            // Pivot, P * B  [ P * (L**H \ (T \ (L \P**T * B) )) ]

            DO K = N, 1, -1
               KP = IPIV( K )
               IF( KP.NE.K ) CALL ZSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END DO
         }

      }

      RETURN

      // End of ZHETRS_AA

      }
