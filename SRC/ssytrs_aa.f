      SUBROUTINE SSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

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
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                K, KP, LWKMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGTSV, SSWAP, SLACPY, STRSM, XERBLA
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
         xerbla('SSYTRS_AA', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = SROUNDUP_LWORK( LWKMIN )
         RETURN
      }

      // Quick return if possible

      IF( MIN( N, NRHS ).EQ.0 ) RETURN

      if ( UPPER ) {

         // Solve A*X = B, where A = U**T*T*U.

         // 1) Forward substitution with U**T

         if ( N.GT.1 ) {

            // Pivot, P**T * B -> B

            K = 1
            DO WHILE ( K.LE.N )
               KP = IPIV( K )
               IF( KP.NE.K ) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               K = K + 1
            END DO

            // Compute U**T \ B -> B    [ (U**T \P**T * B) ]

            strsm('L', 'U', 'T', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B( 2, 1 ), LDB);
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (U**T \P**T * B) ]

         slacpy('F', 1, N, A(1, 1), LDA+1, WORK(N), 1);
         if ( N.GT.1 ) {
             slacpy('F', 1, N-1, A(1, 2), LDA+1, WORK(1), 1);
             slacpy('F', 1, N-1, A(1, 2), LDA+1, WORK(2*N), 1);
         }
         sgtsv(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO);

         // 3) Backward substitution with U

         if ( N.GT.1 ) {


            // Compute U \ B -> B   [ U \ (T \ (U**T \P**T * B) ) ]

            strsm('L', 'U', 'N', 'U', N-1, NRHS, ONE, A( 1, 2 ), LDA, B(2, 1), LDB);

            // Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]

            K = N
            DO WHILE ( K.GE.1 )
               KP = IPIV( K )
               IF( KP.NE.K ) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               K = K - 1
            END DO
         }

      } else {

         // Solve A*X = B, where A = L*T*L**T.

         // 1) Forward substitution with L

         if ( N.GT.1 ) {

            // Pivot, P**T * B -> B

            K = 1
            DO WHILE ( K.LE.N )
               KP = IPIV( K )
               IF( KP.NE.K ) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               K = K + 1
            END DO

            // Compute L \ B -> B    [ (L \P**T * B) ]

            strsm('L', 'L', 'N', 'U', N-1, NRHS, ONE, A( 2, 1), LDA, B(2, 1), LDB);
         }

         // 2) Solve with triangular matrix T

         // Compute T \ B -> B   [ T \ (L \P**T * B) ]

         slacpy('F', 1, N, A(1, 1), LDA+1, WORK(N), 1);
         if ( N.GT.1 ) {
             slacpy('F', 1, N-1, A(2, 1), LDA+1, WORK(1), 1);
             slacpy('F', 1, N-1, A(2, 1), LDA+1, WORK(2*N), 1);
         }
         sgtsv(N, NRHS, WORK(1), WORK(N), WORK(2*N), B, LDB, INFO);

         // 3) Backward substitution with L**T

         if ( N.GT.1 ) {

            // Compute L**T \ B -> B   [ L**T \ (T \ (L \P**T * B) ) ]

            strsm('L', 'L', 'T', 'U', N-1, NRHS, ONE, A( 2, 1 ), LDA, B( 2, 1 ), LDB);

            // Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]

            K = N
            DO WHILE ( K.GE.1 )
               KP = IPIV( K )
               IF( KP.NE.K ) CALL SSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               K = K - 1
            END DO
         }

      }

      RETURN

      // End of SSYTRS_AA

      }
