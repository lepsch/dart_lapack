      SUBROUTINE SSYTRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LTB, LDB, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IPIV2( * );
      REAL               A( LDA, * ), TB( * ), B( LDB, * )
      // ..
*
*  =====================================================================
*
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                LDTB, NB;
      bool               UPPER;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGBTRS, SLASWP, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LTB.LT.( 4*N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTRS_AA_2STAGE', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
      // Read NB and compute LDTB
*
      NB = INT( TB( 1 ) )
      LDTB = LTB/N
*
      IF( UPPER ) THEN
*
         // Solve A*X = B, where A = U**T*T*U.
*
         IF( N.GT.NB ) THEN
*
            // Pivot, P**T * B -> B
*
            CALL SLASWP( NRHS, B, LDB, NB+1, N, IPIV, 1 )
*
            // Compute (U**T \ B) -> B    [ (U**T \P**T * B) ]
*
            CALL STRSM( 'L', 'U', 'T', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB)
*
         END IF
*
         // Compute T \ B -> B   [ T \ (U**T \P**T * B) ]
*
         CALL SGBTRS( 'N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO)
         IF( N.GT.NB ) THEN
*
            // Compute (U \ B) -> B   [ U \ (T \ (U**T \P**T * B) ) ]
*
            CALL STRSM( 'L', 'U', 'N', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB)
*
            // Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]
*
            CALL SLASWP( NRHS, B, LDB, NB+1, N, IPIV, -1 )
*
         END IF
*
      ELSE
*
         // Solve A*X = B, where A = L*T*L**T.
*
         IF( N.GT.NB ) THEN
*
            // Pivot, P**T * B -> B
*
            CALL SLASWP( NRHS, B, LDB, NB+1, N, IPIV, 1 )
*
            // Compute (L \ B) -> B    [ (L \P**T * B) ]
*
            CALL STRSM( 'L', 'L', 'N', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB)
*
         END IF
*
         // Compute T \ B -> B   [ T \ (L \P**T * B) ]
*
         CALL SGBTRS( 'N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO)
         IF( N.GT.NB ) THEN
*
            // Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]
*
            CALL STRSM( 'L', 'L', 'T', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB)
*
            // Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]
*
            CALL SLASWP( NRHS, B, LDB, NB+1, N, IPIV, -1 )
*
         END IF
      END IF
*
      RETURN
*
      // End of SSYTRS_AA_2STAGE
*
      END
