      SUBROUTINE ZHETRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LTB, LDB, INFO
*     ..
*     .. Array Arguments ..
      int                IPIV( * ), IPIV2( * )
      COMPLEX*16         A( LDA, * ), TB( * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      int                LDTB, NB
      bool               UPPER;
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGBTRS, ZLASWP, ZTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
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
         CALL XERBLA( 'ZHETRS_AA_2STAGE', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
*     Read NB and compute LDTB
*
      NB = INT( TB( 1 ) )
      LDTB = LTB/N
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U**H*T*U.
*
         IF( N.GT.NB ) THEN
*
*           Pivot, P**T * B -> B
*
            CALL ZLASWP( NRHS, B, LDB, NB+1, N, IPIV, 1 )
*
*           Compute (U**H \ B) -> B    [ (U**H \P**T * B) ]
*
            CALL ZTRSM( 'L', 'U', 'C', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB)
*
         END IF
*
*        Compute T \ B -> B   [ T \ (U**H \P**T * B) ]
*
         CALL ZGBTRS( 'N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO)
         IF( N.GT.NB ) THEN
*
*           Compute (U \ B) -> B   [ U \ (T \ (U**H \P**T * B) ) ]
*
            CALL ZTRSM( 'L', 'U', 'N', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB)
*
*           Pivot, P * B -> B  [ P * (U \ (T \ (U**H \P**T * B) )) ]
*
            CALL ZLASWP( NRHS, B, LDB, NB+1, N, IPIV, -1 )
*
         END IF
*
      ELSE
*
*        Solve A*X = B, where A = L*T*L**H.
*
         IF( N.GT.NB ) THEN
*
*           Pivot, P**T * B -> B
*
            CALL ZLASWP( NRHS, B, LDB, NB+1, N, IPIV, 1 )
*
*           Compute (L \ B) -> B    [ (L \P**T * B) ]
*
            CALL ZTRSM( 'L', 'L', 'N', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB)
*
         END IF
*
*        Compute T \ B -> B   [ T \ (L \P**T * B) ]
*
         CALL ZGBTRS( 'N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO)
         IF( N.GT.NB ) THEN
*
*           Compute (L**H \ B) -> B   [ L**H \ (T \ (L \P**T * B) ) ]
*
            CALL ZTRSM( 'L', 'L', 'C', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB)
*
*           Pivot, P * B -> B  [ P * (L**H \ (T \ (L \P**T * B) )) ]
*
            CALL ZLASWP( NRHS, B, LDB, NB+1, N, IPIV, -1 )
*
         END IF
      END IF
*
      RETURN
*
*     End of ZHETRS_AA_2STAGE
*
      END
