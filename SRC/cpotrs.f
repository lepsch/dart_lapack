      SUBROUTINE CPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
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
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPOTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B where A = U**H *U.
*
*        Solve U**H *X = B, overwriting B with X.
*
         CALL CTRSM( 'Left', 'Upper', 'Conjugate transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL CTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A*X = B where A = L*L**H.
*
*        Solve L*X = B, overwriting B with X.
*
         CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )
*
*        Solve L**H *X = B, overwriting B with X.
*
         CALL CTRSM( 'Left', 'Lower', 'Conjugate transpose', 'Non-unit', N, NRHS, ONE, A, LDA, B, LDB )
      END IF
*
      RETURN
*
*     End of CPOTRS
*
      END
