      SUBROUTINE CSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N
      REAL               ANORM, RCOND
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      COMPLEX            A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               UPPER;
      int                I, KASE
      REAL               AINVNM
*     ..
*     .. Local Arrays ..
      int                ISAVE( 3 )
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACN2, CSYTRS, XERBLA
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
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSYCON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.LE.ZERO ) THEN
         RETURN
      END IF
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         DO 10 I = N, 1, -1
            IF( IPIV( I ).GT.0 .AND. A( I, I ).EQ.ZERO ) RETURN
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         DO 20 I = 1, N
            IF( IPIV( I ).GT.0 .AND. A( I, I ).EQ.ZERO ) RETURN
   20    CONTINUE
      END IF
*
*     Estimate the 1-norm of the inverse.
*
      KASE = 0
   30 CONTINUE
      CALL CLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
*
*        Multiply by inv(L*D*L**T) or inv(U*D*U**T).
*
         CALL CSYTRS( UPLO, N, 1, A, LDA, IPIV, WORK, N, INFO )
         GO TO 30
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM
*
      RETURN
*
*     End of CSYCON
*
      END
