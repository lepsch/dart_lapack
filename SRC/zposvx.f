      SUBROUTINE ZPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ), S( * )
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU;
      int                I, INFEQU, J
      DOUBLE PRECISION   AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM
*     ..
*     .. External Functions ..
      bool               LSAME;
      DOUBLE PRECISION   DLAMCH, ZLANHE
      EXTERNAL           LSAME, DLAMCH, ZLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACPY, ZLAQHE, ZPOCON, ZPOEQU, ZPORFS, ZPOTRF, ZPOTRS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         RCEQU = .FALSE.
      ELSE
         RCEQU = LSAME( EQUED, 'Y' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      END IF
*
*     Test the input parameters.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT. ( RCEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -9
      ELSE
         IF( RCEQU ) THEN
            SMIN = BIGNUM
            SMAX = ZERO
            DO 10 J = 1, N
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
   10       CONTINUE
            IF( SMIN.LE.ZERO ) THEN
               INFO = -10
            ELSE IF( N.GT.0 ) THEN
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            ELSE
               SCOND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -12
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -14
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPOSVX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*        Compute row and column scalings to equilibrate the matrix A.
*
         CALL ZPOEQU( N, A, LDA, S, SCOND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*           Equilibrate the matrix.
*
            CALL ZLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )
            RCEQU = LSAME( EQUED, 'Y' )
         END IF
      END IF
*
*     Scale the right hand side.
*
      IF( RCEQU ) THEN
         DO 30 J = 1, NRHS
            DO 20 I = 1, N
               B( I, J ) = S( I )*B( I, J )
   20       CONTINUE
   30    CONTINUE
      END IF
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the Cholesky factorization A = U**H *U or A = L*L**H.
*
         CALL ZLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL ZPOTRF( UPLO, N, AF, LDAF, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.GT.0 )THEN
            RCOND = ZERO
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A.
*
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
*
*     Compute the reciprocal of the condition number of A.
*
      CALL ZPOCON( UPLO, N, AF, LDAF, ANORM, RCOND, WORK, RWORK, INFO )
*
*     Compute the solution matrix X.
*
      CALL ZLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL ZPOTRS( UPLO, N, NRHS, AF, LDAF, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL ZPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
*     Transform the solution matrix X to a solution of the original
*     system.
*
      IF( RCEQU ) THEN
         DO 50 J = 1, NRHS
            DO 40 I = 1, N
               X( I, J ) = S( I )*X( I, J )
   40       CONTINUE
   50    CONTINUE
         DO 60 J = 1, NRHS
            FERR( J ) = FERR( J ) / SCOND
   60    CONTINUE
      END IF
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      RETURN
*
*     End of ZPOSVX
*
      END
