      SUBROUTINE CPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      REAL               BERR( * ), FERR( * ), RWORK( * ), S( * )
      COMPLEX            AFP( * ), AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU;
      int                I, INFEQU, J;
      REAL               AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHP, SLAMCH
      // EXTERNAL LSAME, CLANHP, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CLAQHP, CPPCON, CPPEQU, CPPRFS, CPPTRF, CPPTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         RCEQU = .FALSE.
      ELSE
         RCEQU = LSAME( EQUED, 'Y' )
         SMLNUM = SLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      END IF
*
      // Test the input parameters.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT. ( RCEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -7
      ELSE
         IF( RCEQU ) THEN
            SMIN = BIGNUM
            SMAX = ZERO
            DO 10 J = 1, N
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
   10       CONTINUE
            IF( SMIN.LE.ZERO ) THEN
               INFO = -8
            ELSE IF( N.GT.0 ) THEN
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            ELSE
               SCOND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -10
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -12
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPPSVX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
         // Compute row and column scalings to equilibrate the matrix A.
*
         CALL CPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
            // Equilibrate the matrix.
*
            CALL CLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED )
            RCEQU = LSAME( EQUED, 'Y' )
         END IF
      END IF
*
      // Scale the right-hand side.
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
         // Compute the Cholesky factorization A = U**H * U or A = L * L**H.
*
         CALL CCOPY( N*( N+1 ) / 2, AP, 1, AFP, 1 )
         CALL CPPTRF( UPLO, N, AFP, INFO )
*
         // Return if INFO is non-zero.
*
         IF( INFO.GT.0 )THEN
            RCOND = ZERO
            RETURN
         END IF
      END IF
*
      // Compute the norm of the matrix A.
*
      ANORM = CLANHP( 'I', UPLO, N, AP, RWORK )
*
      // Compute the reciprocal of the condition number of A.
*
      CALL CPPCON( UPLO, N, AFP, ANORM, RCOND, WORK, RWORK, INFO )
*
      // Compute the solution matrix X.
*
      CALL CLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL CPPTRS( UPLO, N, NRHS, AFP, X, LDX, INFO )
*
      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.
*
      CALL CPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
*
      // Transform the solution matrix X to a solution of the original
      // system.
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
      // Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1
*
      RETURN
*
      // End of CPPSVX
*
      END
