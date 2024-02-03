      SUBROUTINE ZPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), FERR( * ), RWORK( * ), S( * );
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU, UPPER;
      int                I, INFEQU, J, J1, J2;
      double             AMAX, ANORM, BIGNUM, SCOND, SMAX, SMIN, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHB;
      // EXTERNAL LSAME, DLAMCH, ZLANHB
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZLACPY, ZLAQHB, ZPBCON, ZPBEQU, ZPBRFS, ZPBTRF, ZPBTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      UPPER = LSAME( UPLO, 'U' )
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         RCEQU = .FALSE.
      } else {
         RCEQU = LSAME( EQUED, 'Y' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      END IF

      // Test the input parameters.

      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -7
      ELSE IF( LDAFB.LT.KD+1 ) THEN
         INFO = -9
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT. ( RCEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -10
      } else {
         IF( RCEQU ) THEN
            SMIN = BIGNUM
            SMAX = ZERO
            DO 10 J = 1, N
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
   10       CONTINUE
            IF( SMIN.LE.ZERO ) THEN
               INFO = -11
            ELSE IF( N.GT.0 ) THEN
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            } else {
               SCOND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -13
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -15
            END IF
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPBSVX', -INFO )
         RETURN
      END IF

      IF( EQUIL ) THEN

         // Compute row and column scalings to equilibrate the matrix A.

         CALL ZPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN

            // Equilibrate the matrix.

            CALL ZLAQHB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )
            RCEQU = LSAME( EQUED, 'Y' )
         END IF
      END IF

      // Scale the right-hand side.

      IF( RCEQU ) THEN
         DO 30 J = 1, NRHS
            DO 20 I = 1, N
               B( I, J ) = S( I )*B( I, J )
   20       CONTINUE
   30    CONTINUE
      END IF

      IF( NOFACT .OR. EQUIL ) THEN

         // Compute the Cholesky factorization A = U**H *U or A = L*L**H.

         IF( UPPER ) THEN
            DO 40 J = 1, N
               J1 = MAX( J-KD, 1 )
               CALL ZCOPY( J-J1+1, AB( KD+1-J+J1, J ), 1, AFB( KD+1-J+J1, J ), 1 )
   40       CONTINUE
         } else {
            DO 50 J = 1, N
               J2 = MIN( J+KD, N )
               CALL ZCOPY( J2-J+1, AB( 1, J ), 1, AFB( 1, J ), 1 )
   50       CONTINUE
         END IF

         CALL ZPBTRF( UPLO, N, KD, AFB, LDAFB, INFO )

         // Return if INFO is non-zero.

         IF( INFO.GT.0 )THEN
            RCOND = ZERO
            RETURN
         END IF
      END IF

      // Compute the norm of the matrix A.

      ANORM = ZLANHB( '1', UPLO, N, KD, AB, LDAB, RWORK )

      // Compute the reciprocal of the condition number of A.

      CALL ZPBCON( UPLO, N, KD, AFB, LDAFB, ANORM, RCOND, WORK, RWORK, INFO )

      // Compute the solution matrix X.

      CALL ZLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL ZPBTRS( UPLO, N, KD, NRHS, AFB, LDAFB, X, LDX, INFO )

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      CALL ZPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )

      // Transform the solution matrix X to a solution of the original
      // system.

      IF( RCEQU ) THEN
         DO 70 J = 1, NRHS
            DO 60 I = 1, N
               X( I, J ) = S( I )*X( I, J )
   60       CONTINUE
   70    CONTINUE
         DO 80 J = 1, NRHS
            FERR( J ) = FERR( J ) / SCOND
   80    CONTINUE
      END IF

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1

      RETURN

      // End of ZPBSVX

      }
