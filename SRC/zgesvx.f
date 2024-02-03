      SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, TRANS;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      String             NORM;
      int                I, INFEQU, J;
      double             AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, RPVGRW, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE, ZLANTR;
      // EXTERNAL LSAME, DLAMCH, ZLANGE, ZLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGECON, ZGEEQU, ZGERFS, ZGETRF, ZGETRS, ZLACPY, ZLAQGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      } else {
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      END IF

      // Test the input parameters.

      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT. ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -10
      } else {
         IF( ROWEQU ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 10 J = 1, N
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
   10       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -11
            ELSE IF( N.GT.0 ) THEN
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               ROWCND = ONE
            END IF
         END IF
         IF( COLEQU .AND. INFO.EQ.0 ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 20 J = 1, N
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
   20       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -12
            ELSE IF( N.GT.0 ) THEN
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               COLCND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -14
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -16
            END IF
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGESVX', -INFO )
         RETURN
      END IF

      IF( EQUIL ) THEN

         // Compute row and column scalings to equilibrate the matrix A.

         CALL ZGEEQU( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN

            // Equilibrate the matrix.

            CALL ZLAQGE( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED )
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         END IF
      END IF

      // Scale the right hand side.

      IF( NOTRAN ) THEN
         IF( ROWEQU ) THEN
            DO 40 J = 1, NRHS
               DO 30 I = 1, N
                  B( I, J ) = R( I )*B( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( COLEQU ) THEN
         DO 60 J = 1, NRHS
            DO 50 I = 1, N
               B( I, J ) = C( I )*B( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF

      IF( NOFACT .OR. EQUIL ) THEN

         // Compute the LU factorization of A.

         CALL ZLACPY( 'Full', N, N, A, LDA, AF, LDAF )
         CALL ZGETRF( N, N, AF, LDAF, IPIV, INFO )

         // Return if INFO is non-zero.

         IF( INFO.GT.0 ) THEN

            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = ZLANTR( 'M', 'U', 'N', INFO, INFO, AF, LDAF, RWORK )
            IF( RPVGRW.EQ.ZERO ) THEN
               RPVGRW = ONE
            } else {
               RPVGRW = ZLANGE( 'M', N, INFO, A, LDA, RWORK ) / RPVGRW
            END IF
            RWORK( 1 ) = RPVGRW
            RCOND = ZERO
            RETURN
         END IF
      END IF

      // Compute the norm of the matrix A and the
      // reciprocal pivot growth factor RPVGRW.

      IF( NOTRAN ) THEN
         NORM = '1'
      } else {
         NORM = 'I'
      END IF
      ANORM = ZLANGE( NORM, N, N, A, LDA, RWORK )
      RPVGRW = ZLANTR( 'M', 'U', 'N', N, N, AF, LDAF, RWORK )
      IF( RPVGRW.EQ.ZERO ) THEN
         RPVGRW = ONE
      } else {
         RPVGRW = ZLANGE( 'M', N, N, A, LDA, RWORK ) / RPVGRW
      END IF

      // Compute the reciprocal of the condition number of A.

      CALL ZGECON( NORM, N, AF, LDAF, ANORM, RCOND, WORK, RWORK, INFO )

      // Compute the solution matrix X.

      CALL ZLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL ZGETRS( TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      CALL ZGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )

      // Transform the solution matrix X to a solution of the original
      // system.

      IF( NOTRAN ) THEN
         IF( COLEQU ) THEN
            DO 80 J = 1, NRHS
               DO 70 I = 1, N
                  X( I, J ) = C( I )*X( I, J )
   70          CONTINUE
   80       CONTINUE
            DO 90 J = 1, NRHS
               FERR( J ) = FERR( J ) / COLCND
   90       CONTINUE
         END IF
      ELSE IF( ROWEQU ) THEN
         DO 110 J = 1, NRHS
            DO 100 I = 1, N
               X( I, J ) = R( I )*X( I, J )
  100       CONTINUE
  110    CONTINUE
         DO 120 J = 1, NRHS
            FERR( J ) = FERR( J ) / ROWCND
  120    CONTINUE
      END IF

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.DLAMCH( 'Epsilon' ) ) INFO = N + 1

      RWORK( 1 ) = RPVGRW
      RETURN

      // End of ZGESVX

      }
