      SUBROUTINE DSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double             RCOND, RPVGRW;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX, * ), WORK( * )       double             S( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * );
      // ..
*
*  ==================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I;
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I;
      int                CMP_ERR_I, PIV_GROWTH_I;
      PARAMETER          ( FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 )
      PARAMETER          ( RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 )
      PARAMETER          ( CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 )
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU;
      int                INFEQU, J;
      double             AMAX, BIGNUM, SMIN, SMAX, SCOND, SMLNUM;
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, DLAMCH, DLA_SYRPVGRW
      bool               LSAME;
      double             DLAMCH, DLA_SYRPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYEQUB, DSYTRF, DSYTRS, DLACPY, DLAQSY, XERBLA, DLASCL2, DSYRFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         RCEQU = .FALSE.
      ELSE
         RCEQU = LSAME( EQUED, 'Y' )
      ENDIF
*
      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in DSYRFSX.
*
      RPVGRW = ZERO
*
      // Test the input parameters.  PARAMS is not tested until DSYRFSX.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT. LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME(UPLO, 'U') .AND. .NOT.LSAME(UPLO, 'L') ) THEN
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
         INFO = -10
      ELSE
         IF ( RCEQU ) THEN
            SMIN = BIGNUM
            SMAX = ZERO
            DO 10 J = 1, N
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
 10         CONTINUE
            IF( SMIN.LE.ZERO ) THEN
               INFO = -11
            ELSE IF( N.GT.0 ) THEN
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            ELSE
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
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYSVXX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
      // Compute row and column scalings to equilibrate the matrix A.
*
         CALL DSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
      // Equilibrate the matrix.
*
            CALL DLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )
            RCEQU = LSAME( EQUED, 'Y' )
         END IF
      END IF
*
      // Scale the right-hand side.
*
      IF( RCEQU ) CALL DLASCL2( N, NRHS, S, B, LDB )
*
      IF( NOFACT .OR. EQUIL ) THEN
*
         // Compute the LDL^T or UDU^T factorization of A.
*
         CALL DLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL DSYTRF( UPLO, N, AF, LDAF, IPIV, WORK, 5*MAX(1,N), INFO )
*
         // Return if INFO is non-zero.
*
         IF( INFO.GT.0 ) THEN
*
            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.
*
            IF ( N.GT.0 ) RPVGRW = DLA_SYRPVGRW(UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, WORK )
            RETURN
         END IF
      END IF
*
      // Compute the reciprocal pivot growth factor RPVGRW.
*
      IF ( N.GT.0 ) RPVGRW = DLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, WORK )
*
      // Compute the solution matrix X.
*
      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DSYTRS( UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
*
      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.
*
      CALL DSYRFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )
*
      // Scale solutions.
*
      IF ( RCEQU ) THEN
         CALL DLASCL2 ( N, NRHS, S, X, LDX )
      END IF
*
      RETURN
*
      // End of DSYSVXX
*
      END
