      SUBROUTINE SPOSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS
      REAL               RCOND, RPVGRW
*     ..
*     .. Array Arguments ..
      int                IWORK( * )
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX, * ), WORK( * )       REAL               S( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * )
*     ..
*
*  ==================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I
      int                CMP_ERR_I, PIV_GROWTH_I
      PARAMETER          ( FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 )
      PARAMETER          ( RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 )
      PARAMETER          ( CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 )
*     ..
*     .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU;
      int                INFEQU, J
      REAL               AMAX, BIGNUM, SMIN, SMAX, SCOND, SMLNUM
*     ..
*     .. External Functions ..
      EXTERNAL           LSAME, SLAMCH, SLA_PORPVGRW
      bool               LSAME;
      REAL               SLAMCH, SLA_PORPVGRW
*     ..
*     .. External Subroutines ..
      EXTERNAL           SPOEQUB, SPOTRF, SPOTRS, SLACPY, SLAQSY, XERBLA, SLASCL2, SPORFSX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      SMLNUM = SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         RCEQU = .FALSE.
      ELSE
         RCEQU = LSAME( EQUED, 'Y' )
      ENDIF
*
*     Default is failure.  If an input parameter is wrong or
*     factorization fails, make everything look horrible.  Only the
*     pivot growth is set here, the rest is initialized in SPORFSX.
*
      RPVGRW = ZERO
*
*     Test the input parameters.  PARAMS is not tested until SPORFSX.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT. LSAME( FACT, 'F' ) ) THEN
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
         IF ( RCEQU ) THEN
            SMIN = BIGNUM
            SMAX = ZERO
            DO 10 J = 1, N
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
 10         CONTINUE
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
         CALL XERBLA( 'SPOSVXX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*     Compute row and column scalings to equilibrate the matrix A.
*
         CALL SPOEQUB( N, A, LDA, S, SCOND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*     Equilibrate the matrix.
*
            CALL SLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )
            RCEQU = LSAME( EQUED, 'Y' )
         END IF
      END IF
*
*     Scale the right-hand side.
*
      IF( RCEQU ) CALL SLASCL2( N, NRHS, S, B, LDB )
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the Cholesky factorization of A.
*
         CALL SLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL SPOTRF( UPLO, N, AF, LDAF, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.NE.0 ) THEN
*
*           Pivot in column INFO is exactly 0
*           Compute the reciprocal pivot growth factor of the
*           leading rank-deficient INFO columns of A.
*
            RPVGRW = SLA_PORPVGRW( UPLO, INFO, A, LDA, AF, LDAF, WORK )
            RETURN
         ENDIF
      END IF
*
*     Compute the reciprocal growth factor RPVGRW.
*
      RPVGRW = SLA_PORPVGRW( UPLO, N, A, LDA, AF, LDAF, WORK )
*
*     Compute the solution matrix X.
*
      CALL SLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL SPOTRS( UPLO, N, NRHS, AF, LDAF, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL SPORFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )

*
*     Scale solutions.
*
      IF ( RCEQU ) THEN
         CALL SLASCL2 ( N, NRHS, S, X, LDX )
      END IF
*
      RETURN
*
*     End of SPOSVXX
*
      END
