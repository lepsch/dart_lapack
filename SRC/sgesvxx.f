      SUBROUTINE SGESVXX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             EQUED, FACT, TRANS;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      REAL               RCOND, RPVGRW
*     ..
*     .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX , * ),WORK( * )       REAL               R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * )
*     ..
*
*  ==================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I;
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I;
      int                CMP_ERR_I, PIV_GROWTH_I;
      PARAMETER          ( FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 )
      PARAMETER          ( RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 )
      PARAMETER          ( CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 )
*     ..
*     .. Local Scalars ..
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      int                INFEQU, J;
      REAL               AMAX, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, SMLNUM
*     ..
*     .. External Functions ..
      EXTERNAL           LSAME, SLAMCH, SLA_GERPVGRW
      bool               LSAME;
      REAL               SLAMCH, SLA_GERPVGRW
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEEQUB, SGETRF, SGETRS, SLACPY, SLAQGE, XERBLA, SLASCL2, SGERFSX
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      SMLNUM = SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      ELSE
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
      END IF
*
*     Default is failure.  If an input parameter is wrong or
*     factorization fails, make everything look horrible.  Only the
*     pivot growth is set here, the rest is initialized in SGERFSX.
*
      RPVGRW = ZERO
*
*     Test the input parameters.  PARAMS is not tested until SGERFSX.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT. LSAME( FACT, 'F' ) ) THEN
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
      ELSE
         IF( ROWEQU ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 10 J = 1, N
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
 10         CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -11
            ELSE IF( N.GT.0 ) THEN
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               ROWCND = ONE
            END IF
         END IF
         IF( COLEQU .AND. INFO.EQ.0 ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 20 J = 1, N
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
 20         CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -12
            ELSE IF( N.GT.0 ) THEN
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
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
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGESVXX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*     Compute row and column scalings to equilibrate the matrix A.
*
         CALL SGEEQUB( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*     Equilibrate the matrix.
*
            CALL SLAQGE( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED )
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         END IF
*
*     If the scaling factors are not applied, set them to 1.0.
*
         IF ( .NOT.ROWEQU ) THEN
            DO J = 1, N
               R( J ) = 1.0
            END DO
         END IF
         IF ( .NOT.COLEQU ) THEN
            DO J = 1, N
               C( J ) = 1.0
            END DO
         END IF
      END IF
*
*     Scale the right-hand side.
*
      IF( NOTRAN ) THEN
         IF( ROWEQU ) CALL SLASCL2( N, NRHS, R, B, LDB )
      ELSE
         IF( COLEQU ) CALL SLASCL2( N, NRHS, C, B, LDB )
      END IF
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the LU factorization of A.
*
         CALL SLACPY( 'Full', N, N, A, LDA, AF, LDAF )
         CALL SGETRF( N, N, AF, LDAF, IPIV, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.GT.0 ) THEN
*
*           Pivot in column INFO is exactly 0
*           Compute the reciprocal pivot growth factor of the
*           leading rank-deficient INFO columns of A.
*
            RPVGRW = SLA_GERPVGRW( N, INFO, A, LDA, AF, LDAF )
            RETURN
         END IF
      END IF
*
*     Compute the reciprocal pivot growth factor RPVGRW.
*
      RPVGRW = SLA_GERPVGRW( N, N, A, LDA, AF, LDAF )
*
*     Compute the solution matrix X.
*
      CALL SLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL SGETRS( TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL SGERFSX( TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )
*
*     Scale solutions.
*
      IF ( COLEQU .AND. NOTRAN ) THEN
         CALL SLASCL2 ( N, NRHS, C, X, LDX )
      ELSE IF ( ROWEQU .AND. .NOT.NOTRAN ) THEN
         CALL SLASCL2 ( N, NRHS, R, X, LDX )
      END IF
*
      RETURN
*
*     End of SGESVXX

      END
