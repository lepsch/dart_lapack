      SUBROUTINE ZGBRFSX( TRANS, EQUED, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS, EQUED;
      int                INFO, LDAB, LDAFB, LDB, LDX, N, KL, KU, NRHS, NPARAMS, N_ERR_BNDS
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), X( LDX , * ),WORK( * )       DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * ), RWORK( * )
*     ..
*
*  ==================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   ITREF_DEFAULT, ITHRESH_DEFAULT
      DOUBLE PRECISION   COMPONENTWISE_DEFAULT, RTHRESH_DEFAULT
      DOUBLE PRECISION   DZTHRESH_DEFAULT
      PARAMETER          ( ITREF_DEFAULT = 1.0D+0 )
      PARAMETER          ( ITHRESH_DEFAULT = 10.0D+0 )
      PARAMETER          ( COMPONENTWISE_DEFAULT = 1.0D+0 )
      PARAMETER          ( RTHRESH_DEFAULT = 0.5D+0 )
      PARAMETER          ( DZTHRESH_DEFAULT = 0.25D+0 )
      int                LA_LINRX_ITREF_I, LA_LINRX_ITHRESH_I, LA_LINRX_CWISE_I       PARAMETER          ( LA_LINRX_ITREF_I = 1, LA_LINRX_ITHRESH_I = 2 )
      PARAMETER          ( LA_LINRX_CWISE_I = 3 )
      int                LA_LINRX_TRUST_I, LA_LINRX_ERR_I, LA_LINRX_RCOND_I
      PARAMETER          ( LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2 )
      PARAMETER          ( LA_LINRX_RCOND_I = 3 )
*     ..
*     .. Local Scalars ..
      String   (1)       NORM;
      LOGICAL            ROWEQU, COLEQU, NOTRAN, IGNORE_CWISE
      int                J, TRANS_TYPE, PREC_TYPE, REF_TYPE, N_NORMS, ITHRESH       DOUBLE PRECISION   ANORM, RCOND_TMP, ILLRCOND_THRESH, ERR_LBND, CWISE_WRONG, RTHRESH, UNSTABLE_THRESH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGBCON, ZLA_GBRFSX_EXTENDED
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT, TRANSFER
*     ..
*     .. External Functions ..
      EXTERNAL           LSAME, ILAPREC
      EXTERNAL           DLAMCH, ZLANGB, ZLA_GBRCOND_X, ZLA_GBRCOND_C
      DOUBLE PRECISION   DLAMCH, ZLANGB, ZLA_GBRCOND_X, ZLA_GBRCOND_C
      LOGICAL            LSAME
      int                ILATRANS, ILAPREC
*     ..
*     .. Executable Statements ..
*
*     Check the input parameters.
*
      INFO = 0
      TRANS_TYPE = ILATRANS( TRANS )
      REF_TYPE = INT( ITREF_DEFAULT )
      IF ( NPARAMS .GE. LA_LINRX_ITREF_I ) THEN
         IF ( PARAMS( LA_LINRX_ITREF_I ) .LT. 0.0D+0 ) THEN
            PARAMS( LA_LINRX_ITREF_I ) = ITREF_DEFAULT
         ELSE
            REF_TYPE = PARAMS( LA_LINRX_ITREF_I )
         END IF
      END IF
*
*     Set default parameters.
*
      ILLRCOND_THRESH = DBLE( N ) * DLAMCH( 'Epsilon' )
      ITHRESH = INT( ITHRESH_DEFAULT )
      RTHRESH = RTHRESH_DEFAULT
      UNSTABLE_THRESH = DZTHRESH_DEFAULT
      IGNORE_CWISE = COMPONENTWISE_DEFAULT .EQ. 0.0D+0
*
      IF ( NPARAMS.GE.LA_LINRX_ITHRESH_I ) THEN
         IF ( PARAMS( LA_LINRX_ITHRESH_I ).LT.0.0D+0 ) THEN
            PARAMS( LA_LINRX_ITHRESH_I ) = ITHRESH
         ELSE
            ITHRESH = INT( PARAMS( LA_LINRX_ITHRESH_I ) )
         END IF
      END IF
      IF ( NPARAMS.GE.LA_LINRX_CWISE_I ) THEN
         IF ( PARAMS( LA_LINRX_CWISE_I ).LT.0.0D+0 ) THEN
            IF ( IGNORE_CWISE ) THEN
               PARAMS( LA_LINRX_CWISE_I ) = 0.0D+0
            ELSE
               PARAMS( LA_LINRX_CWISE_I ) = 1.0D+0
            END IF
         ELSE
            IGNORE_CWISE = PARAMS( LA_LINRX_CWISE_I ) .EQ. 0.0D+0
         END IF
      END IF
      IF ( REF_TYPE .EQ. 0 .OR. N_ERR_BNDS .EQ. 0 ) THEN
         N_NORMS = 0
      ELSE IF ( IGNORE_CWISE ) THEN
         N_NORMS = 1
      ELSE
         N_NORMS = 2
      END IF
*
      NOTRAN = LSAME( TRANS, 'N' )
      ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
      COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
*
*     Test input parameters.
*
      IF( TRANS_TYPE.EQ.-1 ) THEN
        INFO = -1
      ELSE IF( .NOT.ROWEQU .AND. .NOT.COLEQU .AND. .NOT.LSAME( EQUED, 'N' ) ) THEN
        INFO = -2
      ELSE IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( KL.LT.0 ) THEN
        INFO = -4
      ELSE IF( KU.LT.0 ) THEN
        INFO = -5
      ELSE IF( NRHS.LT.0 ) THEN
        INFO = -6
      ELSE IF( LDAB.LT.KL+KU+1 ) THEN
        INFO = -8
      ELSE IF( LDAFB.LT.2*KL+KU+1 ) THEN
        INFO = -10
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
        INFO = -13
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
        INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZGBRFSX', -INFO )
        RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         RCOND = 1.0D+0
         DO J = 1, NRHS
            BERR( J ) = 0.0D+0
            IF ( N_ERR_BNDS .GE. 1 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
            END IF
            IF ( N_ERR_BNDS .GE. 2 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 0.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 0.0D+0
            END IF
            IF ( N_ERR_BNDS .GE. 3 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 1.0D+0
            END IF
         END DO
         RETURN
      END IF
*
*     Default to failure.
*
      RCOND = 0.0D+0
      DO J = 1, NRHS
         BERR( J ) = 1.0D+0
         IF ( N_ERR_BNDS .GE. 1 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
         END IF
         IF ( N_ERR_BNDS .GE. 2 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
         END IF
         IF ( N_ERR_BNDS .GE. 3 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = 0.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 0.0D+0
         END IF
      END DO
*
*     Compute the norm of A and the reciprocal of the condition
*     number of A.
*
      IF( NOTRAN ) THEN
         NORM = 'I'
      ELSE
         NORM = '1'
      END IF
      ANORM = ZLANGB( NORM, N, KL, KU, AB, LDAB, RWORK )
      CALL ZGBCON( NORM, N, KL, KU, AFB, LDAFB, IPIV, ANORM, RCOND, WORK, RWORK, INFO )
*
*     Perform refinement on each right-hand side
*
      IF ( REF_TYPE .NE. 0 .AND. INFO .EQ. 0 ) THEN

         PREC_TYPE = ILAPREC( 'E' )

         IF ( NOTRAN ) THEN
            CALL ZLA_GBRFSX_EXTENDED( PREC_TYPE, TRANS_TYPE,  N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, COLEQU, C, B, LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, WORK, RWORK, WORK(N+1), TRANSFER (RWORK(1:2*N), (/ (ZERO, ZERO) /), N), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, IGNORE_CWISE, INFO )
         ELSE
            CALL ZLA_GBRFSX_EXTENDED( PREC_TYPE, TRANS_TYPE,  N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, ROWEQU, R, B, LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, WORK, RWORK, WORK(N+1), TRANSFER (RWORK(1:2*N), (/ (ZERO, ZERO) /), N), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, IGNORE_CWISE, INFO )
         END IF
      END IF

      ERR_LBND = MAX( 10.0D+0, SQRT( DBLE( N ) ) ) * DLAMCH( 'Epsilon' )
      IF (N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 1) THEN
*
*     Compute scaled normwise condition number cond(A*C).
*
         IF ( COLEQU .AND. NOTRAN ) THEN
            RCOND_TMP = ZLA_GBRCOND_C( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, C, .TRUE., INFO, WORK, RWORK )
         ELSE IF ( ROWEQU .AND. .NOT. NOTRAN ) THEN
            RCOND_TMP = ZLA_GBRCOND_C( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, R, .TRUE., INFO, WORK, RWORK )
         ELSE
            RCOND_TMP = ZLA_GBRCOND_C( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, C, .FALSE., INFO, WORK, RWORK )
         END IF
         DO J = 1, NRHS
*
*     Cap the error at 1.0.
*
            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I .AND. ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .GT. 1.0D+0) ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
*
*     Threshold the error (see LAWN).
*
            IF ( RCOND_TMP .LT. ILLRCOND_THRESH ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 0.0D+0
               IF ( INFO .LE. N ) INFO = N + J
            ELSE IF ( ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .LT. ERR_LBND ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
            END IF
*
*     Save the condition number.
*
            IF ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            END IF

         END DO
      END IF

      IF (N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 2) THEN
*
*     Compute componentwise condition number cond(A*diag(Y(:,J))) for
*     each right-hand side using the current solution as an estimate of
*     the true solution.  If the componentwise error estimate is too
*     large, then the solution is a lousy estimate of truth and the
*     estimated RCOND may be too optimistic.  To avoid misleading users,
*     the inverse condition number is set to 0.0 when the estimated
*     cwise error is at least CWISE_WRONG.
*
         CWISE_WRONG = SQRT( DLAMCH( 'Epsilon' ) )
         DO J = 1, NRHS
            IF (ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .LT. CWISE_WRONG ) THEN                RCOND_TMP = ZLA_GBRCOND_X( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, X( 1, J ), INFO, WORK, RWORK )
            ELSE
               RCOND_TMP = 0.0D+0
            END IF
*
*     Cap the error at 1.0.
*
            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I .AND. ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .GT. 1.0D+0 ) ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
*
*     Threshold the error (see LAWN).
*
            IF ( RCOND_TMP .LT. ILLRCOND_THRESH ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 0.0D+0
               IF ( PARAMS( LA_LINRX_CWISE_I ) .EQ. 1.0D+0 .AND. INFO.LT.N + J ) INFO = N + J             ELSE IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .LT. ERR_LBND ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
            END IF
*
*     Save the condition number.
*
            IF ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            END IF

         END DO
      END IF
*
      RETURN
*
*     End of ZGBRFSX
*
      END
