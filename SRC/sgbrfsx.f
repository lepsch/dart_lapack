      SUBROUTINE SGBRFSX( TRANS, EQUED, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS, EQUED;
      int                INFO, LDAB, LDAFB, LDB, LDX, N, KL, KU, NRHS, NPARAMS, N_ERR_BNDS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), X( LDX , * ),WORK( * )       REAL               R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * )
      // ..

*  ==================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      REAL               ITREF_DEFAULT, ITHRESH_DEFAULT, COMPONENTWISE_DEFAULT
      REAL               RTHRESH_DEFAULT, DZTHRESH_DEFAULT
      const              ITREF_DEFAULT = 1.0 ;
      const              ITHRESH_DEFAULT = 10.0 ;
      const              COMPONENTWISE_DEFAULT = 1.0 ;
      const              RTHRESH_DEFAULT = 0.5 ;
      const              DZTHRESH_DEFAULT = 0.25 ;
      int                LA_LINRX_ITREF_I, LA_LINRX_ITHRESH_I, LA_LINRX_CWISE_I;
      const              LA_LINRX_ITREF_I = 1, LA_LINRX_ITHRESH_I = 2 ;
      const              LA_LINRX_CWISE_I = 3 ;
      int                LA_LINRX_TRUST_I, LA_LINRX_ERR_I, LA_LINRX_RCOND_I;
      const              LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2 ;
      const              LA_LINRX_RCOND_I = 3 ;
      // ..
      // .. Local Scalars ..
      String   (1)       NORM;
      bool               ROWEQU, COLEQU, NOTRAN;
      int                J, TRANS_TYPE, PREC_TYPE, REF_TYPE;
      int                N_NORMS;
      REAL               ANORM, RCOND_TMP
      REAL               ILLRCOND_THRESH, ERR_LBND, CWISE_WRONG
      bool               IGNORE_CWISE;
      int                ITHRESH;
      REAL               RTHRESH, UNSTABLE_THRESH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, SGBCON
      // EXTERNAL SLA_GBRFSX_EXTENDED
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, ILATRANS, ILAPREC
      // EXTERNAL SLAMCH, SLANGB, SLA_GBRCOND
      REAL               SLAMCH, SLANGB, SLA_GBRCOND
      bool               LSAME;
      int                ILATRANS, ILAPREC;
      // ..
      // .. Executable Statements ..

      // Check the input parameters.

      INFO = 0
      TRANS_TYPE = ILATRANS( TRANS )
      REF_TYPE = INT( ITREF_DEFAULT )
      IF ( NPARAMS .GE. LA_LINRX_ITREF_I ) THEN
         IF ( PARAMS( LA_LINRX_ITREF_I ) .LT. 0.0 ) THEN
            PARAMS( LA_LINRX_ITREF_I ) = ITREF_DEFAULT
         } else {
            REF_TYPE = PARAMS( LA_LINRX_ITREF_I )
         END IF
      END IF

      // Set default parameters.

      ILLRCOND_THRESH = REAL( N ) * SLAMCH( 'Epsilon' )
      ITHRESH = INT( ITHRESH_DEFAULT )
      RTHRESH = RTHRESH_DEFAULT
      UNSTABLE_THRESH = DZTHRESH_DEFAULT
      IGNORE_CWISE = COMPONENTWISE_DEFAULT .EQ. 0.0

      IF ( NPARAMS.GE.LA_LINRX_ITHRESH_I ) THEN
         IF ( PARAMS( LA_LINRX_ITHRESH_I ).LT.0.0 ) THEN
            PARAMS( LA_LINRX_ITHRESH_I ) = ITHRESH
         } else {
            ITHRESH = INT( PARAMS( LA_LINRX_ITHRESH_I ) )
         END IF
      END IF
      IF ( NPARAMS.GE.LA_LINRX_CWISE_I ) THEN
         IF ( PARAMS( LA_LINRX_CWISE_I ).LT.0.0 ) THEN
            IF ( IGNORE_CWISE ) THEN
               PARAMS( LA_LINRX_CWISE_I ) = 0.0
            } else {
               PARAMS( LA_LINRX_CWISE_I ) = 1.0
            END IF
         } else {
            IGNORE_CWISE = PARAMS( LA_LINRX_CWISE_I ) .EQ. 0.0
         END IF
      END IF
      IF ( REF_TYPE .EQ. 0 .OR. N_ERR_BNDS .EQ. 0 ) THEN
         N_NORMS = 0
      ELSE IF ( IGNORE_CWISE ) THEN
         N_NORMS = 1
      } else {
         N_NORMS = 2
      END IF

      NOTRAN = LSAME( TRANS, 'N' )
      ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
      COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )

      // Test input parameters.

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
        CALL XERBLA( 'SGBRFSX', -INFO )
        RETURN
      END IF

      // Quick return if possible.

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         RCOND = 1.0
         DO J = 1, NRHS
            BERR( J ) = 0.0
            IF ( N_ERR_BNDS .GE. 1 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0
            END IF
            IF ( N_ERR_BNDS .GE. 2 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 0.0
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 0.0
            END IF
            IF ( N_ERR_BNDS .GE. 3 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = 1.0
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 1.0
            END IF
         END DO
         RETURN
      END IF

      // Default to failure.

      RCOND = 0.0
      DO J = 1, NRHS
         BERR( J ) = 1.0
         IF ( N_ERR_BNDS .GE. 1 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0
            ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0
         END IF
         IF ( N_ERR_BNDS .GE. 2 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0
            ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0
         END IF
         IF ( N_ERR_BNDS .GE. 3 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = 0.0
            ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 0.0
         END IF
      END DO

      // Compute the norm of A and the reciprocal of the condition
      // number of A.

      IF( NOTRAN ) THEN
         NORM = 'I'
      } else {
         NORM = '1'
      END IF
      ANORM = SLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
      CALL SGBCON( NORM, N, KL, KU, AFB, LDAFB, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

      // Perform refinement on each right-hand side

      IF ( REF_TYPE .NE. 0 .AND. INFO .EQ. 0 ) THEN

         PREC_TYPE = ILAPREC( 'D' )

         IF ( NOTRAN ) THEN
            CALL SLA_GBRFSX_EXTENDED( PREC_TYPE, TRANS_TYPE,  N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, COLEQU, C, B, LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, WORK( N+1 ), WORK( 1 ), WORK( 2*N+1 ), WORK( 1 ), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, IGNORE_CWISE, INFO )
         } else {
            CALL SLA_GBRFSX_EXTENDED( PREC_TYPE, TRANS_TYPE,  N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, ROWEQU, R, B, LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, WORK( N+1 ), WORK( 1 ), WORK( 2*N+1 ), WORK( 1 ), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, IGNORE_CWISE, INFO )
         END IF
      END IF

      ERR_LBND = MAX( 10.0, SQRT( REAL( N ) ) ) * SLAMCH( 'Epsilon' )
      IF ( N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 1 ) THEN

      // Compute scaled normwise condition number cond(A*C).

         IF ( COLEQU .AND. NOTRAN ) THEN
            RCOND_TMP = SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, -1, C, INFO, WORK, IWORK )
         ELSE IF ( ROWEQU .AND. .NOT. NOTRAN ) THEN
            RCOND_TMP = SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, -1, R, INFO, WORK, IWORK )
         } else {
            RCOND_TMP = SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, 0, R, INFO, WORK, IWORK )
         END IF
         DO J = 1, NRHS

      // Cap the error at 1.0.

            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I .AND. ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .GT. 1.0 ) ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0

      // Threshold the error (see LAWN).

            IF ( RCOND_TMP .LT. ILLRCOND_THRESH ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 0.0
               IF ( INFO .LE. N ) INFO = N + J
            ELSE IF ( ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .LT. ERR_LBND ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0
            END IF

      // Save the condition number.

            IF ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            END IF

         END DO
      END IF

      IF (N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 2) THEN

      // Compute componentwise condition number cond(A*diag(Y(:,J))) for
      // each right-hand side using the current solution as an estimate of
     t // he true solution.  If the componentwise error estimate is too
      // large, then the solution is a lousy estimate of truth and the
      // estimated RCOND may be too optimistic.  To avoid misleading users,
     t // he inverse condition number is set to 0.0 when the estimated
      // cwise error is at least CWISE_WRONG.

         CWISE_WRONG = SQRT( SLAMCH( 'Epsilon' ) )
         DO J = 1, NRHS
            IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .LT. CWISE_WRONG ) THEN                RCOND_TMP = SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, 1, X( 1, J ), INFO, WORK, IWORK )
            } else {
               RCOND_TMP = 0.0
            END IF

      // Cap the error at 1.0.

            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I .AND. ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .GT. 1.0 ) ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0

      // Threshold the error (see LAWN).

            IF ( RCOND_TMP .LT. ILLRCOND_THRESH ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 0.0
               IF ( PARAMS( LA_LINRX_CWISE_I ) .EQ. 1.0 .AND. INFO.LT.N + J ) INFO = N + J             ELSE IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .LT. ERR_LBND ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0
            END IF

      // Save the condition number.

            IF ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            END IF

         END DO
      END IF

      RETURN

      // End of SGBRFSX

      }
