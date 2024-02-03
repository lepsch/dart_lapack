      SUBROUTINE SGBSVXX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, TRANS;
      int                INFO, LDAB, LDAFB, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      REAL               RCOND, RPVGRW
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), X( LDX , * ),WORK( * )       REAL               R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * )
      // ..

*  ==================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I;
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I;
      int                CMP_ERR_I, PIV_GROWTH_I;
      const              FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 ;
      const              RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 ;
      const              CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 ;
      // ..
      // .. Local Scalars ..
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      int                INFEQU, I, J, KL, KU;
      REAL               AMAX, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, SMLNUM
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, SLAMCH, SLA_GBRPVGRW
      bool               LSAME;
      REAL               SLAMCH, SLA_GBRPVGRW
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGBEQUB, SGBTRF, SGBTRS, SLACPY, SLAQGB, XERBLA, SLASCL2, SGBRFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      SMLNUM = SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      if ( NOFACT .OR. EQUIL ) {
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      } else {
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
      }

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in SGBRFSX.

      RPVGRW = ZERO

      // Test the input parameters.  PARAMS is not tested until SGBRFSX.

      if ( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT. LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( KL.LT.0 ) {
         INFO = -4
      } else if ( KU.LT.0 ) {
         INFO = -5
      } else if ( NRHS.LT.0 ) {
         INFO = -6
      } else if ( LDAB.LT.KL+KU+1 ) {
         INFO = -8
      } else if ( LDAFB.LT.2*KL+KU+1 ) {
         INFO = -10
      } else if ( LSAME( FACT, 'F' ) .AND. .NOT. ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) {
         INFO = -12
      } else {
         if ( ROWEQU ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
            } // 10
            if ( RCMIN.LE.ZERO ) {
               INFO = -13
            } else if ( N.GT.0 ) {
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               ROWCND = ONE
            }
         }
         if ( COLEQU .AND. INFO.EQ.0 ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            for (J = 1; J <= N; J++) { // 20
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
            } // 20
            if ( RCMIN.LE.ZERO ) {
               INFO = -14
            } else if ( N.GT.0 ) {
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               COLCND = ONE
            }
         }
         if ( INFO.EQ.0 ) {
            if ( LDB.LT.MAX( 1, N ) ) {
               INFO = -15
            } else if ( LDX.LT.MAX( 1, N ) ) {
               INFO = -16
            }
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SGBSVXX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         sgbequb(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU.EQ.0 ) {

      // Equilibrate the matrix.

            slaqgb(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         }

      // If the scaling factors are not applied, set them to 1.0.

         if ( .NOT.ROWEQU ) {
            for (J = 1; J <= N; J++) {
               R( J ) = 1.0
            END DO
         }
         if ( .NOT.COLEQU ) {
            for (J = 1; J <= N; J++) {
               C( J ) = 1.0
            END DO
         }
      }

      // Scale the right hand side.

      if ( NOTRAN ) {
         IF( ROWEQU ) CALL SLASCL2(N, NRHS, R, B, LDB)
      } else {
         IF( COLEQU ) CALL SLASCL2(N, NRHS, C, B, LDB)
      }

      if ( NOFACT .OR. EQUIL ) {

         // Compute the LU factorization of A.

         DO 40, J = 1, N
            DO 30, I = KL+1, 2*KL+KU+1
               AFB( I, J ) = AB( I-KL, J )
            } // 30
         } // 40
         sgbtrf(N, N, KL, KU, AFB, LDAFB, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = SLA_GBRPVGRW( N, KL, KU, INFO, AB, LDAB, AFB, LDAFB )
            RETURN
         }
      }

      // Compute the reciprocal pivot growth factor RPVGRW.

      RPVGRW = SLA_GBRPVGRW( N, KL, KU, N, AB, LDAB, AFB, LDAFB )

      // Compute the solution matrix X.

      slacpy('Full', N, NRHS, B, LDB, X, LDX );
      sgbtrs(TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      sgbrfsx(TRANS, EQUED, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO );

      // Scale solutions.

      if ( COLEQU .AND. NOTRAN ) {
         slascl2(N, NRHS, C, X, LDX );
      } else if ( ROWEQU .AND. .NOT.NOTRAN ) {
         slascl2(N, NRHS, R, X, LDX );
      }

      RETURN

      // End of SGBSVXX

      }
