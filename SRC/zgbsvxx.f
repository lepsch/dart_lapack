      SUBROUTINE ZGBSVXX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, TRANS;
      int                INFO, LDAB, LDAFB, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double             RCOND, RPVGRW;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), X( LDX , * ),WORK( * );
      double             R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * ), RWORK( * );
      // ..

*  ==================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
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
      double             AMAX, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, SMLNUM;
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, DLAMCH, ZLA_GBRPVGRW
      bool               LSAME;
      double             DLAMCH, ZLA_GBRPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGBEQUB, ZGBTRF, ZGBTRS, ZLACPY, ZLAQGB, XERBLA, ZLASCL2, ZGBRFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      if ( NOFACT || EQUIL ) {
         EQUED = 'N'
         ROWEQU = false;
         COLEQU = false;
      } else {
         ROWEQU = LSAME( EQUED, 'R' ) || LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) || LSAME( EQUED, 'B' )
      }

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in ZGBRFSX.

      RPVGRW = ZERO

      // Test the input parameters.  PARAMS is not tested until ZGERFSX.

      if ( !NOFACT && !EQUIL && !LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( !NOTRAN && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( KL < 0 ) {
         INFO = -4
      } else if ( KU < 0 ) {
         INFO = -5
      } else if ( NRHS < 0 ) {
         INFO = -6
      } else if ( LDAB < KL+KU+1 ) {
         INFO = -8
      } else if ( LDAFB < 2*KL+KU+1 ) {
         INFO = -10
      } else if ( LSAME( FACT, 'F' ) && !( ROWEQU || COLEQU || LSAME( EQUED, 'N' ) ) ) {
         INFO = -12
      } else {
         if ( ROWEQU ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
            } // 10
            if ( RCMIN <= ZERO ) {
               INFO = -13
            } else if ( N > 0 ) {
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               ROWCND = ONE
            }
         }
         if ( COLEQU && INFO == 0 ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            for (J = 1; J <= N; J++) { // 20
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
            } // 20
            if ( RCMIN <= ZERO ) {
               INFO = -14
            } else if ( N > 0 ) {
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               COLCND = ONE
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < MAX( 1, N ) ) {
               INFO = -15
            } else if ( LDX < MAX( 1, N ) ) {
               INFO = -16
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZGBSVXX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         zgbequb(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

      // Equilibrate the matrix.

            zlaqgb(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = LSAME( EQUED, 'R' ) || LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) || LSAME( EQUED, 'B' )
         }

      // If the scaling factors are not applied, set them to 1.0.

         if ( !ROWEQU ) {
            for (J = 1; J <= N; J++) {
               R( J ) = 1.0;
            }
         }
         if ( !COLEQU ) {
            for (J = 1; J <= N; J++) {
               C( J ) = 1.0;
            }
         }
      }

      // Scale the right-hand side.

      if ( NOTRAN ) {
         if (ROWEQU) CALL ZLASCL2( N, NRHS, R, B, LDB );
      } else {
         if (COLEQU) CALL ZLASCL2( N, NRHS, C, B, LDB );
      }

      if ( NOFACT || EQUIL ) {

         // Compute the LU factorization of A.

         DO 40, J = 1, N
            DO 30, I = KL+1, 2*KL+KU+1
               AFB( I, J ) = AB( I-KL, J )
            } // 30
         } // 40
         zgbtrf(N, N, KL, KU, AFB, LDAFB, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = ZLA_GBRPVGRW( N, KL, KU, INFO, AB, LDAB, AFB, LDAFB )
            RETURN
         }
      }

      // Compute the reciprocal pivot growth factor RPVGRW.

      RPVGRW = ZLA_GBRPVGRW( N, KL, KU, N, AB, LDAB, AFB, LDAFB )

      // Compute the solution matrix X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zgbtrs(TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      zgbrfsx(TRANS, EQUED, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO );


      // Scale solutions.

      if ( COLEQU && NOTRAN ) {
         zlascl2(N, NRHS, C, X, LDX );
      } else if ( ROWEQU && !NOTRAN ) {
         zlascl2(N, NRHS, R, X, LDX );
      }

      RETURN

      // End of ZGBSVXX

      }
