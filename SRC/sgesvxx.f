      SUBROUTINE SGESVXX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, TRANS;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      REAL               RCOND, RPVGRW
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX , * ),WORK( * )       REAL               R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * )
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
      int                INFEQU, J;
      REAL               AMAX, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, SMLNUM
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, SLAMCH, SLA_GERPVGRW
      bool               LSAME;
      REAL               SLAMCH, SLA_GERPVGRW
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEEQUB, SGETRF, SGETRS, SLACPY, SLAQGE, XERBLA, SLASCL2, SGERFSX
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
         ROWEQU = false;
         COLEQU = false;
      } else {
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
      }

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in SGERFSX.

      RPVGRW = ZERO

      // Test the input parameters.  PARAMS is not tested until SGERFSX.

      if ( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT. LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LSAME( FACT, 'F' ) .AND. .NOT. ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) {
         INFO = -10
      } else {
         if ( ROWEQU ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
            } // 10
            if ( RCMIN.LE.ZERO ) {
               INFO = -11
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
               INFO = -12
            } else if ( N.GT.0 ) {
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               COLCND = ONE
            }
         }
         if ( INFO.EQ.0 ) {
            if ( LDB.LT.MAX( 1, N ) ) {
               INFO = -14
            } else if ( LDX.LT.MAX( 1, N ) ) {
               INFO = -16
            }
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SGESVXX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         sgeequb(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU.EQ.0 ) {

      // Equilibrate the matrix.

            slaqge(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         }

      // If the scaling factors are not applied, set them to 1.0.

         if ( .NOT.ROWEQU ) {
            for (J = 1; J <= N; J++) {
               R( J ) = 1.0
            }
         }
         if ( .NOT.COLEQU ) {
            for (J = 1; J <= N; J++) {
               C( J ) = 1.0
            }
         }
      }

      // Scale the right-hand side.

      if ( NOTRAN ) {
         if (ROWEQU) CALL SLASCL2( N, NRHS, R, B, LDB );
      } else {
         if (COLEQU) CALL SLASCL2( N, NRHS, C, B, LDB );
      }

      if ( NOFACT .OR. EQUIL ) {

         // Compute the LU factorization of A.

         slacpy('Full', N, N, A, LDA, AF, LDAF );
         sgetrf(N, N, AF, LDAF, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = SLA_GERPVGRW( N, INFO, A, LDA, AF, LDAF )
            RETURN
         }
      }

      // Compute the reciprocal pivot growth factor RPVGRW.

      RPVGRW = SLA_GERPVGRW( N, N, A, LDA, AF, LDAF )

      // Compute the solution matrix X.

      slacpy('Full', N, NRHS, B, LDB, X, LDX );
      sgetrs(TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      sgerfsx(TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO );

      // Scale solutions.

      if ( COLEQU .AND. NOTRAN ) {
         slascl2(N, NRHS, C, X, LDX );
      } else if ( ROWEQU .AND. .NOT.NOTRAN ) {
         slascl2(N, NRHS, R, X, LDX );
      }

      RETURN

      // End of SGESVXX

      }
