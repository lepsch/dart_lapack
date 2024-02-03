      SUBROUTINE DPOSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double             RCOND, RPVGRW;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX, * ), WORK( * );;
      double             S( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * );
      // ..

*  ==================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I;
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I;
      int                CMP_ERR_I, PIV_GROWTH_I;
      const              FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 ;
      const              RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 ;
      const              CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU;
      int                INFEQU, J;
      double             AMAX, BIGNUM, SMIN, SMAX, SCOND, SMLNUM;
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, DLAMCH, DLA_PORPVGRW
      bool               LSAME;
      double             DLAMCH, DLA_PORPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DPOEQUB, DPOTRF, DPOTRS, DLACPY, DLAQSY, XERBLA, DLASCL2, DPORFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      if ( NOFACT .OR. EQUIL ) {
         EQUED = 'N'
         RCEQU = .FALSE.
      } else {
         RCEQU = LSAME( EQUED, 'Y' )
      ENDIF

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in DPORFSX.

      RPVGRW = ZERO

      // Test the input parameters.  PARAMS is not tested until DPORFSX.

      if ( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT. LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LSAME( FACT, 'F' ) .AND. .NOT. ( RCEQU .OR. LSAME( EQUED, 'N' ) ) ) {
         INFO = -9
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM
            SMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
            } // 10
            if ( SMIN.LE.ZERO ) {
               INFO = -10
            } else if ( N.GT.0 ) {
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            } else {
               SCOND = ONE
            }
         }
         if ( INFO.EQ.0 ) {
            if ( LDB.LT.MAX( 1, N ) ) {
               INFO = -12
            } else if ( LDX.LT.MAX( 1, N ) ) {
               INFO = -14
            }
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('DPOSVXX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         dpoequb(N, A, LDA, S, SCOND, AMAX, INFEQU );
         if ( INFEQU.EQ.0 ) {

      // Equilibrate the matrix.

            dlaqsy(UPLO, N, A, LDA, S, SCOND, AMAX, EQUED );
            RCEQU = LSAME( EQUED, 'Y' )
         }
      }

      // Scale the right-hand side.

      IF( RCEQU ) CALL DLASCL2( N, NRHS, S, B, LDB )

      if ( NOFACT .OR. EQUIL ) {

         // Compute the Cholesky factorization of A.

         dlacpy(UPLO, N, N, A, LDA, AF, LDAF );
         dpotrf(UPLO, N, AF, LDAF, INFO );

         // Return if INFO is non-zero.

         if ( INFO.NE.0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = DLA_PORPVGRW( UPLO, INFO, A, LDA, AF, LDAF, WORK )
            RETURN
         ENDIF
      }

      // Compute the reciprocal growth factor RPVGRW.

      RPVGRW = DLA_PORPVGRW( UPLO, N, A, LDA, AF, LDAF, WORK )

      // Compute the solution matrix X.

      dlacpy('Full', N, NRHS, B, LDB, X, LDX );
      dpotrs(UPLO, N, NRHS, AF, LDAF, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      dporfsx(UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO );


      // Scale solutions.

      if ( RCEQU ) {
         dlascl2(N, NRHS, S, X, LDX );
      }

      RETURN

      // End of DPOSVXX

      }
