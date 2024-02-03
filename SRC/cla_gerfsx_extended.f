      SUBROUTINE CLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A, LDA, AF, LDAF, IPIV, COLEQU, C, B, LDB, Y, LDY, BERR_OUT, N_NORMS, ERRS_N, ERRS_C, RES, AYB, DY, Y_TAIL, RCOND, ITHRESH, RTHRESH, DZ_UB, IGNORE_CWISE, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE, TRANS_TYPE, N_NORMS;
      bool               COLEQU, IGNORE_CWISE;
      int                ITHRESH;
      REAL               RTHRESH, DZ_UB
      // ..
      // .. Array Arguments
      int                IPIV( * );
      COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ), ERRS_N( NRHS, * ), ERRS_C( NRHS, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      String             TRANS;
      int                CNT, I, J,  X_STATE, Z_STATE, Y_PREC_STATE;
      REAL               YK, DYK, YMIN, NORMY, NORMX, NORMDX, DXRAT, DZRAT, PREVNORMDX, PREV_DZ_Z, DXRATMAX, DZRATMAX, DX_X, DZ_Z, FINAL_DX_X, FINAL_DZ_Z, EPS, HUGEVAL, INCR_THRESH
      bool               INCR_PREC;
      COMPLEX            ZDUM
      // ..
      // .. Parameters ..
      int                UNSTABLE_STATE, WORKING_STATE, CONV_STATE, NOPROG_STATE, BASE_RESIDUAL, EXTRA_RESIDUAL, EXTRA_Y;
      const              UNSTABLE_STATE = 0, WORKING_STATE = 1, CONV_STATE = 2, NOPROG_STATE = 3 ;
      const              BASE_RESIDUAL = 0, EXTRA_RESIDUAL = 1, EXTRA_Y = 2 ;
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I;
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I;
      int                CMP_ERR_I, PIV_GROWTH_I;
      const              FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 ;
      const              RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 ;
      const              CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 ;
      int                LA_LINRX_ITREF_I, LA_LINRX_ITHRESH_I, LA_LINRX_CWISE_I;
      const              LA_LINRX_ITREF_I = 1, LA_LINRX_ITHRESH_I = 2 ;
      const              LA_LINRX_CWISE_I = 3 ;
      int                LA_LINRX_TRUST_I, LA_LINRX_ERR_I, LA_LINRX_RCOND_I;
      const              LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2 ;
      const              LA_LINRX_RCOND_I = 3 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGETRS, CGEMV, BLAS_CGEMV_X, BLAS_CGEMV2_X, CLA_GEAMV, CLA_WWADDW, SLAMCH, CHLA_TRANSTYPE, CLA_LIN_BERR
      REAL               SLAMCH
      String             CHLA_TRANSTYPE;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      IF ( INFO.NE.0 ) RETURN
      TRANS = CHLA_TRANSTYPE(TRANS_TYPE)
      EPS = SLAMCH( 'Epsilon' )
      HUGEVAL = SLAMCH( 'Overflow' )
      // Force HUGEVAL to Inf
      HUGEVAL = HUGEVAL * HUGEVAL
      // Using HUGEVAL may lead to spurious underflows.
      INCR_THRESH = REAL( N ) * EPS

      for (J = 1; J <= NRHS; J++) {
         Y_PREC_STATE = EXTRA_RESIDUAL
         if ( Y_PREC_STATE .EQ. EXTRA_Y ) {
            for (I = 1; I <= N; I++) {
               Y_TAIL( I ) = 0.0
            END DO
         }

         DXRAT = 0.0
         DXRATMAX = 0.0
         DZRAT = 0.0
         DZRATMAX = 0.0
         FINAL_DX_X = HUGEVAL
         FINAL_DZ_Z = HUGEVAL
         PREVNORMDX = HUGEVAL
         PREV_DZ_Z = HUGEVAL
         DZ_Z = HUGEVAL
         DX_X = HUGEVAL

         X_STATE = WORKING_STATE
         Z_STATE = UNSTABLE_STATE
         INCR_PREC = .FALSE.

         for (CNT = 1; CNT <= ITHRESH; CNT++) {

          // Compute residual RES = B_s - op(A_s) * Y,
              // op(A) = A, A**T, or A**H depending on TRANS (and type).

            ccopy(N, B( 1, J ), 1, RES, 1 );
            if ( Y_PREC_STATE .EQ. BASE_RESIDUAL ) {
               cgemv(TRANS, N, N, (-1.0E+0,0.0E+0), A, LDA, Y( 1, J ), 1, (1.0E+0,0.0E+0), RES, 1);
            } else if (Y_PREC_STATE .EQ. EXTRA_RESIDUAL) {
               blas_cgemv_x(TRANS_TYPE, N, N, (-1.0E+0,0.0E+0), A, LDA, Y( 1, J ), 1, (1.0E+0,0.0E+0), RES, 1, PREC_TYPE );
            } else {
               blas_cgemv2_x(TRANS_TYPE, N, N, (-1.0E+0,0.0E+0), A, LDA, Y(1, J), Y_TAIL, 1, (1.0E+0,0.0E+0), RES, 1, PREC_TYPE);
            }

          // XXX: RES is no longer needed.
            ccopy(N, RES, 1, DY, 1 );
            cgetrs(TRANS, N, 1, AF, LDAF, IPIV, DY, N, INFO );

          // Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.

            NORMX = 0.0E+0
            NORMY = 0.0E+0
            NORMDX = 0.0E+0
            DZ_Z = 0.0E+0
            YMIN = HUGEVAL

            for (I = 1; I <= N; I++) {
               YK = CABS1( Y( I, J ) )
               DYK = CABS1( DY( I ) )

               if ( YK .NE. 0.0E+0 ) {
                  DZ_Z = MAX( DZ_Z, DYK / YK )
               } else if ( DYK .NE. 0.0 ) {
                  DZ_Z = HUGEVAL
               }

               YMIN = MIN( YMIN, YK )

               NORMY = MAX( NORMY, YK )

               if ( COLEQU ) {
                  NORMX = MAX( NORMX, YK * C( I ) )
                  NORMDX = MAX( NORMDX, DYK * C( I ) )
               } else {
                  NORMX = NORMY
                  NORMDX = MAX(NORMDX, DYK)
               }
            END DO

            if ( NORMX .NE. 0.0 ) {
               DX_X = NORMDX / NORMX
            } else if ( NORMDX .EQ. 0.0 ) {
               DX_X = 0.0
            } else {
               DX_X = HUGEVAL
            }

            DXRAT = NORMDX / PREVNORMDX
            DZRAT = DZ_Z / PREV_DZ_Z

          // Check termination criteria

            IF (.NOT.IGNORE_CWISE .AND. YMIN*RCOND .LT. INCR_THRESH*NORMY .AND. Y_PREC_STATE .LT. EXTRA_Y ) INCR_PREC = .TRUE.
             IF ( X_STATE .EQ. NOPROG_STATE .AND. DXRAT .LE. RTHRESH ) X_STATE = WORKING_STATE
            if ( X_STATE .EQ. WORKING_STATE ) {
               if (DX_X .LE. EPS) {
                  X_STATE = CONV_STATE
               } else if ( DXRAT .GT. RTHRESH ) {
                  if ( Y_PREC_STATE .NE. EXTRA_Y ) {
                     INCR_PREC = .TRUE.
                  } else {
                     X_STATE = NOPROG_STATE
                  }
               } else {
                  IF ( DXRAT .GT. DXRATMAX ) DXRATMAX = DXRAT
               }
               IF ( X_STATE .GT. WORKING_STATE ) FINAL_DX_X = DX_X
            }
             IF ( Z_STATE .EQ. UNSTABLE_STATE .AND. DZ_Z .LE. DZ_UB ) Z_STATE = WORKING_STATE             IF ( Z_STATE .EQ. NOPROG_STATE .AND. DZRAT .LE. RTHRESH ) Z_STATE = WORKING_STATE
            if ( Z_STATE .EQ. WORKING_STATE ) {
               if ( DZ_Z .LE. EPS ) {
                  Z_STATE = CONV_STATE
               } else if ( DZ_Z .GT. DZ_UB ) {
                  Z_STATE = UNSTABLE_STATE
                  DZRATMAX = 0.0
                  FINAL_DZ_Z = HUGEVAL
               } else if ( DZRAT .GT. RTHRESH ) {
                  if ( Y_PREC_STATE .NE. EXTRA_Y ) {
                     INCR_PREC = .TRUE.
                  } else {
                     Z_STATE = NOPROG_STATE
                  }
               } else {
                  IF ( DZRAT .GT. DZRATMAX ) DZRATMAX = DZRAT
               }
               IF ( Z_STATE .GT. WORKING_STATE ) FINAL_DZ_Z = DZ_Z
            }

            // Exit if both normwise and componentwise stopped working,
            // but if componentwise is unstable, let it go at least two
            // iterations.

            if ( X_STATE.NE.WORKING_STATE ) {
               IF ( IGNORE_CWISE ) GOTO 666
               IF ( Z_STATE.EQ.NOPROG_STATE .OR. Z_STATE.EQ.CONV_STATE ) GOTO 666
               IF ( Z_STATE.EQ.UNSTABLE_STATE .AND. CNT.GT.1 ) GOTO 666
            }

            if ( INCR_PREC ) {
               INCR_PREC = .FALSE.
               Y_PREC_STATE = Y_PREC_STATE + 1
               for (I = 1; I <= N; I++) {
                  Y_TAIL( I ) = 0.0
               END DO
            }

            PREVNORMDX = NORMDX
            PREV_DZ_Z = DZ_Z

            // Update solution.

            if ( Y_PREC_STATE .LT. EXTRA_Y ) {
               caxpy(N, (1.0E+0,0.0E+0), DY, 1, Y(1,J), 1 );
            } else {
               cla_wwaddw(N, Y( 1, J ), Y_TAIL, DY );
            }

         END DO
         // Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT.
 666     CONTINUE

      // Set final_* when cnt hits ithresh

         IF ( X_STATE .EQ. WORKING_STATE ) FINAL_DX_X = DX_X
         IF ( Z_STATE .EQ. WORKING_STATE ) FINAL_DZ_Z = DZ_Z

      // Compute error bounds

         if (N_NORMS .GE. 1) {
            ERRS_N( J, LA_LINRX_ERR_I ) = FINAL_DX_X / (1 - DXRATMAX)

         }
         if ( N_NORMS .GE. 2 ) {
            ERRS_C( J, LA_LINRX_ERR_I ) = FINAL_DZ_Z / (1 - DZRATMAX)
         }

      // Compute componentwise relative backward error from formula
          // max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
      // where abs(Z) is the componentwise absolute value of the matrix
      // or vector Z.

         // Compute residual RES = B_s - op(A_s) * Y,
             // op(A) = A, A**T, or A**H depending on TRANS (and type).

         ccopy(N, B( 1, J ), 1, RES, 1 );
         cgemv(TRANS, N, N, (-1.0E+0,0.0E+0), A, LDA, Y(1,J), 1, (1.0E+0,0.0E+0), RES, 1 );

         for (I = 1; I <= N; I++) {
            AYB( I ) = CABS1( B( I, J ) )
         END DO

      // Compute abs(op(A_s))*abs(Y) + abs(B_s).

         cla_geamv(TRANS_TYPE, N, N, 1.0E+0, A, LDA, Y(1, J), 1, 1.0E+0, AYB, 1 );

         cla_lin_berr(N, N, 1, RES, AYB, BERR_OUT( J ) );

      // End of loop for each RHS.

      END DO

      RETURN

      // End of CLA_GERFSX_EXTENDED

      }
