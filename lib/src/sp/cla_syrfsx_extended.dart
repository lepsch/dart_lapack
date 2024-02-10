      void cla_syrfsx_extended(PREC_TYPE, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, COLEQU, C, B, LDB, Y, LDY, BERR_OUT, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, RES, AYB, DY, Y_TAIL, RCOND, ITHRESH, RTHRESH, DZ_UB, IGNORE_CWISE, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE, N_NORMS, ITHRESH;
      String             UPLO;
      bool               COLEQU, IGNORE_CWISE;
      double               RTHRESH, DZ_UB;
      int                IPIV( * );
      Complex            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )       double               C( * ), AYB( * ), RCOND, BERR_OUT( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                UPLO2, CNT, I, J, X_STATE, Z_STATE, Y_PREC_STATE;
      double               YK, DYK, YMIN, NORMY, NORMX, NORMDX, DXRAT, DZRAT, PREVNORMDX, PREV_DZ_Z, DXRATMAX, DZRATMAX, DX_X, DZ_Z, FINAL_DX_X, FINAL_DZ_Z, EPS, HUGEVAL, INCR_THRESH;
      bool               INCR_PREC, UPPER;
      Complex            ZDUM;
      // ..
      // .. Parameters ..
      int                UNSTABLE_STATE, WORKING_STATE, CONV_STATE, NOPROG_STATE, BASE_RESIDUAL, EXTRA_RESIDUAL, EXTRA_Y;
      const              UNSTABLE_STATE = 0, WORKING_STATE = 1, CONV_STATE = 2, NOPROG_STATE = 3 );
      PARAMETER          ( BASE_RESIDUAL = 0, EXTRA_RESIDUAL = 1, EXTRA_Y = 2 ;
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
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL ILAUPLO
      int                ILAUPLO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CSYTRS, CSYMV, BLAS_CSYMV_X, BLAS_CSYMV2_X, CLA_SYAMV, CLA_WWADDW, CLA_LIN_BERR
      double               SLAMCH;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, AIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -8;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -13;
      } else if ( LDY < max( 1, N ) ) {
         INFO = -15;
      }
      if ( INFO != 0 ) {
         xerbla('CLA_SYRFSX_EXTENDED', -INFO );
         return;
      }
      EPS = SLAMCH( 'Epsilon' );
      HUGEVAL = SLAMCH( 'Overflow' );
      // Force HUGEVAL to Inf
      HUGEVAL = HUGEVAL * HUGEVAL;
      // Using HUGEVAL may lead to spurious underflows.
      INCR_THRESH = double( N ) * EPS;

      if ( lsame ( UPLO, 'L' ) ) {
         UPLO2 = ILAUPLO( 'L' );
      } else {
         UPLO2 = ILAUPLO( 'U' );
      }

      for (J = 1; J <= NRHS; J++) {
         Y_PREC_STATE = EXTRA_RESIDUAL;
         if ( Y_PREC_STATE == EXTRA_Y ) {
            for (I = 1; I <= N; I++) {
               Y_TAIL[I] = 0.0;
            }
         }

         DXRAT = 0.0;
         DXRATMAX = 0.0;
         DZRAT = 0.0;
         DZRATMAX = 0.0;
         FINAL_DX_X = HUGEVAL;
         FINAL_DZ_Z = HUGEVAL;
         PREVNORMDX = HUGEVAL;
         PREV_DZ_Z = HUGEVAL;
         DZ_Z = HUGEVAL;
         DX_X = HUGEVAL;

         X_STATE = WORKING_STATE;
         Z_STATE = UNSTABLE_STATE;
         INCR_PREC = false;

         for (CNT = 1; CNT <= ITHRESH; CNT++) {

          // Compute residual RES = B_s - op(A_s) * Y,
              // op(A) = A, A**T, or A**H depending on TRANS (and type).

            ccopy(N, B( 1, J ), 1, RES, 1 );
            if ( Y_PREC_STATE == BASE_RESIDUAL ) {
               csymv(UPLO, N, CMPLX(-1.0), A, LDA, Y(1,J), 1, CMPLX(1.0), RES, 1 );
            } else if ( Y_PREC_STATE == EXTRA_RESIDUAL ) {
               blas_csymv_x(UPLO2, N, CMPLX(-1.0), A, LDA, Y( 1, J ), 1, CMPLX(1.0), RES, 1, PREC_TYPE );
            } else {
               blas_csymv2_x(UPLO2, N, CMPLX(-1.0), A, LDA, Y(1, J), Y_TAIL, 1, CMPLX(1.0), RES, 1, PREC_TYPE);
            }

          // XXX: RES is no longer needed.
            ccopy(N, RES, 1, DY, 1 );
            csytrs(UPLO, N, 1, AF, LDAF, IPIV, DY, N, INFO );

          // Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.

            NORMX = 0.0;
            NORMY = 0.0;
            NORMDX = 0.0;
            DZ_Z = 0.0;
            YMIN = HUGEVAL;

            for (I = 1; I <= N; I++) {
               YK = CABS1( Y( I, J ) );
               DYK = CABS1( DY( I ) );

               if ( YK != 0.0 ) {
                  DZ_Z = max( DZ_Z, DYK / YK );
               } else if ( DYK != 0.0 ) {
                  DZ_Z = HUGEVAL;
               }

               YMIN = min( YMIN, YK );

               NORMY = max( NORMY, YK );

               if ( COLEQU ) {
                  NORMX = max( NORMX, YK * C( I ) );
                  NORMDX = max( NORMDX, DYK * C( I ) );
               } else {
                  NORMX = NORMY;
                  NORMDX = max( NORMDX, DYK );
               }
            }

            if ( NORMX != 0.0 ) {
               DX_X = NORMDX / NORMX;
            } else if ( NORMDX == 0.0 ) {
               DX_X = 0.0;
            } else {
               DX_X = HUGEVAL;
            }

            DXRAT = NORMDX / PREVNORMDX;
            DZRAT = DZ_Z / PREV_DZ_Z;

          // Check termination criteria.

            if (YMIN*RCOND < INCR_THRESH*NORMY && Y_PREC_STATE < EXTRA_Y) INCR_PREC = true ;
             if (X_STATE == NOPROG_STATE && DXRAT <= RTHRESH) X_STATE = WORKING_STATE;
            if ( X_STATE == WORKING_STATE ) {
               if ( DX_X <= EPS ) {
                  X_STATE = CONV_STATE;
               } else if ( DXRAT > RTHRESH ) {
                  if ( Y_PREC_STATE != EXTRA_Y ) {
                     INCR_PREC = true;
                  } else {
                     X_STATE = NOPROG_STATE;
                  }
               } else {
                  if (DXRAT > DXRATMAX) DXRATMAX = DXRAT;
               }
               if (X_STATE > WORKING_STATE) FINAL_DX_X = DX_X;
            }
             if (Z_STATE == UNSTABLE_STATE && DZ_Z <= DZ_UB) Z_STATE = WORKING_STATE;
            IF ( Z_STATE == NOPROG_STATE && DZRAT <= RTHRESH ) Z_STATE = WORKING_STATE;
            if ( Z_STATE == WORKING_STATE ) {
               if ( DZ_Z <= EPS ) {
                  Z_STATE = CONV_STATE;
               } else if ( DZ_Z > DZ_UB ) {
                  Z_STATE = UNSTABLE_STATE;
                  DZRATMAX = 0.0;
                  FINAL_DZ_Z = HUGEVAL;
               } else if ( DZRAT > RTHRESH ) {
                  if ( Y_PREC_STATE != EXTRA_Y ) {
                     INCR_PREC = true;
                  } else {
                     Z_STATE = NOPROG_STATE;
                  }
               } else {
                  if (DZRAT > DZRATMAX) DZRATMAX = DZRAT;
               }
               if (Z_STATE > WORKING_STATE) FINAL_DZ_Z = DZ_Z;
            }
             if ( X_STATE != WORKING_STATE && ( IGNORE_CWISE || Z_STATE != WORKING_STATE ) ) GOTO 666;

            if ( INCR_PREC ) {
               INCR_PREC = false;
               Y_PREC_STATE = Y_PREC_STATE + 1;
               for (I = 1; I <= N; I++) {
                  Y_TAIL[I] = 0.0;
               }
            }

            PREVNORMDX = NORMDX;
            PREV_DZ_Z = DZ_Z;

            // Update solution.

            if ( Y_PREC_STATE < EXTRA_Y ) {
               caxpy(N, CMPLX(1.0), DY, 1, Y(1,J), 1 );
            } else {
               cla_wwaddw(N, Y(1,J), Y_TAIL, DY );
            }

         }
         // Target of "IF (Z_STOP && X_STOP)".  Sun's f77 won't EXIT.
         } // 666

      // Set final_* when cnt hits ithresh.

         if (X_STATE == WORKING_STATE) FINAL_DX_X = DX_X;
         if (Z_STATE == WORKING_STATE) FINAL_DZ_Z = DZ_Z;

      // Compute error bounds.

         if ( N_NORMS >= 1 ) {
            ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = FINAL_DX_X / (1 - DXRATMAX);
         }
         if ( N_NORMS >= 2 ) {
            ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = FINAL_DZ_Z / (1 - DZRATMAX);
         }

      // Compute componentwise relative backward error from formula
          // max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
      // where abs(Z) is the componentwise absolute value of the matrix
      // or vector Z.

         // Compute residual RES = B_s - op(A_s) * Y,
             // op(A) = A, A**T, or A**H depending on TRANS (and type).

         ccopy(N, B( 1, J ), 1, RES, 1 );
         csymv(UPLO, N, CMPLX(-1.0), A, LDA, Y(1,J), 1, CMPLX(1.0), RES, 1 );

         for (I = 1; I <= N; I++) {
            AYB[I] = CABS1( B( I, J ) );
         }

      // Compute abs(op(A_s))*abs(Y) + abs(B_s).

         cla_syamv(UPLO2, N, 1.0, A, LDA, Y(1, J), 1, 1.0, AYB, 1 );

         cla_lin_berr(N, N, 1, RES, AYB, BERR_OUT( J ) );

      // End of loop for each RHS.

      }

      }
