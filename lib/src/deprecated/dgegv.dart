      void dgegv(final int JOBVL, final int JOBVR, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int ALPHAR, final int ALPHAI, final int BETA, final Matrix<double> VL_, final int LDVL, final Matrix<double> VR_, final int LDVR, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.having();
  final B = B_.having();
  final VL = VL_.having();
  final VR = VR_.having();
  final WORK = WORK_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               ILIMIT, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, ITAU, IWORK, JC, JR, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      double             ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM, BNRM1, BNRM2, EPS, ONEPLS, SAFMAX, SAFMIN, SALFAI, SALFAR, SBETA, SCALE, TEMP;
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, DLASCL, DLASET, DORGQR, DORMQR, DTGEVC, XERBLA
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, DLANGE;
      // EXTERNAL lsame, ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, MAX

      // Decode the input arguments

      if ( lsame( JOBVL, 'N' ) ) {
         IJOBVL = 1;
         ILVL = false;
      } else if ( lsame( JOBVL, 'V' ) ) {
         IJOBVL = 2;
         ILVL = true;
      } else {
         IJOBVL = -1;
         ILVL = false;
      }

      if ( lsame( JOBVR, 'N' ) ) {
         IJOBVR = 1;
         ILVR = false;
      } else if ( lsame( JOBVR, 'V' ) ) {
         IJOBVR = 2;
         ILVR = true;
      } else {
         IJOBVR = -1;
         ILVR = false;
      }
      ILV = ILVL || ILVR;

      // Test the input arguments

      LWKMIN = max( 8*N, 1 );
      LWKOPT = LWKMIN;
      WORK[1] = LWKOPT;
      LQUERY = ( LWORK == -1 );
      INFO = 0;
      if ( IJOBVL <= 0 ) {
         INFO = -1;
      } else if ( IJOBVR <= 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDVL < 1 || ( ILVL && LDVL < N ) ) {
         INFO = -12;
      } else if ( LDVR < 1 || ( ILVR && LDVR < N ) ) {
         INFO = -14;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -16;
      }

      if ( INFO == 0 ) {
         NB1 = ilaenv( 1, 'DGEQRF', ' ', N, N, -1, -1 );
         NB2 = ilaenv( 1, 'DORMQR', ' ', N, N, N, -1 );
         NB3 = ilaenv( 1, 'DORGQR', ' ', N, N, N, -1 );
         NB = max( NB1, NB2, NB3 );
         LOPT = 2*N + max( 6*N, N*( NB+1 ) );
         WORK[1] = LOPT;
      }

      if ( INFO != 0 ) {
         xerbla('DGEGV ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = dlamch( 'E' )*dlamch( 'B' );
      SAFMIN = dlamch( 'S' );
      SAFMIN += SAFMIN;
      SAFMAX = ONE / SAFMIN;
      ONEPLS = ONE + ( 4*EPS );

      // Scale A

      ANRM = dlange( 'M', N, N, A, LDA, WORK );
      ANRM1 = ANRM;
      ANRM2 = ONE;
      if ( ANRM < ONE ) {
         if ( SAFMAX*ANRM < ONE ) {
            ANRM1 = SAFMIN;
            ANRM2 = SAFMAX*ANRM;
         }
      }

      if ( ANRM > ZERO ) {
         dlascl('G', -1, -1, ANRM, ONE, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 10;
            return;
         }
      }

      // Scale B

      BNRM = dlange( 'M', N, N, B, LDB, WORK );
      BNRM1 = BNRM;
      BNRM2 = ONE;
      if ( BNRM < ONE ) {
         if ( SAFMAX*BNRM < ONE ) {
            BNRM1 = SAFMIN;
            BNRM2 = SAFMAX*BNRM;
         }
      }

      if ( BNRM > ZERO ) {
         dlascl('G', -1, -1, BNRM, ONE, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 10;
            return;
         }
      }

      // Permute the matrix to make it more nearly triangular
      // Workspace layout:  (8*N words -- "work" requires 6*N words)
      //    left_permutation, right_permutation, work...

      ILEFT = 1;
      IRIGHT = N + 1;
      IWORK = IRIGHT + N;
      dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWORK ), IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 1;
         GO TO 120;
      }

      // Reduce B to triangular form, and initialize VL and/or VR
      // Workspace layout:  ("work..." must have at least N words)
      //    left_permutation, right_permutation, tau, work...

      IROWS = IHI + 1 - ILO;
      if ( ILV ) {
         ICOLS = N + 1 - ILO;
      } else {
         ICOLS = IROWS;
      }
      ITAU = IWORK;
      IWORK = ITAU + IROWS;
      dgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO >= 0 ) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 2;
         GO TO 120;
      }

      dormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 3;
         GO TO 120;
      }

      if ( ILVL ) {
         dlaset('Full', N, N, ZERO, ONE, VL, LDVL );
         dlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         dorgqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO );
         if (IINFO >= 0) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
         if ( IINFO != 0 ) {
            INFO = N + 4;
            GO TO 120;
         }
      }

      if (ILVR) dlaset( 'Full', N, N, ZERO, ONE, VR, LDVR );

      // Reduce to generalized Hessenberg form

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         dgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IINFO );
      } else {
         dgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IINFO );
      }
      if ( IINFO != 0 ) {
         INFO = N + 5;
         GO TO 120;
      }

      // Perform QZ algorithm
      // Workspace layout:  ("work..." must have at least 1 word)
      //    left_permutation, right_permutation, work...

      IWORK = ITAU;
      if ( ILV ) {
         CHTEMP = 'S';
      } else {
         CHTEMP = 'E';
      }
      dhgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         if ( IINFO > 0 && IINFO <= N ) {
            INFO = IINFO;
         } else if ( IINFO > N && IINFO <= 2*N ) {
            INFO = IINFO - N;
         } else {
            INFO = N + 6;
         }
         GO TO 120;
      }

      if ( ILV ) {

         // Compute Eigenvectors  (DTGEVC requires 6*N words of workspace)

         if ( ILVL ) {
            if ( ILVR ) {
               CHTEMP = 'B';
            } else {
               CHTEMP = 'L';
            }
         } else {
            CHTEMP = 'R';
         }

         dtgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWORK ), IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 7;
            GO TO 120;
         }

         // Undo balancing on VL and VR, rescale

         if ( ILVL ) {
            dggbak('P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VL, LDVL, IINFO );
            if ( IINFO != 0 ) {
               INFO = N + 8;
               GO TO 120;
            }
            for (JC = 1; JC <= N; JC++) { // 50
               if( ALPHAI( JC ) < ZERO ) GO TO 50;
               TEMP = ZERO;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 10
                     TEMP = max( TEMP, ( VL( JR, JC ) ).abs() );
                  } // 10
               } else {
                  for (JR = 1; JR <= N; JR++) { // 20
                     TEMP = max( TEMP, ( VL( JR, JC ) ).abs()+ ( VL( JR, JC+1 ) ).abs() );
                  } // 20
               }
               if (TEMP < SAFMIN) GO TO 50;
               TEMP = ONE / TEMP;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 30
                     VL[JR][JC] = VL( JR, JC )*TEMP;
                  } // 30
               } else {
                  for (JR = 1; JR <= N; JR++) { // 40
                     VL[JR][JC] = VL( JR, JC )*TEMP;
                     VL[JR][JC+1] = VL( JR, JC+1 )*TEMP;
                  } // 40
               }
            } // 50
         }
         if ( ILVR ) {
            dggbak('P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VR, LDVR, IINFO );
            if ( IINFO != 0 ) {
               INFO = N + 9;
               GO TO 120;
            }
            for (JC = 1; JC <= N; JC++) { // 100
               if( ALPHAI( JC ) < ZERO ) GO TO 100;
               TEMP = ZERO;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 60
                     TEMP = max( TEMP, ( VR( JR, JC ) ).abs() );
                  } // 60
               } else {
                  for (JR = 1; JR <= N; JR++) { // 70
                     TEMP = max( TEMP, ( VR( JR, JC ) ).abs()+ ( VR( JR, JC+1 ) ).abs() );
                  } // 70
               }
               if (TEMP < SAFMIN) GO TO 100;
               TEMP = ONE / TEMP;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 80
                     VR[JR][JC] = VR( JR, JC )*TEMP;
                  } // 80
               } else {
                  for (JR = 1; JR <= N; JR++) { // 90
                     VR[JR][JC] = VR( JR, JC )*TEMP;
                     VR[JR][JC+1] = VR( JR, JC+1 )*TEMP;
                  } // 90
               }
            } // 100
         }

         // End of eigenvector calculation

      }

      // Undo scaling in alpha, beta

      // Note: this does not give the alpha and beta for the unscaled
      // problem.

      // Un-scaling is limited to avoid underflow in alpha and beta
      // if they are significant.

      for (JC = 1; JC <= N; JC++) { // 110
         ABSAR = ( ALPHAR( JC ) ).abs();
         ABSAI = ( ALPHAI( JC ) ).abs();
         ABSB = ( BETA( JC ) ).abs();
         SALFAR = ANRM*ALPHAR( JC );
         SALFAI = ANRM*ALPHAI( JC );
         SBETA = BNRM*BETA( JC );
         ILIMIT = false;
         SCALE = ONE;

         // Check for significant underflow in ALPHAI

         if ( ( SALFAI ).abs() < SAFMIN && ABSAI >= max( SAFMIN, EPS*ABSAR, EPS*ABSB ) ) {
            ILIMIT = true;
            SCALE = ( ONEPLS*SAFMIN / ANRM1 ) / max( ONEPLS*SAFMIN, ANRM2*ABSAI );

         } else if ( SALFAI == ZERO ) {

            // If insignificant underflow in ALPHAI, then make the
            // conjugate eigenvalue real.

            if ( ALPHAI( JC ) < ZERO && JC > 1 ) {
               ALPHAI[JC-1] = ZERO;
            } else if ( ALPHAI( JC ) > ZERO && JC < N ) {
               ALPHAI[JC+1] = ZERO;
            }
         }

         // Check for significant underflow in ALPHAR

         if ( ( SALFAR ).abs() < SAFMIN && ABSAR >= max( SAFMIN, EPS*ABSAI, EPS*ABSB ) ) {
            ILIMIT = true;
            SCALE = max( SCALE, ( ONEPLS*SAFMIN / ANRM1 ) / max( ONEPLS*SAFMIN, ANRM2*ABSAR ) );
         }

         // Check for significant underflow in BETA

         if ( ( SBETA ).abs() < SAFMIN && ABSB >= max( SAFMIN, EPS*ABSAR, EPS*ABSAI ) ) {
            ILIMIT = true;
            SCALE = max( SCALE, ( ONEPLS*SAFMIN / BNRM1 ) / max( ONEPLS*SAFMIN, BNRM2*ABSB ) );
         }

         // Check for possible overflow when limiting scaling

         if ( ILIMIT ) {
            TEMP = ( SCALE*SAFMIN )*max( ( SALFAR ).abs(), ( SALFAI ).abs(), ( SBETA ).abs() )             IF( TEMP > ONE ) SCALE = SCALE / TEMP             IF( SCALE < ONE ) ILIMIT = false;
         }

         // Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary.

         if ( ILIMIT ) {
            SALFAR = ( SCALE*ALPHAR( JC ) )*ANRM;
            SALFAI = ( SCALE*ALPHAI( JC ) )*ANRM;
            SBETA = ( SCALE*BETA( JC ) )*BNRM;
         }
         ALPHAR[JC] = SALFAR;
         ALPHAI[JC] = SALFAI;
         BETA[JC] = SBETA;
      } // 110

      } // 120
      WORK[1] = LWKOPT;

      }
