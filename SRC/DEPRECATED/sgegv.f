      SUBROUTINE SGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO );

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               ILIMIT, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, ITAU, IWORK, JC, JR, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      REAL               ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM, BNRM1, BNRM2, EPS, ONEPLS, SAFMAX, SAFMIN, SALFAI, SALFAR, SBETA, SCALE, TEMP;
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRF, SGGBAK, SGGBAL, SGGHRD, SHGEQZ, SLACPY, SLASCL, SLASET, SORGQR, SORMQR, STGEVC, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANGE;
      // EXTERNAL ILAENV, LSAME, SLAMCH, SLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, MAX
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

      if ( LSAME( JOBVL, 'N' ) ) {
         IJOBVL = 1;
         ILVL = false;
      } else if ( LSAME( JOBVL, 'V' ) ) {
         IJOBVL = 2;
         ILVL = true;
      } else {
         IJOBVL = -1;
         ILVL = false;
      }

      if ( LSAME( JOBVR, 'N' ) ) {
         IJOBVR = 1;
         ILVR = false;
      } else if ( LSAME( JOBVR, 'V' ) ) {
         IJOBVR = 2;
         ILVR = true;
      } else {
         IJOBVR = -1;
         ILVR = false;
      }
      ILV = ILVL || ILVR;

      // Test the input arguments

      LWKMIN = MAX( 8*N, 1 );
      LWKOPT = LWKMIN;
      WORK( 1 ) = LWKOPT;
      LQUERY = ( LWORK == -1 );
      INFO = 0;
      if ( IJOBVL <= 0 ) {
         INFO = -1;
      } else if ( IJOBVR <= 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7;
      } else if ( LDVL < 1 || ( ILVL && LDVL < N ) ) {
         INFO = -12;
      } else if ( LDVR < 1 || ( ILVR && LDVR < N ) ) {
         INFO = -14;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -16;
      }

      if ( INFO == 0 ) {
         NB1 = ILAENV( 1, 'SGEQRF', ' ', N, N, -1, -1 );
         NB2 = ILAENV( 1, 'SORMQR', ' ', N, N, N, -1 );
         NB3 = ILAENV( 1, 'SORGQR', ' ', N, N, N, -1 );
         NB = MAX( NB1, NB2, NB3 );
         LOPT = 2*N + MAX( 6*N, N*(NB+1) );
         WORK( 1 ) = LOPT;
      }

      if ( INFO != 0 ) {
         xerbla('SGEGV ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Get machine constants

      EPS = SLAMCH( 'E' )*SLAMCH( 'B' );
      SAFMIN = SLAMCH( 'S' );
      SAFMIN = SAFMIN + SAFMIN;
      SAFMAX = ONE / SAFMIN;
      ONEPLS = ONE + ( 4*EPS );

      // Scale A

      ANRM = SLANGE( 'M', N, N, A, LDA, WORK );
      ANRM1 = ANRM;
      ANRM2 = ONE;
      if ( ANRM < ONE ) {
         if ( SAFMAX*ANRM < ONE ) {
            ANRM1 = SAFMIN;
            ANRM2 = SAFMAX*ANRM;
         }
      }

      if ( ANRM > ZERO ) {
         slascl('G', -1, -1, ANRM, ONE, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 10;
            return;
         }
      }

      // Scale B

      BNRM = SLANGE( 'M', N, N, B, LDB, WORK );
      BNRM1 = BNRM;
      BNRM2 = ONE;
      if ( BNRM < ONE ) {
         if ( SAFMAX*BNRM < ONE ) {
            BNRM1 = SAFMIN;
            BNRM2 = SAFMAX*BNRM;
         }
      }

      if ( BNRM > ZERO ) {
         slascl('G', -1, -1, BNRM, ONE, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 10;
            return;
         }
      }

      // Permute the matrix to make it more nearly triangular
      // Workspace layout:  (8*N words -- "work" requires 6*N words)
         // left_permutation, right_permutation, work...

      ILEFT = 1;
      IRIGHT = N + 1;
      IWORK = IRIGHT + N;
      sggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWORK ), IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 1;
         GO TO 120;
      }

      // Reduce B to triangular form, and initialize VL and/or VR
      // Workspace layout:  ("work..." must have at least N words)
         // left_permutation, right_permutation, tau, work...

      IROWS = IHI + 1 - ILO;
      if ( ILV ) {
         ICOLS = N + 1 - ILO;
      } else {
         ICOLS = IROWS;
      }
      ITAU = IWORK;
      IWORK = ITAU + IROWS;
      sgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO >= 0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 2;
         GO TO 120;
      }

      sormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 3;
         GO TO 120;
      }

      if ( ILVL ) {
         slaset('Full', N, N, ZERO, ONE, VL, LDVL );
         slacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         sorgqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO );
         if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
         if ( IINFO != 0 ) {
            INFO = N + 4;
            GO TO 120;
         }
      }

      if (ILVR) CALL SLASET( 'Full', N, N, ZERO, ONE, VR, LDVR );

      // Reduce to generalized Hessenberg form

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         sgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IINFO );
      } else {
         sgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IINFO );
      }
      if ( IINFO != 0 ) {
         INFO = N + 5;
         GO TO 120;
      }

      // Perform QZ algorithm
      // Workspace layout:  ("work..." must have at least 1 word)
         // left_permutation, right_permutation, work...

      IWORK = ITAU;
      if ( ILV ) {
         CHTEMP = 'S';
      } else {
         CHTEMP = 'E';
      }
      shgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
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

         // Compute Eigenvectors  (STGEVC requires 6*N words of workspace)

         if ( ILVL ) {
            if ( ILVR ) {
               CHTEMP = 'B';
            } else {
               CHTEMP = 'L';
            }
         } else {
            CHTEMP = 'R';
         }

         stgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWORK ), IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 7;
            GO TO 120;
         }

         // Undo balancing on VL and VR, rescale

         if ( ILVL ) {
            sggbak('P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VL, LDVL, IINFO );
            if ( IINFO != 0 ) {
               INFO = N + 8;
               GO TO 120;
            }
            for (JC = 1; JC <= N; JC++) { // 50
               IF( ALPHAI( JC ) < ZERO ) GO TO 50;
               TEMP = ZERO;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 10
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) );
                  } // 10
               } else {
                  for (JR = 1; JR <= N; JR++) { // 20
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ ABS( VL( JR, JC+1 ) ) );
                  } // 20
               }
               if (TEMP < SAFMIN) GO TO 50;
               TEMP = ONE / TEMP;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 30
                     VL( JR, JC ) = VL( JR, JC )*TEMP;
                  } // 30
               } else {
                  for (JR = 1; JR <= N; JR++) { // 40
                     VL( JR, JC ) = VL( JR, JC )*TEMP;
                     VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP;
                  } // 40
               }
            } // 50
         }
         if ( ILVR ) {
            sggbak('P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VR, LDVR, IINFO );
            if ( IINFO != 0 ) {
               INFO = N + 9;
               GO TO 120;
            }
            for (JC = 1; JC <= N; JC++) { // 100
               IF( ALPHAI( JC ) < ZERO ) GO TO 100;
               TEMP = ZERO;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 60
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) );
                  } // 60
               } else {
                  for (JR = 1; JR <= N; JR++) { // 70
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+ ABS( VR( JR, JC+1 ) ) );
                  } // 70
               }
               if (TEMP < SAFMIN) GO TO 100;
               TEMP = ONE / TEMP;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 80
                     VR( JR, JC ) = VR( JR, JC )*TEMP;
                  } // 80
               } else {
                  for (JR = 1; JR <= N; JR++) { // 90
                     VR( JR, JC ) = VR( JR, JC )*TEMP;
                     VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP;
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
         ABSAR = ABS( ALPHAR( JC ) );
         ABSAI = ABS( ALPHAI( JC ) );
         ABSB = ABS( BETA( JC ) );
         SALFAR = ANRM*ALPHAR( JC );
         SALFAI = ANRM*ALPHAI( JC );
         SBETA = BNRM*BETA( JC );
         ILIMIT = false;
         SCALE = ONE;

         // Check for significant underflow in ALPHAI

         if ( ABS( SALFAI ) < SAFMIN && ABSAI >= MAX( SAFMIN, EPS*ABSAR, EPS*ABSB ) ) {
            ILIMIT = true;
            SCALE = ( ONEPLS*SAFMIN / ANRM1 ) / MAX( ONEPLS*SAFMIN, ANRM2*ABSAI );

         } else if ( SALFAI == ZERO ) {

            // If insignificant underflow in ALPHAI, then make the
            // conjugate eigenvalue real.

            if ( ALPHAI( JC ) < ZERO && JC > 1 ) {
               ALPHAI( JC-1 ) = ZERO;
            } else if ( ALPHAI( JC ) > ZERO && JC < N ) {
               ALPHAI( JC+1 ) = ZERO;
            }
         }

         // Check for significant underflow in ALPHAR

         if ( ABS( SALFAR ) < SAFMIN && ABSAR >= MAX( SAFMIN, EPS*ABSAI, EPS*ABSB ) ) {
            ILIMIT = true;
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / ANRM1 ) / MAX( ONEPLS*SAFMIN, ANRM2*ABSAR ) );
         }

         // Check for significant underflow in BETA

         if ( ABS( SBETA ) < SAFMIN && ABSB >= MAX( SAFMIN, EPS*ABSAR, EPS*ABSAI ) ) {
            ILIMIT = true;
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / BNRM1 ) / MAX( ONEPLS*SAFMIN, BNRM2*ABSB ) );
         }

         // Check for possible overflow when limiting scaling

         if ( ILIMIT ) {
            TEMP = ( SCALE*SAFMIN )*MAX( ABS( SALFAR ), ABS( SALFAI ), ABS( SBETA ) )             IF( TEMP > ONE ) SCALE = SCALE / TEMP             IF( SCALE < ONE ) ILIMIT = false;
         }

         // Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary.

         if ( ILIMIT ) {
            SALFAR = ( SCALE*ALPHAR( JC ) )*ANRM;
            SALFAI = ( SCALE*ALPHAI( JC ) )*ANRM;
            SBETA = ( SCALE*BETA( JC ) )*BNRM;
         }
         ALPHAR( JC ) = SALFAR;
         ALPHAI( JC ) = SALFAI;
         BETA( JC ) = SBETA;
      } // 110

      } // 120
      WORK( 1 ) = LWKOPT;

      return;

      // End of SGEGV

      }
