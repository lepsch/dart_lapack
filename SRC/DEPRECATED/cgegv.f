      SUBROUTINE CGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILIMIT, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, IRWORK, ITAU, IWORK, JC, JR, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      REAL               ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM, BNRM1, BNRM2, EPS, SAFMAX, SAFMIN, SALFAI, SALFAR, SBETA, SCALE, TEMP;
      COMPLEX            X;
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CGGBAK, CGGBAL, CGGHRD, CHGEQZ, CLACPY, CLASCL, CLASET, CTGEVC, CUNGQR, CUNMQR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               CLANGE, SLAMCH;
      // EXTERNAL ILAENV, LSAME, CLANGE, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, INT, MAX, REAL
      // ..
      // .. Statement Functions ..
      REAL               ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( REAL( X ) ) + ABS( AIMAG( X ) );
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

      LWKMIN = MAX( 2*N, 1 );
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
         INFO = -11;
      } else if ( LDVR < 1 || ( ILVR && LDVR < N ) ) {
         INFO = -13;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -15;
      }

      if ( INFO == 0 ) {
         NB1 = ILAENV( 1, 'CGEQRF', ' ', N, N, -1, -1 );
         NB2 = ILAENV( 1, 'CUNMQR', ' ', N, N, N, -1 );
         NB3 = ILAENV( 1, 'CUNGQR', ' ', N, N, N, -1 );
         NB = MAX( NB1, NB2, NB3 );
         LOPT = MAX( 2*N, N*(NB+1) );
         WORK( 1 ) = LOPT;
      }

      if ( INFO != 0 ) {
         xerbla('CGEGV ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = SLAMCH( 'E' )*SLAMCH( 'B' );
      SAFMIN = SLAMCH( 'S' );
      SAFMIN = SAFMIN + SAFMIN;
      SAFMAX = ONE / SAFMIN;

      // Scale A

      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK );
      ANRM1 = ANRM;
      ANRM2 = ONE;
      if ( ANRM < ONE ) {
         if ( SAFMAX*ANRM < ONE ) {
            ANRM1 = SAFMIN;
            ANRM2 = SAFMAX*ANRM;
         }
      }

      if ( ANRM > ZERO ) {
         clascl('G', -1, -1, ANRM, ONE, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 10;
            return;
         }
      }

      // Scale B

      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK );
      BNRM1 = BNRM;
      BNRM2 = ONE;
      if ( BNRM < ONE ) {
         if ( SAFMAX*BNRM < ONE ) {
            BNRM1 = SAFMIN;
            BNRM2 = SAFMAX*BNRM;
         }
      }

      if ( BNRM > ZERO ) {
         clascl('G', -1, -1, BNRM, ONE, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 10;
            return;
         }
      }

      // Permute the matrix to make it more nearly triangular
      // Also "balance" the matrix.

      ILEFT = 1;
      IRIGHT = N + 1;
      IRWORK = IRIGHT + N;
      cggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWORK ), IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 1;
         GO TO 80;
      }

      // Reduce B to triangular form, and initialize VL and/or VR

      IROWS = IHI + 1 - ILO;
      if ( ILV ) {
         ICOLS = N + 1 - ILO;
      } else {
         ICOLS = IROWS;
      }
      ITAU = 1;
      IWORK = ITAU + IROWS;
      cgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO >= 0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 2;
         GO TO 80;
      }

      cunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 3;
         GO TO 80;
      }

      if ( ILVL ) {
         claset('Full', N, N, CZERO, CONE, VL, LDVL );
         clacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         cungqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO );
         if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
         if ( IINFO != 0 ) {
            INFO = N + 4;
            GO TO 80;
         }
      }

      if (ILVR) CALL CLASET( 'Full', N, N, CZERO, CONE, VR, LDVR );

      // Reduce to generalized Hessenberg form

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         cgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IINFO );
      } else {
         cgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IINFO );
      }
      if ( IINFO != 0 ) {
         INFO = N + 5;
         GO TO 80;
      }

      // Perform QZ algorithm

      IWORK = ITAU;
      if ( ILV ) {
         CHTEMP = 'S';
      } else {
         CHTEMP = 'E';
      }
      chgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWORK ), LWORK+1-IWORK, RWORK( IRWORK ), IINFO );
      if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         if ( IINFO > 0 && IINFO <= N ) {
            INFO = IINFO;
         } else if ( IINFO > N && IINFO <= 2*N ) {
            INFO = IINFO - N;
         } else {
            INFO = N + 6;
         }
         GO TO 80;
      }

      if ( ILV ) {

         // Compute Eigenvectors

         if ( ILVL ) {
            if ( ILVR ) {
               CHTEMP = 'B';
            } else {
               CHTEMP = 'L';
            }
         } else {
            CHTEMP = 'R';
         }

         ctgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWORK ), RWORK( IRWORK ), IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 7;
            GO TO 80;
         }

         // Undo balancing on VL and VR, rescale

         if ( ILVL ) {
            cggbak('P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VL, LDVL, IINFO );
            if ( IINFO != 0 ) {
               INFO = N + 8;
               GO TO 80;
            }
            for (JC = 1; JC <= N; JC++) { // 30
               TEMP = ZERO;
               for (JR = 1; JR <= N; JR++) { // 10
                  TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) );
               } // 10
               if (TEMP < SAFMIN) GO TO 30;
               TEMP = ONE / TEMP;
               for (JR = 1; JR <= N; JR++) { // 20
                  VL( JR, JC ) = VL( JR, JC )*TEMP;
               } // 20
            } // 30
         }
         if ( ILVR ) {
            cggbak('P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VR, LDVR, IINFO );
            if ( IINFO != 0 ) {
               INFO = N + 9;
               GO TO 80;
            }
            for (JC = 1; JC <= N; JC++) { // 60
               TEMP = ZERO;
               for (JR = 1; JR <= N; JR++) { // 40
                  TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) );
               } // 40
               if (TEMP < SAFMIN) GO TO 60;
               TEMP = ONE / TEMP;
               for (JR = 1; JR <= N; JR++) { // 50
                  VR( JR, JC ) = VR( JR, JC )*TEMP;
               } // 50
            } // 60
         }

         // End of eigenvector calculation

      }

      // Undo scaling in alpha, beta

      // Note: this does not give the alpha and beta for the unscaled
      // problem.

      // Un-scaling is limited to avoid underflow in alpha and beta
      // if they are significant.

      for (JC = 1; JC <= N; JC++) { // 70
         ABSAR = ABS( REAL( ALPHA( JC ) ) );
         ABSAI = ABS( AIMAG( ALPHA( JC ) ) );
         ABSB = ABS( REAL( BETA( JC ) ) );
         SALFAR = ANRM*REAL( ALPHA( JC ) );
         SALFAI = ANRM*AIMAG( ALPHA( JC ) );
         SBETA = BNRM*REAL( BETA( JC ) );
         ILIMIT = false;
         SCALE = ONE;

         // Check for significant underflow in imaginary part of ALPHA

         if ( ABS( SALFAI ) < SAFMIN && ABSAI >= MAX( SAFMIN, EPS*ABSAR, EPS*ABSB ) ) {
            ILIMIT = true;
            SCALE = ( SAFMIN / ANRM1 ) / MAX( SAFMIN, ANRM2*ABSAI );
         }

         // Check for significant underflow in real part of ALPHA

         if ( ABS( SALFAR ) < SAFMIN && ABSAR >= MAX( SAFMIN, EPS*ABSAI, EPS*ABSB ) ) {
            ILIMIT = true;
            SCALE = MAX( SCALE, ( SAFMIN / ANRM1 ) / MAX( SAFMIN, ANRM2*ABSAR ) );
         }

         // Check for significant underflow in BETA

         if ( ABS( SBETA ) < SAFMIN && ABSB >= MAX( SAFMIN, EPS*ABSAR, EPS*ABSAI ) ) {
            ILIMIT = true;
            SCALE = MAX( SCALE, ( SAFMIN / BNRM1 ) / MAX( SAFMIN, BNRM2*ABSB ) );
         }

         // Check for possible overflow when limiting scaling

         if ( ILIMIT ) {
            TEMP = ( SCALE*SAFMIN )*MAX( ABS( SALFAR ), ABS( SALFAI ), ABS( SBETA ) )             IF( TEMP > ONE ) SCALE = SCALE / TEMP             IF( SCALE < ONE ) ILIMIT = false;
         }

         // Recompute un-scaled ALPHA, BETA if necessary.

         if ( ILIMIT ) {
            SALFAR = ( SCALE*REAL( ALPHA( JC ) ) )*ANRM;
            SALFAI = ( SCALE*AIMAG( ALPHA( JC ) ) )*ANRM;
            SBETA = ( SCALE*BETA( JC ) )*BNRM;
         }
         ALPHA( JC ) = CMPLX( SALFAR, SALFAI );
         BETA( JC ) = SBETA;
      } // 70

      } // 80
      WORK( 1 ) = LWKOPT;

      return;

      // End of CGEGV

      }
