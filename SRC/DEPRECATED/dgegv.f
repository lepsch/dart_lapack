      SUBROUTINE DGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               ILIMIT, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, ITAU, IWORK, JC, JR, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      double             ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM, BNRM1, BNRM2, EPS, ONEPLS, SAFMAX, SAFMIN, SALFAI, SALFAR, SBETA, SCALE, TEMP;
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, DLASCL, DLASET, DORGQR, DORMQR, DTGEVC, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, MAX
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

      if ( LSAME( JOBVL, 'N' ) ) {
         IJOBVL = 1
         ILVL = .FALSE.
      } else if ( LSAME( JOBVL, 'V' ) ) {
         IJOBVL = 2
         ILVL = .TRUE.
      } else {
         IJOBVL = -1
         ILVL = .FALSE.
      }

      if ( LSAME( JOBVR, 'N' ) ) {
         IJOBVR = 1
         ILVR = .FALSE.
      } else if ( LSAME( JOBVR, 'V' ) ) {
         IJOBVR = 2
         ILVR = .TRUE.
      } else {
         IJOBVR = -1
         ILVR = .FALSE.
      }
      ILV = ILVL .OR. ILVR

      // Test the input arguments

      LWKMIN = MAX( 8*N, 1 )
      LWKOPT = LWKMIN
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      INFO = 0
      if ( IJOBVL.LE.0 ) {
         INFO = -1
      } else if ( IJOBVR.LE.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) {
         INFO = -12
      } else if ( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) {
         INFO = -14
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -16
      }

      if ( INFO.EQ.0 ) {
         NB1 = ILAENV( 1, 'DGEQRF', ' ', N, N, -1, -1 )
         NB2 = ILAENV( 1, 'DORMQR', ' ', N, N, N, -1 )
         NB3 = ILAENV( 1, 'DORGQR', ' ', N, N, N, -1 )
         NB = MAX( NB1, NB2, NB3 )
         LOPT = 2*N + MAX( 6*N, N*( NB+1 ) )
         WORK( 1 ) = LOPT
      }

      if ( INFO.NE.0 ) {
         xerbla('DGEGV ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      // Get machine constants

      EPS = DLAMCH( 'E' )*DLAMCH( 'B' )
      SAFMIN = DLAMCH( 'S' )
      SAFMIN = SAFMIN + SAFMIN
      SAFMAX = ONE / SAFMIN
      ONEPLS = ONE + ( 4*EPS )

      // Scale A

      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ANRM1 = ANRM
      ANRM2 = ONE
      if ( ANRM.LT.ONE ) {
         if ( SAFMAX*ANRM.LT.ONE ) {
            ANRM1 = SAFMIN
            ANRM2 = SAFMAX*ANRM
         }
      }

      if ( ANRM.GT.ZERO ) {
         dlascl('G', -1, -1, ANRM, ONE, N, N, A, LDA, IINFO );
         if ( IINFO.NE.0 ) {
            INFO = N + 10
            RETURN
         }
      }

      // Scale B

      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      BNRM1 = BNRM
      BNRM2 = ONE
      if ( BNRM.LT.ONE ) {
         if ( SAFMAX*BNRM.LT.ONE ) {
            BNRM1 = SAFMIN
            BNRM2 = SAFMAX*BNRM
         }
      }

      if ( BNRM.GT.ZERO ) {
         dlascl('G', -1, -1, BNRM, ONE, N, N, B, LDB, IINFO );
         if ( IINFO.NE.0 ) {
            INFO = N + 10
            RETURN
         }
      }

      // Permute the matrix to make it more nearly triangular
      // Workspace layout:  (8*N words -- "work" requires 6*N words)
         // left_permutation, right_permutation, work...

      ILEFT = 1
      IRIGHT = N + 1
      IWORK = IRIGHT + N
      dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWORK ), IINFO );
      if ( IINFO.NE.0 ) {
         INFO = N + 1
         GO TO 120
      }

      // Reduce B to triangular form, and initialize VL and/or VR
      // Workspace layout:  ("work..." must have at least N words)
         // left_permutation, right_permutation, tau, work...

      IROWS = IHI + 1 - ILO
      if ( ILV ) {
         ICOLS = N + 1 - ILO
      } else {
         ICOLS = IROWS
      }
      ITAU = IWORK
      IWORK = ITAU + IROWS
      dgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO.NE.0 ) {
         INFO = N + 2
         GO TO 120
      }

      dormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO.GE.0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO.NE.0 ) {
         INFO = N + 3
         GO TO 120
      }

      if ( ILVL ) {
         dlaset('Full', N, N, ZERO, ONE, VL, LDVL );
         dlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         dorgqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO );
         if (IINFO.GE.0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
         if ( IINFO.NE.0 ) {
            INFO = N + 4
            GO TO 120
         }
      }

      if (ILVR) CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR );

      // Reduce to generalized Hessenberg form

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         dgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IINFO );
      } else {
         dgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IINFO );
      }
      if ( IINFO.NE.0 ) {
         INFO = N + 5
         GO TO 120
      }

      // Perform QZ algorithm
      // Workspace layout:  ("work..." must have at least 1 word)
         // left_permutation, right_permutation, work...

      IWORK = ITAU
      if ( ILV ) {
         CHTEMP = 'S'
      } else {
         CHTEMP = 'E'
      }
      dhgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO.GE.0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO.NE.0 ) {
         if ( IINFO.GT.0 .AND. IINFO.LE.N ) {
            INFO = IINFO
         } else if ( IINFO.GT.N .AND. IINFO.LE.2*N ) {
            INFO = IINFO - N
         } else {
            INFO = N + 6
         }
         GO TO 120
      }

      if ( ILV ) {

         // Compute Eigenvectors  (DTGEVC requires 6*N words of workspace)

         if ( ILVL ) {
            if ( ILVR ) {
               CHTEMP = 'B'
            } else {
               CHTEMP = 'L'
            }
         } else {
            CHTEMP = 'R'
         }

         dtgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWORK ), IINFO );
         if ( IINFO.NE.0 ) {
            INFO = N + 7
            GO TO 120
         }

         // Undo balancing on VL and VR, rescale

         if ( ILVL ) {
            dggbak('P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VL, LDVL, IINFO );
            if ( IINFO.NE.0 ) {
               INFO = N + 8
               GO TO 120
            }
            for (JC = 1; JC <= N; JC++) { // 50
               IF( ALPHAI( JC ).LT.ZERO ) GO TO 50
               TEMP = ZERO
               if ( ALPHAI( JC ).EQ.ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 10
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
                  } // 10
               } else {
                  for (JR = 1; JR <= N; JR++) { // 20
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ ABS( VL( JR, JC+1 ) ) )
                  } // 20
               }
               if (TEMP.LT.SAFMIN) GO TO 50;
               TEMP = ONE / TEMP
               if ( ALPHAI( JC ).EQ.ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 30
                     VL( JR, JC ) = VL( JR, JC )*TEMP
                  } // 30
               } else {
                  for (JR = 1; JR <= N; JR++) { // 40
                     VL( JR, JC ) = VL( JR, JC )*TEMP
                     VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
                  } // 40
               }
            } // 50
         }
         if ( ILVR ) {
            dggbak('P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VR, LDVR, IINFO );
            if ( IINFO.NE.0 ) {
               INFO = N + 9
               GO TO 120
            }
            for (JC = 1; JC <= N; JC++) { // 100
               IF( ALPHAI( JC ).LT.ZERO ) GO TO 100
               TEMP = ZERO
               if ( ALPHAI( JC ).EQ.ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 60
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
                  } // 60
               } else {
                  for (JR = 1; JR <= N; JR++) { // 70
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+ ABS( VR( JR, JC+1 ) ) )
                  } // 70
               }
               if (TEMP.LT.SAFMIN) GO TO 100;
               TEMP = ONE / TEMP
               if ( ALPHAI( JC ).EQ.ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 80
                     VR( JR, JC ) = VR( JR, JC )*TEMP
                  } // 80
               } else {
                  for (JR = 1; JR <= N; JR++) { // 90
                     VR( JR, JC ) = VR( JR, JC )*TEMP
                     VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
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
         ABSAR = ABS( ALPHAR( JC ) )
         ABSAI = ABS( ALPHAI( JC ) )
         ABSB = ABS( BETA( JC ) )
         SALFAR = ANRM*ALPHAR( JC )
         SALFAI = ANRM*ALPHAI( JC )
         SBETA = BNRM*BETA( JC )
         ILIMIT = .FALSE.
         SCALE = ONE

         // Check for significant underflow in ALPHAI

         if ( ABS( SALFAI ).LT.SAFMIN .AND. ABSAI.GE. MAX( SAFMIN, EPS*ABSAR, EPS*ABSB ) ) {
            ILIMIT = .TRUE.
            SCALE = ( ONEPLS*SAFMIN / ANRM1 ) / MAX( ONEPLS*SAFMIN, ANRM2*ABSAI )

         } else if ( SALFAI.EQ.ZERO ) {

            // If insignificant underflow in ALPHAI, then make the
            // conjugate eigenvalue real.

            if ( ALPHAI( JC ).LT.ZERO .AND. JC.GT.1 ) {
               ALPHAI( JC-1 ) = ZERO
            } else if ( ALPHAI( JC ).GT.ZERO .AND. JC.LT.N ) {
               ALPHAI( JC+1 ) = ZERO
            }
         }

         // Check for significant underflow in ALPHAR

         if ( ABS( SALFAR ).LT.SAFMIN .AND. ABSAR.GE. MAX( SAFMIN, EPS*ABSAI, EPS*ABSB ) ) {
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / ANRM1 ) / MAX( ONEPLS*SAFMIN, ANRM2*ABSAR ) )
         }

         // Check for significant underflow in BETA

         if ( ABS( SBETA ).LT.SAFMIN .AND. ABSB.GE. MAX( SAFMIN, EPS*ABSAR, EPS*ABSAI ) ) {
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / BNRM1 ) / MAX( ONEPLS*SAFMIN, BNRM2*ABSB ) )
         }

         // Check for possible overflow when limiting scaling

         if ( ILIMIT ) {
            TEMP = ( SCALE*SAFMIN )*MAX( ABS( SALFAR ), ABS( SALFAI ), ABS( SBETA ) )             IF( TEMP.GT.ONE ) SCALE = SCALE / TEMP             IF( SCALE.LT.ONE ) ILIMIT = .FALSE.
         }

         // Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary.

         if ( ILIMIT ) {
            SALFAR = ( SCALE*ALPHAR( JC ) )*ANRM
            SALFAI = ( SCALE*ALPHAI( JC ) )*ANRM
            SBETA = ( SCALE*BETA( JC ) )*BNRM
         }
         ALPHAR( JC ) = SALFAR
         ALPHAI( JC ) = SALFAI
         BETA( JC ) = SBETA
      } // 110

      } // 120
      WORK( 1 ) = LWKOPT

      RETURN

      // End of DGEGV

      }
