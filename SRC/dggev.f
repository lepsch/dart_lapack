      SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

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
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, ITAU, IWRK, JC, JR, MAXWRK, MINWRK;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP;
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
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

      if ( LSAME( JOBVL, 'N' ) ) {
         IJOBVL = 1
         ILVL = false;
      } else if ( LSAME( JOBVL, 'V' ) ) {
         IJOBVL = 2
         ILVL = true;
      } else {
         IJOBVL = -1
         ILVL = false;
      }

      if ( LSAME( JOBVR, 'N' ) ) {
         IJOBVR = 1
         ILVR = false;
      } else if ( LSAME( JOBVR, 'V' ) ) {
         IJOBVR = 2
         ILVR = true;
      } else {
         IJOBVR = -1
         ILVR = false;
      }
      ILV = ILVL .OR. ILVR

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
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
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV. The workspace is
        // computed assuming ILO = 1 and IHI = N, the worst case.)

      if ( INFO.EQ.0 ) {
         MINWRK = MAX( 1, 8*N )
         MAXWRK = MAX( 1, N*( 7 + ILAENV( 1, 'DGEQRF', ' ', N, 1, N, 0 ) ) )          MAXWRK = MAX( MAXWRK, N*( 7 + ILAENV( 1, 'DORMQR', ' ', N, 1, N, 0 ) ) )
         if ( ILVL ) {
            MAXWRK = MAX( MAXWRK, N*( 7 + ILAENV( 1, 'DORGQR', ' ', N, 1, N, -1 ) ) )
         }
         WORK( 1 ) = MAXWRK

         if (LWORK.LT.MINWRK .AND. .NOT.LQUERY) INFO = -16;
      }

      if ( INFO.NE.0 ) {
         xerbla('DGGEV ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ILASCL = false;
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = true;
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = true;
      }
      if (ILASCL) CALL DLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      ILBSCL = false;
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = true;
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = true;
      }
      if (ILBSCL) CALL DLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute the matrices A, B to isolate eigenvalues if possible
      // (Workspace: need 6*N)

      ILEFT = 1
      IRIGHT = N + 1
      IWRK = IRIGHT + N
      dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)
      // (Workspace: need N, prefer N*NB)

      IROWS = IHI + 1 - ILO
      if ( ILV ) {
         ICOLS = N + 1 - ILO
      } else {
         ICOLS = IROWS
      }
      ITAU = IWRK
      IWRK = ITAU + IROWS
      dgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to matrix A
      // (Workspace: need N, prefer N*NB)

      dormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL
      // (Workspace: need N, prefer N*NB)

      if ( ILVL ) {
         dlaset('Full', N, N, ZERO, ONE, VL, LDVL );
         if ( IROWS.GT.1 ) {
            dlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         dorgqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VR

      if (ILVR) CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR );

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         dgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IERR );
      } else {
         dgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur forms and Schur vectors)
      // (Workspace: need N)

      IWRK = ITAU
      if ( ILV ) {
         CHTEMP = 'S'
      } else {
         CHTEMP = 'E'
      }
      dhgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR );
      if ( IERR.NE.0 ) {
         if ( IERR.GT.0 .AND. IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N .AND. IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 110
      }

      // Compute Eigenvectors
      // (Workspace: need 6*N)

      if ( ILV ) {
         if ( ILVL ) {
            if ( ILVR ) {
               CHTEMP = 'B'
            } else {
               CHTEMP = 'L'
            }
         } else {
            CHTEMP = 'R'
         }
         dtgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), IERR );
         if ( IERR.NE.0 ) {
            INFO = N + 2
            GO TO 110
         }

         // Undo balancing on VL and VR and normalization
         // (Workspace: none needed)

         if ( ILVL ) {
            dggbak('P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VL, LDVL, IERR );
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
               if (TEMP.LT.SMLNUM) GO TO 50;
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
            dggbak('P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VR, LDVR, IERR );
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
               if (TEMP.LT.SMLNUM) GO TO 100;
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

      // Undo scaling if necessary

      } // 110

      if ( ILASCL ) {
         dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
         dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
      }

      if ( ILBSCL ) {
         dlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      WORK( 1 ) = MAXWRK
      RETURN

      // End of DGGEV

      }
