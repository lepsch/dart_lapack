      SUBROUTINE CGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, IRWRK, ITAU, IWRK, JC, JR, LWKMIN, LWKOPT;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP
      COMPLEX            X
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
      REAL               CLANGE, SLAMCH, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, CLANGE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               ABS1
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( REAL( X ) ) + ABS( AIMAG( X ) )
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
         INFO = -11
      } else if ( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) {
         INFO = -13
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV. The workspace is
        // computed assuming ILO = 1 and IHI = N, the worst case.)

      if ( INFO.EQ.0 ) {
         LWKMIN = MAX( 1, 2*N )
         LWKOPT = MAX( 1, N + N*ILAENV( 1, 'CGEQRF', ' ', N, 1, N, 0 ) )
         LWKOPT = MAX( LWKOPT, N + N*ILAENV( 1, 'CUNMQR', ' ', N, 1, N, 0 ) )
         if ( ILVL ) {
            LWKOPT = MAX( LWKOPT, N + N*ILAENV( 1, 'CUNGQR', ' ', N, 1, N, -1 ) )
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

         if (LWORK.LT.LWKMIN .AND. .NOT.LQUERY) INFO = -15;
      }

      if ( INFO.NE.0 ) {
         xerbla('CGGEV ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      // Get machine constants

      EPS = SLAMCH( 'E' )*SLAMCH( 'B' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = false;
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = true;
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = true;
      }
      if (ILASCL) CALL CLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = false;
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = true;
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = true;
      }
      if (ILBSCL) CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute the matrices A, B to isolate eigenvalues if possible
      // (Real Workspace: need 6*N)

      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      cggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)
      // (Complex Workspace: need N, prefer N*NB)

      IROWS = IHI + 1 - ILO
      if ( ILV ) {
         ICOLS = N + 1 - ILO
      } else {
         ICOLS = IROWS
      }
      ITAU = 1
      IWRK = ITAU + IROWS
      cgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to matrix A
      // (Complex Workspace: need N, prefer N*NB)

      cunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL
      // (Complex Workspace: need N, prefer N*NB)

      if ( ILVL ) {
         claset('Full', N, N, CZERO, CONE, VL, LDVL );
         if ( IROWS.GT.1 ) {
            clacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         cungqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VR

      if (ILVR) CALL CLASET( 'Full', N, N, CZERO, CONE, VR, LDVR );

      // Reduce to generalized Hessenberg form

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         cgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IERR );
      } else {
         cgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur form and Schur vectors)
      // (Complex Workspace: need N)
      // (Real Workspace: need N)

      IWRK = ITAU
      if ( ILV ) {
         CHTEMP = 'S'
      } else {
         CHTEMP = 'E'
      }
      chgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, RWORK( IRWRK ), IERR );
      if ( IERR.NE.0 ) {
         if ( IERR.GT.0 .AND. IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N .AND. IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 70
      }

      // Compute Eigenvectors
      // (Real Workspace: need 2*N)
      // (Complex Workspace: need 2*N)

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

         ctgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), RWORK( IRWRK ), IERR );
         if ( IERR.NE.0 ) {
            INFO = N + 2
            GO TO 70
         }

         // Undo balancing on VL and VR and normalization
         // (Workspace: none needed)

         if ( ILVL ) {
            cggbak('P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VL, LDVL, IERR );
            for (JC = 1; JC <= N; JC++) { // 30
               TEMP = ZERO
               for (JR = 1; JR <= N; JR++) { // 10
                  TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
               } // 10
               if (TEMP.LT.SMLNUM) GO TO 30;
               TEMP = ONE / TEMP
               for (JR = 1; JR <= N; JR++) { // 20
                  VL( JR, JC ) = VL( JR, JC )*TEMP
               } // 20
            } // 30
         }
         if ( ILVR ) {
            cggbak('P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VR, LDVR, IERR );
            for (JC = 1; JC <= N; JC++) { // 60
               TEMP = ZERO
               for (JR = 1; JR <= N; JR++) { // 40
                  TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
               } // 40
               if (TEMP.LT.SMLNUM) GO TO 60;
               TEMP = ONE / TEMP
               for (JR = 1; JR <= N; JR++) { // 50
                  VR( JR, JC ) = VR( JR, JC )*TEMP
               } // 50
            } // 60
         }
      }

      // Undo scaling if necessary

      } // 70

      if (ILASCL) CALL CLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR );

      if (ILBSCL) CALL CLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of CGGEV

      }
