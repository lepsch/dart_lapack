      SUBROUTINE CGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, ILO, IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      REAL               ABNRM, BBNRM
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      REAL               LSCALE( * ), RCONDE( * ), RCONDV( * ), RSCALE( * ), RWORK( * )       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY, NOSCL, WANTSB, WANTSE, WANTSN, WANTSV;
      String             CHTEMP;
      int                I, ICOLS, IERR, IJOBVL, IJOBVR, IN, IROWS, ITAU, IWRK, IWRK1, J, JC, JR, M, MAXWRK, MINWRK;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP;
      COMPLEX            X
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CGGBAK, CGGBAL, CGGHRD, CHGEQZ, CLACPY, CLASCL, CLASET, CTGEVC, CTGSNA, CUNGQR, CUNMQR, SLASCL, XERBLA
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
      ILV = ILVL || ILVR

      NOSCL  = LSAME( BALANC, 'N' ) || LSAME( BALANC, 'P' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( !( NOSCL || LSAME( BALANC,'S' ) || LSAME( BALANC, 'B' ) ) ) {
         INFO = -1
      } else if ( IJOBVL <= 0 ) {
         INFO = -2
      } else if ( IJOBVR <= 0 ) {
         INFO = -3
      } else if ( !( WANTSN || WANTSE || WANTSB || WANTSV ) ) {
         INFO = -4
      } else if ( N < 0 ) {
         INFO = -5
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9
      } else if ( LDVL < 1 || ( ILVL && LDVL < N ) ) {
         INFO = -13
      } else if ( LDVR < 1 || ( ILVR && LDVR < N ) ) {
         INFO = -15
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV. The workspace is
        // computed assuming ILO = 1 and IHI = N, the worst case.)

      if ( INFO == 0 ) {
         if ( N == 0 ) {
            MINWRK = 1
            MAXWRK = 1
         } else {
            MINWRK = 2*N
            if ( WANTSE ) {
               MINWRK = 4*N
            } else if ( WANTSV || WANTSB ) {
               MINWRK = 2*N*( N + 1)
            }
            MAXWRK = MINWRK
            MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'CGEQRF', ' ', N, 1, N, 0 ) )             MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'CUNMQR', ' ', N, 1, N, 0 ) )
            if ( ILVL ) {
               MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'CUNGQR', ' ', N, 1, N, 0 ) )
            }
         }
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -25
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGGEVX', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = true;
      } else if ( ANRM > BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = true;
      }
      if (ILASCL) CALL CLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = false;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = true;
      } else if ( BNRM > BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = true;
      }
      if (ILBSCL) CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute and/or balance the matrix pair (A,B)
      // (Real Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)

      cggbal(BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, RWORK, IERR );

      // Compute ABNRM and BBNRM

      ABNRM = CLANGE( '1', N, N, A, LDA, RWORK( 1 ) )
      if ( ILASCL ) {
         RWORK( 1 ) = ABNRM
         slascl('G', 0, 0, ANRMTO, ANRM, 1, 1, RWORK( 1 ), 1, IERR );
         ABNRM = RWORK( 1 )
      }

      BBNRM = CLANGE( '1', N, N, B, LDB, RWORK( 1 ) )
      if ( ILBSCL ) {
         RWORK( 1 ) = BBNRM
         slascl('G', 0, 0, BNRMTO, BNRM, 1, 1, RWORK( 1 ), 1, IERR );
         BBNRM = RWORK( 1 )
      }

      // Reduce B to triangular form (QR decomposition of B)
      // (Complex Workspace: need N, prefer N*NB )

      IROWS = IHI + 1 - ILO
      if ( ILV || !WANTSN ) {
         ICOLS = N + 1 - ILO
      } else {
         ICOLS = IROWS
      }
      ITAU = 1
      IWRK = ITAU + IROWS
      cgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the unitary transformation to A
      // (Complex Workspace: need N, prefer N*NB)

      cunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL and/or VR
      // (Workspace: need N, prefer N*NB)

      if ( ILVL ) {
         claset('Full', N, N, CZERO, CONE, VL, LDVL );
         if ( IROWS > 1 ) {
            clacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         cungqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      if (ILVR) CALL CLASET( 'Full', N, N, CZERO, CONE, VR, LDVR );

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      if ( ILV || !WANTSN ) {

         // Eigenvectors requested -- work on whole matrix.

         cgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IERR );
      } else {
         cgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur forms and Schur vectors)
      // (Complex Workspace: need N)
      // (Real Workspace: need N)

      IWRK = ITAU
      if ( ILV || !WANTSN ) {
         CHTEMP = 'S'
      } else {
         CHTEMP = 'E'
      }

      chgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, RWORK, IERR );
      if ( IERR != 0 ) {
         if ( IERR > 0 && IERR <= N ) {
            INFO = IERR
         } else if ( IERR > N && IERR <= 2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 90
      }

      // Compute Eigenvectors and estimate condition numbers if desired
      // CTGEVC: (Complex Workspace: need 2*N )
              // (Real Workspace:    need 2*N )
      // CTGSNA: (Complex Workspace: need 2*N*N if SENSE='V' or 'B')
              // (Integer Workspace: need N+2 )

      if ( ILV || !WANTSN ) {
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

            ctgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), RWORK, IERR );
            if ( IERR != 0 ) {
               INFO = N + 2
               GO TO 90
            }
         }

         if ( !WANTSN ) {

            // compute eigenvectors (CTGEVC) and estimate condition
            // numbers (CTGSNA). Note that the definition of the condition
            // number is not invariant under transformation (u,v) to
            // (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
            // Schur form (S,T), Q and Z are orthogonal matrices. In order
            // to avoid using extra 2*N*N workspace, we have to
            // re-calculate eigenvectors and estimate the condition numbers
            // one at a time.

            for (I = 1; I <= N; I++) { // 20

               for (J = 1; J <= N; J++) { // 10
                  BWORK( J ) = false;
               } // 10
               BWORK( I ) = true;

               IWRK = N + 1
               IWRK1 = IWRK + N

               if ( WANTSE || WANTSB ) {
                  ctgevc('B', 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, 1, M, WORK( IWRK1 ), RWORK, IERR );
                  if ( IERR != 0 ) {
                     INFO = N + 2
                     GO TO 90
                  }
               }

               ctgsna(SENSE, 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, RCONDE( I ), RCONDV( I ), 1, M, WORK( IWRK1 ), LWORK-IWRK1+1, IWORK, IERR );

            } // 20
         }
      }

      // Undo balancing on VL and VR and normalization
      // (Workspace: none needed)

      if ( ILVL ) {
         cggbak(BALANC, 'L', N, ILO, IHI, LSCALE, RSCALE, N, VL, LDVL, IERR );

         for (JC = 1; JC <= N; JC++) { // 50
            TEMP = ZERO
            for (JR = 1; JR <= N; JR++) { // 30
               TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
            } // 30
            if (TEMP < SMLNUM) GO TO 50;
            TEMP = ONE / TEMP
            for (JR = 1; JR <= N; JR++) { // 40
               VL( JR, JC ) = VL( JR, JC )*TEMP
            } // 40
         } // 50
      }

      if ( ILVR ) {
         cggbak(BALANC, 'R', N, ILO, IHI, LSCALE, RSCALE, N, VR, LDVR, IERR );
         for (JC = 1; JC <= N; JC++) { // 80
            TEMP = ZERO
            for (JR = 1; JR <= N; JR++) { // 60
               TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
            } // 60
            if (TEMP < SMLNUM) GO TO 80;
            TEMP = ONE / TEMP
            for (JR = 1; JR <= N; JR++) { // 70
               VR( JR, JC ) = VR( JR, JC )*TEMP
            } // 70
         } // 80
      }

      // Undo scaling if necessary

      } // 90

      if (ILASCL) CALL CLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR );

      if (ILBSCL) CALL CLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      RETURN

      // End of CGGEVX

      }
