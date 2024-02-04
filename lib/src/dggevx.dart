      void dggevx(BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, BWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      double             ABNRM, BBNRM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), LSCALE( * ), RCONDE( * ), RCONDV( * ), RSCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY, NOSCL, PAIR, WANTSB, WANTSE, WANTSN, WANTSV;
      String             CHTEMP;
      int                I, ICOLS, IERR, IJOBVL, IJOBVR, IN, IROWS, ITAU, IWRK, IWRK1, J, JC, JR, M, MAXWRK, MINWRK, MM;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP;
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, DLASCL, DLASET, DORGQR, DORMQR, DTGEVC, DTGSNA, XERBLA
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, DLANGE;
      // EXTERNAL lsame, ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

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

      NOSCL  = lsame( BALANC, 'N' ) || lsame( BALANC, 'P' );
      WANTSN = lsame( SENSE, 'N' );
      WANTSE = lsame( SENSE, 'E' );
      WANTSV = lsame( SENSE, 'V' );
      WANTSB = lsame( SENSE, 'B' );

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( !( lsame( BALANC, 'N' ) || lsame( BALANC, 'S' ) || lsame( BALANC, 'P' ) || lsame( BALANC, 'B' ) ) ) {
         INFO = -1;
      } else if ( IJOBVL <= 0 ) {
         INFO = -2;
      } else if ( IJOBVR <= 0 ) {
         INFO = -3;
      } else if ( !( WANTSN || WANTSE || WANTSB || WANTSV ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDVL < 1 || ( ILVL && LDVL < N ) ) {
         INFO = -14;
      } else if ( LDVR < 1 || ( ILVR && LDVR < N ) ) {
         INFO = -16;
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
            MINWRK = 1;
            MAXWRK = 1;
         } else {
            if ( NOSCL && !ILV ) {
               MINWRK = 2*N;
            } else {
               MINWRK = 6*N;
            }
            if ( WANTSE || WANTSB ) {
               MINWRK = 10*N;
            }
            if ( WANTSV || WANTSB ) {
               MINWRK = max( MINWRK, 2*N*( N + 4 ) + 16 );
            }
            MAXWRK = MINWRK;
            MAXWRK = max( MAXWRK, N + N*ILAENV( 1, 'DGEQRF', ' ', N, 1, N, 0 ) )             MAXWRK = max( MAXWRK, N + N*ILAENV( 1, 'DORMQR', ' ', N, 1, N, 0 ) );
            if ( ILVL ) {
               MAXWRK = max( MAXWRK, N + N*ILAENV( 1, 'DORGQR', ' ', N, 1, N, 0 ) );
            }
         }
         WORK[1] = MAXWRK;

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -26;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGGEVX', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;


      // Get machine constants

      EPS = DLAMCH( 'P' );
      SMLNUM = DLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = sqrt( SMLNUM ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', N, N, A, LDA, WORK );
      ILASCL = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ANRMTO = SMLNUM;
         ILASCL = true;
      } else if ( ANRM > BIGNUM ) {
         ANRMTO = BIGNUM;
         ILASCL = true;
      }
      if (ILASCL) dlascl( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = DLANGE( 'M', N, N, B, LDB, WORK );
      ILBSCL = false;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {
         BNRMTO = SMLNUM;
         ILBSCL = true;
      } else if ( BNRM > BIGNUM ) {
         BNRMTO = BIGNUM;
         ILBSCL = true;
      }
      if (ILBSCL) dlascl( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute and/or balance the matrix pair (A,B)
      // (Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)

      dggbal(BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, IERR );

      // Compute ABNRM and BBNRM

      ABNRM = DLANGE( '1', N, N, A, LDA, WORK( 1 ) );
      if ( ILASCL ) {
         WORK[1] = ABNRM;
         dlascl('G', 0, 0, ANRMTO, ANRM, 1, 1, WORK( 1 ), 1, IERR );
         ABNRM = WORK( 1 );
      }

      BBNRM = DLANGE( '1', N, N, B, LDB, WORK( 1 ) );
      if ( ILBSCL ) {
         WORK[1] = BBNRM;
         dlascl('G', 0, 0, BNRMTO, BNRM, 1, 1, WORK( 1 ), 1, IERR );
         BBNRM = WORK( 1 );
      }

      // Reduce B to triangular form (QR decomposition of B)
      // (Workspace: need N, prefer N*NB )

      IROWS = IHI + 1 - ILO;
      if ( ILV || !WANTSN ) {
         ICOLS = N + 1 - ILO;
      } else {
         ICOLS = IROWS;
      }
      ITAU = 1;
      IWRK = ITAU + IROWS;
      dgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to A
      // (Workspace: need N, prefer N*NB)

      dormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL and/or VR
      // (Workspace: need N, prefer N*NB)

      if ( ILVL ) {
         dlaset('Full', N, N, ZERO, ONE, VL, LDVL );
         if ( IROWS > 1 ) {
            dlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         dorgqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      if (ILVR) dlaset( 'Full', N, N, ZERO, ONE, VR, LDVR );

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      if ( ILV || !WANTSN ) {

         // Eigenvectors requested -- work on whole matrix.

         dgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IERR );
      } else {
         dgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur forms and Schur vectors)
      // (Workspace: need N)

      if ( ILV || !WANTSN ) {
         CHTEMP = 'S';
      } else {
         CHTEMP = 'E';
      }

      dhgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, IERR );
      if ( IERR != 0 ) {
         if ( IERR > 0 && IERR <= N ) {
            INFO = IERR;
         } else if ( IERR > N && IERR <= 2*N ) {
            INFO = IERR - N;
         } else {
            INFO = N + 1;
         }
         GO TO 130;
      }

      // Compute Eigenvectors and estimate condition numbers if desired
      // (Workspace: DTGEVC: need 6*N
                  // DTGSNA: need 2*N*(N+2)+16 if SENSE = 'V' or 'B',
                          // need N otherwise )

      if ( ILV || !WANTSN ) {
         if ( ILV ) {
            if ( ILVL ) {
               if ( ILVR ) {
                  CHTEMP = 'B';
               } else {
                  CHTEMP = 'L';
               }
            } else {
               CHTEMP = 'R';
            }

            dtgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK, IERR );
            if ( IERR != 0 ) {
               INFO = N + 2;
               GO TO 130;
            }
         }

         if ( !WANTSN ) {

            // compute eigenvectors (DTGEVC) and estimate condition
            // numbers (DTGSNA). Note that the definition of the condition
            // number is not invariant under transformation (u,v) to
            // (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
            // Schur form (S,T), Q and Z are orthogonal matrices. In order
            // to avoid using extra 2*N*N workspace, we have to recalculate
            // eigenvectors and estimate one condition numbers at a time.

            PAIR = false;
            for (I = 1; I <= N; I++) { // 20

               if ( PAIR ) {
                  PAIR = false;
                  GO TO 20;
               }
               MM = 1;
               if ( I < N ) {
                  if ( A( I+1, I ) != ZERO ) {
                     PAIR = true;
                     MM = 2;
                  }
               }

               for (J = 1; J <= N; J++) { // 10
                  BWORK[J] = false;
               } // 10
               if ( MM == 1 ) {
                  BWORK[I] = true;
               } else if ( MM == 2 ) {
                  BWORK[I] = true;
                  BWORK[I+1] = true;
               }

               IWRK = MM*N + 1;
               IWRK1 = IWRK + MM*N;

               // Compute a pair of left and right eigenvectors.
               // (compute workspace: need up to 4*N + 6*N)

               if ( WANTSE || WANTSB ) {
                  dtgevc('B', 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, MM, M, WORK( IWRK1 ), IERR );
                  if ( IERR != 0 ) {
                     INFO = N + 2;
                     GO TO 130;
                  }
               }

               dtgsna(SENSE, 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, RCONDE( I ), RCONDV( I ), MM, M, WORK( IWRK1 ), LWORK-IWRK1+1, IWORK, IERR );

            } // 20
         }
      }

      // Undo balancing on VL and VR and normalization
      // (Workspace: none needed)

      if ( ILVL ) {
         dggbak(BALANC, 'L', N, ILO, IHI, LSCALE, RSCALE, N, VL, LDVL, IERR );

         for (JC = 1; JC <= N; JC++) { // 70
            if( ALPHAI( JC ) < ZERO ) GO TO 70;
            TEMP = ZERO;
            if ( ALPHAI( JC ) == ZERO ) {
               for (JR = 1; JR <= N; JR++) { // 30
                  TEMP = max( TEMP, ( VL( JR, JC ) ) ).abs();
               } // 30
            } else {
               for (JR = 1; JR <= N; JR++) { // 40
                  TEMP = max( TEMP, ( VL( JR, JC ) ).abs()+ ( VL( JR, JC+1 ) ) ).abs();
               } // 40
            }
            if (TEMP < SMLNUM) GO TO 70;
            TEMP = ONE / TEMP;
            if ( ALPHAI( JC ) == ZERO ) {
               for (JR = 1; JR <= N; JR++) { // 50
                  VL[JR, JC] = VL( JR, JC )*TEMP;
               } // 50
            } else {
               for (JR = 1; JR <= N; JR++) { // 60
                  VL[JR, JC] = VL( JR, JC )*TEMP;
                  VL[JR, JC+1] = VL( JR, JC+1 )*TEMP;
               } // 60
            }
         } // 70
      }
      if ( ILVR ) {
         dggbak(BALANC, 'R', N, ILO, IHI, LSCALE, RSCALE, N, VR, LDVR, IERR );
         for (JC = 1; JC <= N; JC++) { // 120
            if( ALPHAI( JC ) < ZERO ) GO TO 120;
            TEMP = ZERO;
            if ( ALPHAI( JC ) == ZERO ) {
               for (JR = 1; JR <= N; JR++) { // 80
                  TEMP = max( TEMP, ( VR( JR, JC ) ) ).abs();
               } // 80
            } else {
               for (JR = 1; JR <= N; JR++) { // 90
                  TEMP = max( TEMP, ( VR( JR, JC ) ).abs()+ ( VR( JR, JC+1 ) ) ).abs();
               } // 90
            }
            if (TEMP < SMLNUM) GO TO 120;
            TEMP = ONE / TEMP;
            if ( ALPHAI( JC ) == ZERO ) {
               for (JR = 1; JR <= N; JR++) { // 100
                  VR[JR, JC] = VR( JR, JC )*TEMP;
               } // 100
            } else {
               for (JR = 1; JR <= N; JR++) { // 110
                  VR[JR, JC] = VR( JR, JC )*TEMP;
                  VR[JR, JC+1] = VR( JR, JC+1 )*TEMP;
               } // 110
            }
         } // 120
      }

      // Undo scaling if necessary

      } // 130

      if ( ILASCL ) {
         dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
         dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
      }

      if ( ILBSCL ) {
         dlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      WORK[1] = MAXWRK;
      return;
      }
