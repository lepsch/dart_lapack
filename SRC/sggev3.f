      void sggev3(JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, ITAU, IWRK, JC, JR, LWKOPT, LWKMIN;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP;
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRF, SGGBAK, SGGBAL, SGGHD3, SLAQZ0, SLACPY, SLASCL, SLASET, SORGQR, SORMQR, STGEVC
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- REAL               SLAMCH, SLANGE, SROUNDUP_LWORK;
      // EXTERNAL LSAME, SLAMCH, SLANGE, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
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

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      LWKMIN = max( 1, 8*N );
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

      // Compute workspace

      if ( INFO == 0 ) {
         sgeqrf(N, N, B, LDB, WORK, WORK, -1, IERR );
         LWKOPT = max( LWKMIN, 3*N+INT( WORK( 1 ) ) );
         sormqr('L', 'T', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR );
         LWKOPT = max( LWKOPT, 3*N+INT( WORK( 1 ) ) );
         sgghd3(JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR );
         LWKOPT = max( LWKOPT, 3*N+INT( WORK( 1 ) ) );
         if ( ILVL ) {
            sorgqr(N, N, N, VL, LDVL, WORK, WORK, -1, IERR );
            LWKOPT = max( LWKOPT, 3*N+INT( WORK( 1 ) ) );
            slaqz0('S', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, -1, 0, IERR );
            LWKOPT = max( LWKOPT, 2*N+INT( WORK( 1 ) ) );
         } else {
            slaqz0('E', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, -1, 0, IERR );
            LWKOPT = max( LWKOPT, 2*N+INT( WORK( 1 ) ) );
         }
         if ( N == 0 ) {
            WORK( 1 ) = 1;
         } else {
            WORK( 1 ) = SROUNDUP_LWORK( LWKOPT );
         }
      }

      if ( INFO != 0 ) {
         xerbla('SGGEV3 ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = sqrt( SMLNUM ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', N, N, A, LDA, WORK );
      ILASCL = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ANRMTO = SMLNUM;
         ILASCL = true;
      } else if ( ANRM > BIGNUM ) {
         ANRMTO = BIGNUM;
         ILASCL = true;
      }
      if (ILASCL) slascl( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = SLANGE( 'M', N, N, B, LDB, WORK );
      ILBSCL = false;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {
         BNRMTO = SMLNUM;
         ILBSCL = true;
      } else if ( BNRM > BIGNUM ) {
         BNRMTO = BIGNUM;
         ILBSCL = true;
      }
      if (ILBSCL) slascl( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute the matrices A, B to isolate eigenvalues if possible

      ILEFT = 1;
      IRIGHT = N + 1;
      IWRK = IRIGHT + N;
      sggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)

      IROWS = IHI + 1 - ILO;
      if ( ILV ) {
         ICOLS = N + 1 - ILO;
      } else {
         ICOLS = IROWS;
      }
      ITAU = IWRK;
      IWRK = ITAU + IROWS;
      sgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to matrix A

      sormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL

      if ( ILVL ) {
         slaset('Full', N, N, ZERO, ONE, VL, LDVL );
         if ( IROWS > 1 ) {
            slacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         sorgqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VR

      if (ILVR) slaset( 'Full', N, N, ZERO, ONE, VR, LDVR );

      // Reduce to generalized Hessenberg form

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         sgghd3(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR );
      } else {
         sgghd3('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur forms and Schur vectors)

      IWRK = ITAU;
      if ( ILV ) {
         CHTEMP = 'S';
      } else {
         CHTEMP = 'E';
      }
      slaqz0(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, 0, IERR );
      if ( IERR != 0 ) {
         if ( IERR > 0 && IERR <= N ) {
            INFO = IERR;
         } else if ( IERR > N && IERR <= 2*N ) {
            INFO = IERR - N;
         } else {
            INFO = N + 1;
         }
         GO TO 110;
      }

      // Compute Eigenvectors

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
         stgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), IERR );
         if ( IERR != 0 ) {
            INFO = N + 2;
            GO TO 110;
         }

         // Undo balancing on VL and VR and normalization

         if ( ILVL ) {
            sggbak('P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VL, LDVL, IERR );
            for (JC = 1; JC <= N; JC++) { // 50
               if( ALPHAI( JC ) < ZERO ) GO TO 50;
               TEMP = ZERO;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 10
                     TEMP = max( TEMP, ( VL( JR, JC ) ) ).abs();
                  } // 10
               } else {
                  for (JR = 1; JR <= N; JR++) { // 20
                     TEMP = max( TEMP, ( VL( JR, JC ) ).abs()+ ( VL( JR, JC+1 ) ) ).abs();
                  } // 20
               }
               if (TEMP < SMLNUM) GO TO 50;
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
            sggbak('P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VR, LDVR, IERR );
            for (JC = 1; JC <= N; JC++) { // 100
               if( ALPHAI( JC ) < ZERO ) GO TO 100;
               TEMP = ZERO;
               if ( ALPHAI( JC ) == ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 60
                     TEMP = max( TEMP, ( VR( JR, JC ) ) ).abs();
                  } // 60
               } else {
                  for (JR = 1; JR <= N; JR++) { // 70
                     TEMP = max( TEMP, ( VR( JR, JC ) ).abs()+ ( VR( JR, JC+1 ) ) ).abs();
                  } // 70
               }
               if (TEMP < SMLNUM) GO TO 100;
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

      // Undo scaling if necessary

      } // 110

      if ( ILASCL ) {
         slascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
         slascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
      }

      if ( ILBSCL ) {
         slascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT );
      return;
      }
