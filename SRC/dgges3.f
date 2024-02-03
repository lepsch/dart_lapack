      SUBROUTINE DGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, BWORK, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR, SORT;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * );
      // ..
      // .. Function Arguments ..
      bool               SELCTG;
      // EXTERNAL SELCTG
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, LST2SL, WANTST;
      int                I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IP, IRIGHT, IROWS, ITAU, IWRK, LWKOPT, LWKMIN;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL, PVSR, SAFMAX, SAFMIN, SMLNUM;
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      double             DIF( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DGGBAK, DGGBAL, DGGHD3, DLAQZ0, DLACPY, DLASCL, DLASET, DORGQR, DORMQR, DTGSEN, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

      if ( LSAME( JOBVSL, 'N' ) ) {
         IJOBVL = 1;
         ILVSL = false;
      } else if ( LSAME( JOBVSL, 'V' ) ) {
         IJOBVL = 2;
         ILVSL = true;
      } else {
         IJOBVL = -1;
         ILVSL = false;
      }

      if ( LSAME( JOBVSR, 'N' ) ) {
         IJOBVR = 1;
         ILVSR = false;
      } else if ( LSAME( JOBVSR, 'V' ) ) {
         IJOBVR = 2;
         ILVSR = true;
      } else {
         IJOBVR = -1;
         ILVSR = false;
      }

      WANTST = LSAME( SORT, 'S' );

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( N == 0 ) {
         LWKMIN = 1;
      } else {
         LWKMIN = 6*N+16;
      }

      if ( IJOBVL <= 0 ) {
         INFO = -1;
      } else if ( IJOBVR <= 0 ) {
         INFO = -2;
      } else if ( ( !WANTST ) && ( !LSAME( SORT, 'N' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDVSL < 1 || ( ILVSL && LDVSL < N ) ) {
         INFO = -15;
      } else if ( LDVSR < 1 || ( ILVSR && LDVSR < N ) ) {
         INFO = -17;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -19;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         dgeqrf(N, N, B, LDB, WORK, WORK, -1, IERR );
         LWKOPT = max( LWKMIN, 3*N+INT( WORK( 1 ) ) );
         dormqr('L', 'T', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR );
         LWKOPT = max( LWKOPT, 3*N+INT( WORK( 1 ) ) );
         if ( ILVSL ) {
            dorgqr(N, N, N, VSL, LDVSL, WORK, WORK, -1, IERR );
            LWKOPT = max( LWKOPT, 3*N+INT( WORK( 1 ) ) );
         }
         dgghd3(JOBVSL, JOBVSR, N, 1, N, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, WORK, -1, IERR );
         LWKOPT = max( LWKOPT, 3*N+INT( WORK( 1 ) ) );
         dlaqz0('S', JOBVSL, JOBVSR, N, 1, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, -1, 0, IERR );
         LWKOPT = max( LWKOPT, 2*N+INT( WORK( 1 ) ) );
         if ( WANTST ) {
            dtgsen(0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK, -1, IDUM, 1, IERR );
            LWKOPT = max( LWKOPT, 2*N+INT( WORK( 1 ) ) );
         }
         if ( N == 0 ) {
            WORK( 1 ) = 1;
         } else {
            WORK( 1 ) = LWKOPT;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGGES3 ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         SDIM = 0;
         return;
      }

      // Get machine constants

      EPS = DLAMCH( 'P' );
      SAFMIN = DLAMCH( 'S' );
      SAFMAX = ONE / SAFMIN;
      SMLNUM = sqrt( SAFMIN ) / EPS;
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
      if (ILASCL) CALL DLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

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
      if (ILBSCL) CALL DLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute the matrix to make it more nearly triangular

      ILEFT = 1;
      IRIGHT = N + 1;
      IWRK = IRIGHT + N;
      dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)

      IROWS = IHI + 1 - ILO;
      ICOLS = N + 1 - ILO;
      ITAU = IWRK;
      IWRK = ITAU + IROWS;
      dgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to matrix A

      dormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VSL

      if ( ILVSL ) {
         dlaset('Full', N, N, ZERO, ONE, VSL, LDVSL );
         if ( IROWS > 1 ) {
            dlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL );
         }
         dorgqr(IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VSR

      if (ILVSR) CALL DLASET( 'Full', N, N, ZERO, ONE, VSR, LDVSR );

      // Reduce to generalized Hessenberg form

      dgghd3(JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Perform QZ algorithm, computing Schur vectors if desired

      IWRK = ITAU;
      dlaqz0('S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, 0, IERR );
      if ( IERR != 0 ) {
         if ( IERR > 0 && IERR <= N ) {
            INFO = IERR;
         } else if ( IERR > N && IERR <= 2*N ) {
            INFO = IERR - N;
         } else {
            INFO = N + 1;
         }
         GO TO 50;
      }

      // Sort eigenvalues ALPHA/BETA if desired

      SDIM = 0;
      if ( WANTST ) {

         // Undo scaling on eigenvalues before SELCTGing

         if ( ILASCL ) {
            dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
            dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
         }
         if (ILBSCL) CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );

         // Select eigenvalues

         for (I = 1; I <= N; I++) { // 10
            BWORK( I ) = SELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) );
         } // 10

         dtgsen(0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, IERR );
         if (IERR == 1) INFO = N + 3;

      }

      // Apply back-permutation to VSL and VSR

      if (ILVSL) CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSL, LDVSL, IERR );

      if (ILVSR) CALL DGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSR, LDVSR, IERR );

      // Check if unscaling would cause over/underflow, if so, rescale
      // (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
      // B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)

      if ( ILASCL ) {
         for (I = 1; I <= N; I++) { // 20
            if ( ALPHAI( I ) != ZERO ) {
               if ( ( ALPHAR( I ) / SAFMAX ) > ( ANRMTO / ANRM ) || ( SAFMIN / ALPHAR( I ) ) > ( ANRM / ANRMTO ) ) {
                  WORK( 1 ) = ABS( A( I, I ) / ALPHAR( I ) );
                  BETA( I ) = BETA( I )*WORK( 1 );
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 );
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 );
               } else if ( ( ALPHAI( I ) / SAFMAX ) > ( ANRMTO / ANRM ) || ( SAFMIN / ALPHAI( I ) ) > ( ANRM / ANRMTO ) ) {
                  WORK( 1 ) = ABS( A( I, I+1 ) / ALPHAI( I ) );
                  BETA( I ) = BETA( I )*WORK( 1 );
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 );
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 );
               }
            }
         } // 20
      }

      if ( ILBSCL ) {
         for (I = 1; I <= N; I++) { // 30
            if ( ALPHAI( I ) != ZERO ) {
               if ( ( BETA( I ) / SAFMAX ) > ( BNRMTO / BNRM ) || ( SAFMIN / BETA( I ) ) > ( BNRM / BNRMTO ) ) {
                  WORK( 1 ) = ABS( B( I, I ) / BETA( I ) );
                  BETA( I ) = BETA( I )*WORK( 1 );
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 );
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 );
               }
            }
         } // 30
      }

      // Undo scaling

      if ( ILASCL ) {
         dlascl('H', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR );
         dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
         dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
      }

      if ( ILBSCL ) {
         dlascl('U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR );
         dlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      if ( WANTST ) {

         // Check if reordering is correct

         LASTSL = true;
         LST2SL = true;
         SDIM = 0;
         IP = 0;
         for (I = 1; I <= N; I++) { // 40
            CURSL = SELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) );
            if ( ALPHAI( I ) == ZERO ) {
               if (CURSL) SDIM = SDIM + 1;
               IP = 0;
               if (CURSL && !LASTSL) INFO = N + 2;
            } else {
               if ( IP == 1 ) {

                  // Last eigenvalue of conjugate pair

                  CURSL = CURSL || LASTSL;
                  LASTSL = CURSL;
                  if (CURSL) SDIM = SDIM + 2;
                  IP = -1;
                  if (CURSL && !LST2SL) INFO = N + 2;
               } else {

                  // First eigenvalue of conjugate pair

                  IP = 1;
               }
            }
            LST2SL = LASTSL;
            LASTSL = CURSL;
         } // 40

      }

      } // 50

      WORK( 1 ) = LWKOPT;

      return;

      // End of DGGES3

      }
