      void zgges(JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, BWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBVSL, JOBVSR, SORT;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM;
      bool               BWORK( * );
      double             RWORK( * );
      Complex         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * );
      // ..
      // .. Function Arguments ..
      bool               SELCTG;
      // EXTERNAL SELCTG
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, WANTST;
      int                I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, IRWRK, ITAU, IWRK, LWKMIN, LWKOPT;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL, PVSR, SMLNUM;
      int                IDUM( 1 );
      double             DIF( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD, ZHGEQZ, ZLACPY, ZLASCL, ZLASET, ZTGSEN, ZUNGQR, ZUNMQR
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL lsame, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Decode the input arguments

      if ( lsame( JOBVSL, 'N' ) ) {
         IJOBVL = 1;
         ILVSL = false;
      } else if ( lsame( JOBVSL, 'V' ) ) {
         IJOBVL = 2;
         ILVSL = true;
      } else {
         IJOBVL = -1;
         ILVSL = false;
      }

      if ( lsame( JOBVSR, 'N' ) ) {
         IJOBVR = 1;
         ILVSR = false;
      } else if ( lsame( JOBVSR, 'V' ) ) {
         IJOBVR = 2;
         ILVSR = true;
      } else {
         IJOBVR = -1;
         ILVSR = false;
      }

      WANTST = lsame( SORT, 'S' );

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( IJOBVL <= 0 ) {
         INFO = -1;
      } else if ( IJOBVR <= 0 ) {
         INFO = -2;
      } else if ( ( !WANTST ) && ( !lsame( SORT, 'N' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDVSL < 1 || ( ILVSL && LDVSL < N ) ) {
         INFO = -14;
      } else if ( LDVSR < 1 || ( ILVSR && LDVSR < N ) ) {
         INFO = -16;
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      if ( INFO == 0 ) {
         LWKMIN = max( 1, 2*N );
         LWKOPT = max( 1, N + N*ilaenv( 1, 'ZGEQRF', ' ', N, 1, N, 0 ) );
         LWKOPT = max( LWKOPT, N + N*ilaenv( 1, 'ZUNMQR', ' ', N, 1, N, -1 ) );
         if ( ILVSL ) {
            LWKOPT = max( LWKOPT, N + N*ilaenv( 1, 'ZUNGQR', ' ', N, 1, N, -1 ) );
         }
         WORK[1] = LWKOPT;

         if (LWORK < LWKMIN && !LQUERY) INFO = -18;
      }

      if ( INFO != 0 ) {
         xerbla('ZGGES ', -INFO );
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

      EPS = dlamch( 'P' );
      SMLNUM = dlamch( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = sqrt( SMLNUM ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK );
      ILASCL = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ANRMTO = SMLNUM;
         ILASCL = true;
      } else if ( ANRM > BIGNUM ) {
         ANRMTO = BIGNUM;
         ILASCL = true;
      }

      if (ILASCL) zlascl( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK );
      ILBSCL = false;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {
         BNRMTO = SMLNUM;
         ILBSCL = true;
      } else if ( BNRM > BIGNUM ) {
         BNRMTO = BIGNUM;
         ILBSCL = true;
      }

      if (ILBSCL) zlascl( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute the matrix to make it more nearly triangular
      // (Real Workspace: need 6*N)

      ILEFT = 1;
      IRIGHT = N + 1;
      IRWRK = IRIGHT + N;
      zggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)
      // (Complex Workspace: need N, prefer N*NB)

      IROWS = IHI + 1 - ILO;
      ICOLS = N + 1 - ILO;
      ITAU = 1;
      IWRK = ITAU + IROWS;
      zgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to matrix A
      // (Complex Workspace: need N, prefer N*NB)

      zunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VSL
      // (Complex Workspace: need N, prefer N*NB)

      if ( ILVSL ) {
         zlaset('Full', N, N, CZERO, CONE, VSL, LDVSL );
         if ( IROWS > 1 ) {
            zlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL );
         }
         zungqr(IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VSR

      if (ILVSR) zlaset( 'Full', N, N, CZERO, CONE, VSR, LDVSR );

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      zgghrd(JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IERR );

      SDIM = 0;

      // Perform QZ algorithm, computing Schur vectors if desired
      // (Complex Workspace: need N)
      // (Real Workspace: need N)

      IWRK = ITAU;
      zhgeqz('S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, RWORK( IRWRK ), IERR );
      if ( IERR != 0 ) {
         if ( IERR > 0 && IERR <= N ) {
            INFO = IERR;
         } else if ( IERR > N && IERR <= 2*N ) {
            INFO = IERR - N;
         } else {
            INFO = N + 1;
         }
         GO TO 30;
      }

      // Sort eigenvalues ALPHA/BETA if desired
      // (Workspace: none needed)

      if ( WANTST ) {

         // Undo scaling on eigenvalues before selecting

         if (ILASCL) zlascl( 'G', 0, 0, ANRM, ANRMTO, N, 1, ALPHA, N, IERR );
         IF( ILBSCL ) zlascl( 'G', 0, 0, BNRM, BNRMTO, N, 1, BETA, N, IERR );

         // Select eigenvalues

         for (I = 1; I <= N; I++) { // 10
            BWORK[I] = SELCTG( ALPHA( I ), BETA( I ) );
         } // 10

         ztgsen(0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, IERR );
         if (IERR == 1) INFO = N + 3;

      }

      // Apply back-permutation to VSL and VSR
      // (Workspace: none needed)

      if (ILVSL) zggbak( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSL, LDVSL, IERR );
      IF( ILVSR ) zggbak( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSR, LDVSR, IERR );

      // Undo scaling

      if ( ILASCL ) {
         zlascl('U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR );
         zlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR );
      }

      if ( ILBSCL ) {
         zlascl('U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR );
         zlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      if ( WANTST ) {

         // Check if reordering is correct

         LASTSL = true;
         SDIM = 0;
         for (I = 1; I <= N; I++) { // 20
            CURSL = SELCTG( ALPHA( I ), BETA( I ) );
            if (CURSL) SDIM = SDIM + 1;
            IF( CURSL && !LASTSL ) INFO = N + 2;
            LASTSL = CURSL;
         } // 20

      }

      } // 30

      WORK[1] = LWKOPT;

      return;
      }
