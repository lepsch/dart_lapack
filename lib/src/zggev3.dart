      void zggev3(JOBVL, JOBVR, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, ALPHA, BETA, final Matrix<double> VL, final int LDVL, final Matrix<double> VR, final int LDVR, final Array<double> WORK, final int LWORK, RWORK, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      double             RWORK( * );
      Complex         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, IRWRK, ITAU, IWRK, JC, JR, LWKMIN, LWKOPT;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP;
      Complex         X;
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHD3, ZLAQZ0, ZLACPY, ZLASCL, ZLASET, ZTGEVC, ZUNGQR, ZUNMQR
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL lsame, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, SQRT
      // ..
      // .. Statement Functions ..
      double             ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1[X] = ( X.toDouble() ).abs() + ( DIMAG( X ) ).abs();

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

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      LWKMIN = max( 1, 2*N );
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
         INFO = -11;
      } else if ( LDVR < 1 || ( ILVR && LDVR < N ) ) {
         INFO = -13;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -15;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         zgeqrf(N, N, B, LDB, WORK, WORK, -1, IERR );
         LWKOPT = max( LWKMIN, N+INT( WORK( 1 ) ) );
         zunmqr('L', 'C', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR );
         LWKOPT = max( LWKOPT, N+INT( WORK( 1 ) ) );
         if ( ILVL ) {
            zungqr(N, N, N, VL, LDVL, WORK, WORK, -1, IERR );
            LWKOPT = max( LWKOPT, N+INT( WORK( 1 ) ) );
         }
         if ( ILV ) {
            zgghd3(JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR );
            LWKOPT = max( LWKOPT, N+INT( WORK( 1 ) ) );
            zlaqz0('S', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, -1, RWORK, 0, IERR );
            LWKOPT = max( LWKOPT, N+INT( WORK( 1 ) ) );
         } else {
            zgghd3(JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR );
            LWKOPT = max( LWKOPT, N+INT( WORK( 1 ) ) );
            zlaqz0('E', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, -1, RWORK, 0, IERR );
            LWKOPT = max( LWKOPT, N+INT( WORK( 1 ) ) );
         }
         if ( N == 0 ) {
            WORK[1] = 1;
         } else {
            WORK[1] = DCMPLX( LWKOPT );
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZGGEV3 ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = dlamch( 'E' )*dlamch( 'B' );
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

      // Permute the matrices A, B to isolate eigenvalues if possible

      ILEFT = 1;
      IRIGHT = N + 1;
      IRWRK = IRIGHT + N;
      zggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)

      IROWS = IHI + 1 - ILO;
      if ( ILV ) {
         ICOLS = N + 1 - ILO;
      } else {
         ICOLS = IROWS;
      }
      ITAU = 1;
      IWRK = ITAU + IROWS;
      zgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to matrix A

      zunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL

      if ( ILVL ) {
         zlaset('Full', N, N, CZERO, CONE, VL, LDVL );
         if ( IROWS > 1 ) {
            zlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         zungqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VR

      if (ILVR) zlaset( 'Full', N, N, CZERO, CONE, VR, LDVR );

      // Reduce to generalized Hessenberg form

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         zgghd3(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR );
      } else {
         zgghd3('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur form and Schur vectors)

      IWRK = ITAU;
      if ( ILV ) {
         CHTEMP = 'S';
      } else {
         CHTEMP = 'E';
      }
      zlaqz0(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, RWORK( IRWRK ), 0, IERR );
      if ( IERR != 0 ) {
         if ( IERR > 0 && IERR <= N ) {
            INFO = IERR;
         } else if ( IERR > N && IERR <= 2*N ) {
            INFO = IERR - N;
         } else {
            INFO = N + 1;
         }
         GO TO 70;
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

         ztgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), RWORK( IRWRK ), IERR );
         if ( IERR != 0 ) {
            INFO = N + 2;
            GO TO 70;
         }

         // Undo balancing on VL and VR and normalization

         if ( ILVL ) {
            zggbak('P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VL, LDVL, IERR );
            for (JC = 1; JC <= N; JC++) { // 30
               TEMP = ZERO;
               for (JR = 1; JR <= N; JR++) { // 10
                  TEMP = max( TEMP, ABS1( VL( JR, JC ) ) );
               } // 10
               if (TEMP < SMLNUM) GO TO 30;
               TEMP = ONE / TEMP;
               for (JR = 1; JR <= N; JR++) { // 20
                  VL[JR][JC] = VL( JR, JC )*TEMP;
               } // 20
            } // 30
         }
         if ( ILVR ) {
            zggbak('P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VR, LDVR, IERR );
            for (JC = 1; JC <= N; JC++) { // 60
               TEMP = ZERO;
               for (JR = 1; JR <= N; JR++) { // 40
                  TEMP = max( TEMP, ABS1( VR( JR, JC ) ) );
               } // 40
               if (TEMP < SMLNUM) GO TO 60;
               TEMP = ONE / TEMP;
               for (JR = 1; JR <= N; JR++) { // 50
                  VR[JR][JC] = VR( JR, JC )*TEMP;
               } // 50
            } // 60
         }
      }

      // Undo scaling if necessary

      } // 70

      if (ILASCL) zlascl( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR );

      if (ILBSCL) zlascl( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );

      WORK[1] = DCMPLX( LWKOPT );
      }
