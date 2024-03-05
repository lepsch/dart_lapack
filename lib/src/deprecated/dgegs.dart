      void dgegs(final int JOBVSL, final int JOBVSR, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int ALPHAR, final int ALPHAI, final int BETA, final Matrix<double> VSL_, final int LDVSL, final Matrix<double> VSR_, final int LDVSR, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.having();
  final B = B_.having();
  final VSL = VSL_.having();
  final VSR = VSR_.having();
  final WORK = WORK_.having();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBVSL, JOBVSR;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N;
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               ILASCL, ILBSCL, ILVSL, ILVSR, LQUERY;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, ITAU, IWORK, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SAFMIN, SMLNUM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, DLASCL, DLASET, DORGQR, DORMQR, XERBLA
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, DLANGE;
      // EXTERNAL lsame, ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX

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

      // Test the input arguments

      LWKMIN = max( 4*N, 1 );
      LWKOPT = LWKMIN;
      WORK[1] = LWKOPT;
      LQUERY = ( LWORK == -1 );
      INFO = 0;
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
      } else if ( LDVSL < 1 || ( ILVSL && LDVSL < N ) ) {
         INFO = -12;
      } else if ( LDVSR < 1 || ( ILVSR && LDVSR < N ) ) {
         INFO = -14;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -16;
      }

      if ( INFO == 0 ) {
         NB1 = ilaenv( 1, 'DGEQRF', ' ', N, N, -1, -1 );
         NB2 = ilaenv( 1, 'DORMQR', ' ', N, N, N, -1 );
         NB3 = ilaenv( 1, 'DORGQR', ' ', N, N, N, -1 );
         NB = max( NB1, NB2, NB3 );
         LOPT = 2*N + N*( NB+1 );
         WORK[1] = LOPT;
      }

      if ( INFO != 0 ) {
         xerbla('DGEGS ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = dlamch( 'E' )*dlamch( 'B' );
      SAFMIN = dlamch( 'S' );
      SMLNUM = N*SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = dlange( 'M', N, N, A, LDA, WORK );
      ILASCL = false;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ANRMTO = SMLNUM;
         ILASCL = true;
      } else if ( ANRM > BIGNUM ) {
         ANRMTO = BIGNUM;
         ILASCL = true;
      }

      if ( ILASCL ) {
         dlascl('G', -1, -1, ANRM, ANRMTO, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
      }

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = dlange( 'M', N, N, B, LDB, WORK );
      ILBSCL = false;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {
         BNRMTO = SMLNUM;
         ILBSCL = true;
      } else if ( BNRM > BIGNUM ) {
         BNRMTO = BIGNUM;
         ILBSCL = true;
      }

      if ( ILBSCL ) {
         dlascl('G', -1, -1, BNRM, BNRMTO, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
      }

      // Permute the matrix to make it more nearly triangular
      // Workspace layout:  (2*N words -- "work..." not actually used)
      //    left_permutation, right_permutation, work...

      ILEFT = 1;
      IRIGHT = N + 1;
      IWORK = IRIGHT + N;
      dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWORK ), IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 1;
         GO TO 10;
      }

      // Reduce B to triangular form, and initialize VSL and/or VSR
      // Workspace layout:  ("work..." must have at least N words)
      //    left_permutation, right_permutation, tau, work...

      IROWS = IHI + 1 - ILO;
      ICOLS = N + 1 - ILO;
      ITAU = IWORK;
      IWORK = ITAU + IROWS;
      dgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO >= 0 ) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 2;
         GO TO 10;
      }

      dormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 3;
         GO TO 10;
      }

      if ( ILVSL ) {
         dlaset('Full', N, N, ZERO, ONE, VSL, LDVSL );
         dlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL );
         dorgqr(IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO );
         if (IINFO >= 0) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
         if ( IINFO != 0 ) {
            INFO = N + 4;
            GO TO 10;
         }
      }

      if (ILVSR) dlaset( 'Full', N, N, ZERO, ONE, VSR, LDVSR );

      // Reduce to generalized Hessenberg form

      dgghrd(JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 5;
         GO TO 10;
      }

      // Perform QZ algorithm, computing Schur vectors if desired
      // Workspace layout:  ("work..." must have at least 1 word)
      //    left_permutation, right_permutation, work...

      IWORK = ITAU;
      dhgeqz('S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         if ( IINFO > 0 && IINFO <= N ) {
            INFO = IINFO;
         } else if ( IINFO > N && IINFO <= 2*N ) {
            INFO = IINFO - N;
         } else {
            INFO = N + 6;
         }
         GO TO 10;
      }

      // Apply permutation to VSL and VSR

      if ( ILVSL ) {
         dggbak('P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSL, LDVSL, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 7;
            GO TO 10;
         }
      }
      if ( ILVSR ) {
         dggbak('P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSR, LDVSR, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 8;
            GO TO 10;
         }
      }

      // Undo scaling

      if ( ILASCL ) {
         dlascl('H', -1, -1, ANRMTO, ANRM, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
         dlascl('G', -1, -1, ANRMTO, ANRM, N, 1, ALPHAR, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
         dlascl('G', -1, -1, ANRMTO, ANRM, N, 1, ALPHAI, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
      }

      if ( ILBSCL ) {
         dlascl('U', -1, -1, BNRMTO, BNRM, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
         dlascl('G', -1, -1, BNRMTO, BNRM, N, 1, BETA, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
      }

      } // 10
      WORK[1] = LWKOPT;

      }
