      void zgegs(JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILVSL, ILVSR, LQUERY;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, IRWORK, ITAU, IWORK, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SAFMIN, SMLNUM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD, ZHGEQZ, ZLACPY, ZLASCL, ZLASET, ZUNGQR, ZUNMQR
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ILAENV;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX
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

      // Test the input arguments

      LWKMIN = max( 2*N, 1 );
      LWKOPT = LWKMIN;
      WORK( 1 ) = LWKOPT;
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
         INFO = -11;
      } else if ( LDVSR < 1 || ( ILVSR && LDVSR < N ) ) {
         INFO = -13;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -15;
      }

      if ( INFO == 0 ) {
         NB1 = ILAENV( 1, 'ZGEQRF', ' ', N, N, -1, -1 );
         NB2 = ILAENV( 1, 'ZUNMQR', ' ', N, N, N, -1 );
         NB3 = ILAENV( 1, 'ZUNGQR', ' ', N, N, N, -1 );
         NB = max( NB1, NB2, NB3 );
         LOPT = N*( NB+1 );
         WORK( 1 ) = LOPT;
      }

      if ( INFO != 0 ) {
         xerbla('ZGEGS ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = DLAMCH( 'E' )*DLAMCH( 'B' );
      SAFMIN = DLAMCH( 'S' );
      SMLNUM = N*SAFMIN / EPS;
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

      if ( ILASCL ) {
         zlascl('G', -1, -1, ANRM, ANRMTO, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
      }

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

      if ( ILBSCL ) {
         zlascl('G', -1, -1, BNRM, BNRMTO, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
      }

      // Permute the matrix to make it more nearly triangular

      ILEFT = 1;
      IRIGHT = N + 1;
      IRWORK = IRIGHT + N;
      IWORK = 1;
      zggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWORK ), IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 1;
         GO TO 10;
      }

      // Reduce B to triangular form, and initialize VSL and/or VSR

      IROWS = IHI + 1 - ILO;
      ICOLS = N + 1 - ILO;
      ITAU = IWORK;
      IWORK = ITAU + IROWS;
      zgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO >= 0 ) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 2;
         GO TO 10;
      }

      zunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 3;
         GO TO 10;
      }

      if ( ILVSL ) {
         zlaset('Full', N, N, CZERO, CONE, VSL, LDVSL );
         zlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL );
         zungqr(IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO );
         if (IINFO >= 0) LWKOPT = max( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
         if ( IINFO != 0 ) {
            INFO = N + 4;
            GO TO 10;
         }
      }

      if (ILVSR) zlaset( 'Full', N, N, CZERO, CONE, VSR, LDVSR );

      // Reduce to generalized Hessenberg form

      zgghrd(JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 5;
         GO TO 10;
      }

      // Perform QZ algorithm, computing Schur vectors if desired

      IWORK = ITAU;
      zhgeqz('S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWORK ), LWORK+1-IWORK, RWORK( IRWORK ), IINFO );
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
         zggbak('P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSL, LDVSL, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 7;
            GO TO 10;
         }
      }
      if ( ILVSR ) {
         zggbak('P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSR, LDVSR, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 8;
            GO TO 10;
         }
      }

      // Undo scaling

      if ( ILASCL ) {
         zlascl('U', -1, -1, ANRMTO, ANRM, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
         zlascl('G', -1, -1, ANRMTO, ANRM, N, 1, ALPHA, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
      }

      if ( ILBSCL ) {
         zlascl('U', -1, -1, BNRMTO, BNRM, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
         zlascl('G', -1, -1, BNRMTO, BNRM, N, 1, BETA, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            return;
         }
      }

      } // 10
      WORK( 1 ) = LWKOPT;

      return;
      }
