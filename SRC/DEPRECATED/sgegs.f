      SUBROUTINE SGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, INFO );

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILVSL, ILVSR, LQUERY;
      int                ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, ITAU, IWORK, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SAFMIN, SMLNUM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRF, SGGBAK, SGGBAL, SGGHRD, SHGEQZ, SLACPY, SLASCL, SLASET, SORGQR, SORMQR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANGE;
      // EXTERNAL ILAENV, LSAME, SLAMCH, SLANGE
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

      LWKMIN = MAX( 4*N, 1 );
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
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7;
      } else if ( LDVSL < 1 || ( ILVSL && LDVSL < N ) ) {
         INFO = -12;
      } else if ( LDVSR < 1 || ( ILVSR && LDVSR < N ) ) {
         INFO = -14;
      } else if ( LWORK < LWKMIN && !LQUERY ) {
         INFO = -16;
      }

      if ( INFO == 0 ) {
         NB1 = ILAENV( 1, 'SGEQRF', ' ', N, N, -1, -1 );
         NB2 = ILAENV( 1, 'SORMQR', ' ', N, N, N, -1 );
         NB3 = ILAENV( 1, 'SORGQR', ' ', N, N, N, -1 );
         NB = MAX( NB1, NB2, NB3 );
         LOPT = 2*N+N*(NB+1);
         WORK( 1 ) = LOPT;
      }

      if ( INFO != 0 ) {
         xerbla('SGEGS ', -INFO );
         RETURN;
      } else if ( LQUERY ) {
         RETURN;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Get machine constants

      EPS = SLAMCH( 'E' )*SLAMCH( 'B' );
      SAFMIN = SLAMCH( 'S' );
      SMLNUM = N*SAFMIN / EPS;
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

      if ( ILASCL ) {
         slascl('G', -1, -1, ANRM, ANRMTO, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            RETURN;
         }
      }

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

      if ( ILBSCL ) {
         slascl('G', -1, -1, BNRM, BNRMTO, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            RETURN;
         }
      }

      // Permute the matrix to make it more nearly triangular
      // Workspace layout:  (2*N words -- "work..." not actually used)
         // left_permutation, right_permutation, work...

      ILEFT = 1;
      IRIGHT = N + 1;
      IWORK = IRIGHT + N;
      sggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWORK ), IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 1;
         GO TO 10;
      }

      // Reduce B to triangular form, and initialize VSL and/or VSR
      // Workspace layout:  ("work..." must have at least N words)
         // left_permutation, right_permutation, tau, work...

      IROWS = IHI + 1 - ILO;
      ICOLS = N + 1 - ILO;
      ITAU = IWORK;
      IWORK = ITAU + IROWS;
      sgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO >= 0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 2;
         GO TO 10;
      }

      sormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
      if ( IINFO != 0 ) {
         INFO = N + 3;
         GO TO 10;
      }

      if ( ILVSL ) {
         slaset('Full', N, N, ZERO, ONE, VSL, LDVSL );
         slacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL );
         sorgqr(IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO );
         if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
         if ( IINFO != 0 ) {
            INFO = N + 4;
            GO TO 10;
         }
      }

      if (ILVSR) CALL SLASET( 'Full', N, N, ZERO, ONE, VSR, LDVSR );

      // Reduce to generalized Hessenberg form

      sgghrd(JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IINFO );
      if ( IINFO != 0 ) {
         INFO = N + 5;
         GO TO 10;
      }

      // Perform QZ algorithm, computing Schur vectors if desired
      // Workspace layout:  ("work..." must have at least 1 word)
         // left_permutation, right_permutation, work...

      IWORK = ITAU;
      shgeqz('S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWORK ), LWORK+1-IWORK, IINFO );
      if (IINFO >= 0) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 );
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
         sggbak('P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSL, LDVSL, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 7;
            GO TO 10;
         }
      }
      if ( ILVSR ) {
         sggbak('P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSR, LDVSR, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 8;
            GO TO 10;
         }
      }

      // Undo scaling

      if ( ILASCL ) {
         slascl('H', -1, -1, ANRMTO, ANRM, N, N, A, LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            RETURN;
         }
         slascl('G', -1, -1, ANRMTO, ANRM, N, 1, ALPHAR, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            RETURN;
         }
         slascl('G', -1, -1, ANRMTO, ANRM, N, 1, ALPHAI, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            RETURN;
         }
      }

      if ( ILBSCL ) {
         slascl('U', -1, -1, BNRMTO, BNRM, N, N, B, LDB, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            RETURN;
         }
         slascl('G', -1, -1, BNRMTO, BNRM, N, 1, BETA, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = N + 9;
            RETURN;
         }
      }

      } // 10
      WORK( 1 ) = LWKOPT;

      RETURN;

      // End of SGEGS

      }
