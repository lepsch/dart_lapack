      SUBROUTINE CGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, RCONDE, RCONDV, WORK, LWORK, RWORK, IWORK, LIWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR, SENSE, SORT;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LIWORK, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      REAL               RCONDE( 2 ), RCONDV( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * )
      // ..
      // .. Function Arguments ..
      bool               SELCTG;
      // EXTERNAL SELCTG
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, WANTSB, WANTSE, WANTSN, WANTST, WANTSV;
      int                I, ICOLS, IERR, IHI, IJOB, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, IRWRK, ITAU, IWRK, LIWMIN, LWRK, MAXWRK, MINWRK;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PL, PR, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               DIF( 2 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CGGBAK, CGGBAL, CGGHRD, CHGEQZ, CLACPY, CLASCL, CLASET, CTGSEN, CUNGQR, CUNMQR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               CLANGE, SLAMCH, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, CLANGE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

      if ( LSAME( JOBVSL, 'N' ) ) {
         IJOBVL = 1
         ILVSL = false;
      } else if ( LSAME( JOBVSL, 'V' ) ) {
         IJOBVL = 2
         ILVSL = true;
      } else {
         IJOBVL = -1
         ILVSL = false;
      }

      if ( LSAME( JOBVSR, 'N' ) ) {
         IJOBVR = 1
         ILVSR = false;
      } else if ( LSAME( JOBVSR, 'V' ) ) {
         IJOBVR = 2
         ILVSR = true;
      } else {
         IJOBVR = -1
         ILVSR = false;
      }

      WANTST = LSAME( SORT, 'S' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
      LQUERY = ( LWORK == -1 || LIWORK == -1 )
      if ( WANTSN ) {
         IJOB = 0
      } else if ( WANTSE ) {
         IJOB = 1
      } else if ( WANTSV ) {
         IJOB = 2
      } else if ( WANTSB ) {
         IJOB = 4
      }

      // Test the input arguments

      INFO = 0
      if ( IJOBVL.LE.0 ) {
         INFO = -1
      } else if ( IJOBVR.LE.0 ) {
         INFO = -2
      } else if ( ( .NOT.WANTST ) && ( .NOT.LSAME( SORT, 'N' ) ) ) {
         INFO = -3
      } else if ( .NOT.( WANTSN || WANTSE || WANTSV || WANTSB ) || ( .NOT.WANTST && .NOT.WANTSN ) ) {
         INFO = -5
      } else if ( N.LT.0 ) {
         INFO = -6
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -10
      } else if ( LDVSL.LT.1 || ( ILVSL && LDVSL.LT.N ) ) {
         INFO = -15
      } else if ( LDVSR.LT.1 || ( ILVSR && LDVSR.LT.N ) ) {
         INFO = -17
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      if ( INFO == 0 ) {
         if ( N.GT.0) {
            MINWRK = 2*N
            MAXWRK = N*(1 + ILAENV( 1, 'CGEQRF', ' ', N, 1, N, 0 ) )
            MAXWRK = MAX( MAXWRK, N*( 1 + ILAENV( 1, 'CUNMQR', ' ', N, 1, N, -1 ) ) )
            if ( ILVSL ) {
               MAXWRK = MAX( MAXWRK, N*( 1 + ILAENV( 1, 'CUNGQR', ' ', N, 1, N, -1 ) ) )
            }
            LWRK = MAXWRK
            if (IJOB.GE.1) LWRK = MAX( LWRK, N*N/2 );
         } else {
            MINWRK = 1
            MAXWRK = 1
            LWRK   = 1
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWRK)
         if ( WANTSN || N == 0 ) {
            LIWMIN = 1
         } else {
            LIWMIN = N + 2
         }
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.MINWRK && .NOT.LQUERY ) {
            INFO = -21
         } else if ( LIWORK.LT.LIWMIN && .NOT.LQUERY) {
            INFO = -24
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGGESX', -INFO );
         RETURN
      } else if (LQUERY) {
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         SDIM = 0
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = false;
      if ( ANRM.GT.ZERO && ANRM.LT.SMLNUM ) {
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
      if ( BNRM.GT.ZERO && BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = true;
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = true;
      }
      if (ILBSCL) CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute the matrix to make it more nearly triangular
      // (Real Workspace: need 6*N)

      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      cggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)
      // (Complex Workspace: need N, prefer N*NB)

      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = 1
      IWRK = ITAU + IROWS
      cgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the unitary transformation to matrix A
      // (Complex Workspace: need N, prefer N*NB)

      cunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VSL
      // (Complex Workspace: need N, prefer N*NB)

      if ( ILVSL ) {
         claset('Full', N, N, CZERO, CONE, VSL, LDVSL );
         if ( IROWS.GT.1 ) {
            clacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL );
         }
         cungqr(IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VSR

      if (ILVSR) CALL CLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR );

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      cgghrd(JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IERR );

      SDIM = 0

      // Perform QZ algorithm, computing Schur vectors if desired
      // (Complex Workspace: need N)
      // (Real Workspace:    need N)

      IWRK = ITAU
      chgeqz('S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, RWORK( IRWRK ), IERR );
      if ( IERR != 0 ) {
         if ( IERR.GT.0 && IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N && IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 40
      }

      // Sort eigenvalues ALPHA/BETA and compute the reciprocal of
      // condition number(s)

      if ( WANTST ) {

         // Undo scaling on eigenvalues before SELCTGing

         if (ILASCL) CALL CLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )          IF( ILBSCL ) CALL CLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );

         // Select eigenvalues

         for (I = 1; I <= N; I++) { // 10
            BWORK( I ) = SELCTG( ALPHA( I ), BETA( I ) )
         } // 10

         // Reorder eigenvalues, transform Generalized Schur vectors, and
         // compute reciprocal condition numbers
         // (Complex Workspace: If IJOB >= 1, need MAX(1, 2*SDIM*(N-SDIM))
                             // otherwise, need 1 )

         ctgsen(IJOB, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PL, PR, DIF, WORK( IWRK ), LWORK-IWRK+1, IWORK, LIWORK, IERR );

         if (IJOB.GE.1) MAXWRK = MAX( MAXWRK, 2*SDIM*( N-SDIM ) );
         if ( IERR == -21 ) {

             // not enough complex workspace

            INFO = -21
         } else {
            if ( IJOB == 1 || IJOB == 4 ) {
               RCONDE( 1 ) = PL
               RCONDE( 2 ) = PR
            }
            if ( IJOB == 2 || IJOB == 4 ) {
               RCONDV( 1 ) = DIF( 1 )
               RCONDV( 2 ) = DIF( 2 )
            }
            if (IERR == 1) INFO = N + 3;
         }

      }

      // Apply permutation to VSL and VSR
      // (Workspace: none needed)

      if (ILVSL) CALL CGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSL, LDVSL, IERR );

      if (ILVSR) CALL CGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSR, LDVSR, IERR );

      // Undo scaling

      if ( ILASCL ) {
         clascl('U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR );
         clascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR );
      }

      if ( ILBSCL ) {
         clascl('U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR );
         clascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      if ( WANTST ) {

         // Check if reordering is correct

         LASTSL = true;
         SDIM = 0
         for (I = 1; I <= N; I++) { // 30
            CURSL = SELCTG( ALPHA( I ), BETA( I ) )
            if (CURSL) SDIM = SDIM + 1             IF( CURSL && .NOT.LASTSL ) INFO = N + 2;
            LASTSL = CURSL
         } // 30

      }

      } // 40

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of CGGESX

      }
