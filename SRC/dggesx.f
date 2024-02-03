      SUBROUTINE DGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, INFO )

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
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), RCONDE( 2 ), RCONDV( 2 ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * );
      // ..
      // .. Function Arguments ..
      bool               SELCTG;
      // EXTERNAL SELCTG
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, LST2SL, WANTSB, WANTSE, WANTSN, WANTST, WANTSV;
      int                I, ICOLS, IERR, IHI, IJOB, IJOBVL, IJOBVR, ILEFT, ILO, IP, IRIGHT, IROWS, ITAU, IWRK, LIWMIN, LWRK, MAXWRK, MINWRK;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PL, PR, SAFMAX, SAFMIN, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             DIF( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, DLASCL, DLASET, DORGQR, DORMQR, DTGSEN, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
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
         INFO = -16
      } else if ( LDVSR.LT.1 || ( ILVSR && LDVSR.LT.N ) ) {
         INFO = -18
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      if ( INFO == 0 ) {
         if ( N.GT.0) {
            MINWRK = MAX( 8*N, 6*N + 16 )
            MAXWRK = MINWRK - N + N*ILAENV( 1, 'DGEQRF', ' ', N, 1, N, 0 )             MAXWRK = MAX( MAXWRK, MINWRK - N + N*ILAENV( 1, 'DORMQR', ' ', N, 1, N, -1 ) )
            if ( ILVSL ) {
               MAXWRK = MAX( MAXWRK, MINWRK - N + N*ILAENV( 1, 'DORGQR', ' ', N, 1, N, -1 ) )
            }
            LWRK = MAXWRK
            if (IJOB.GE.1) LWRK = MAX( LWRK, N*N/2 );
         } else {
            MINWRK = 1
            MAXWRK = 1
            LWRK   = 1
         }
         WORK( 1 ) = LWRK
         if ( WANTSN || N == 0 ) {
            LIWMIN = 1
         } else {
            LIWMIN = N + 6
         }
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.MINWRK && .NOT.LQUERY ) {
            INFO = -22
         } else if ( LIWORK.LT.LIWMIN && .NOT.LQUERY ) {
            INFO = -24
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGGESX', -INFO );
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

      EPS = DLAMCH( 'P' )
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SMLNUM = SQRT( SAFMIN ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ILASCL = false;
      if ( ANRM.GT.ZERO && ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = true;
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = true;
      }
      if (ILASCL) CALL DLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR );

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      ILBSCL = false;
      if ( BNRM.GT.ZERO && BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = true;
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = true;
      }
      if (ILBSCL) CALL DLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR );

      // Permute the matrix to make it more nearly triangular
      // (Workspace: need 6*N + 2*N for permutation parameters)

      ILEFT = 1
      IRIGHT = N + 1
      IWRK = IRIGHT + N
      dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)
      // (Workspace: need N, prefer N*NB)

      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = IWRK
      IWRK = ITAU + IROWS
      dgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to matrix A
      // (Workspace: need N, prefer N*NB)

      dormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VSL
      // (Workspace: need N, prefer N*NB)

      if ( ILVSL ) {
         dlaset('Full', N, N, ZERO, ONE, VSL, LDVSL );
         if ( IROWS.GT.1 ) {
            dlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL );
         }
         dorgqr(IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VSR

      if (ILVSR) CALL DLASET( 'Full', N, N, ZERO, ONE, VSR, LDVSR );

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      dgghrd(JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IERR );

      SDIM = 0

      // Perform QZ algorithm, computing Schur vectors if desired
      // (Workspace: need N)

      IWRK = ITAU
      dhgeqz('S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, IERR );
      if ( IERR != 0 ) {
         if ( IERR.GT.0 && IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N && IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 60
      }

      // Sort eigenvalues ALPHA/BETA and compute the reciprocal of
      // condition number(s)
      // (Workspace: If IJOB >= 1, need MAX( 8*(N+1), 2*SDIM*(N-SDIM) )
                  // otherwise, need 8*(N+1) )

      if ( WANTST ) {

         // Undo scaling on eigenvalues before SELCTGing

         if ( ILASCL ) {
            dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
            dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
         }
         if (ILBSCL) CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );

         // Select eigenvalues

         for (I = 1; I <= N; I++) { // 10
            BWORK( I ) = SELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) )
         } // 10

         // Reorder eigenvalues, transform Generalized Schur vectors, and
         // compute reciprocal condition numbers

         dtgsen(IJOB, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PL, PR, DIF, WORK( IWRK ), LWORK-IWRK+1, IWORK, LIWORK, IERR );

         if (IJOB.GE.1) MAXWRK = MAX( MAXWRK, 2*SDIM*( N-SDIM ) );
         if ( IERR == -22 ) {

             // not enough real workspace

            INFO = -22
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

      if (ILVSL) CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSL, LDVSL, IERR );

      if (ILVSR) CALL DGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSR, LDVSR, IERR );

      // Check if unscaling would cause over/underflow, if so, rescale
      // (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
      // B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)

      if ( ILASCL ) {
         for (I = 1; I <= N; I++) { // 20
            if ( ALPHAI( I ) != ZERO ) {
               if ( ( ALPHAR( I ) / SAFMAX ).GT.( ANRMTO / ANRM ) || ( SAFMIN / ALPHAR( I ) ).GT.( ANRM / ANRMTO ) ) {
                  WORK( 1 ) = ABS( A( I, I ) / ALPHAR( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               } else if ( ( ALPHAI( I ) / SAFMAX ).GT. ( ANRMTO / ANRM ) || ( SAFMIN / ALPHAI( I ) ).GT.( ANRM / ANRMTO ) ) {
                  WORK( 1 ) = ABS( A( I, I+1 ) / ALPHAI( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               }
            }
         } // 20
      }

      if ( ILBSCL ) {
         for (I = 1; I <= N; I++) { // 30
            if ( ALPHAI( I ) != ZERO ) {
               if ( ( BETA( I ) / SAFMAX ).GT.( BNRMTO / BNRM ) || ( SAFMIN / BETA( I ) ).GT.( BNRM / BNRMTO ) ) {
                  WORK( 1 ) = ABS( B( I, I ) / BETA( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
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
         SDIM = 0
         IP = 0
         for (I = 1; I <= N; I++) { // 50
            CURSL = SELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) )
            if ( ALPHAI( I ) == ZERO ) {
               if (CURSL) SDIM = SDIM + 1;
               IP = 0
               if (CURSL && .NOT.LASTSL) INFO = N + 2;
            } else {
               if ( IP == 1 ) {

                  // Last eigenvalue of conjugate pair

                  CURSL = CURSL || LASTSL
                  LASTSL = CURSL
                  if (CURSL) SDIM = SDIM + 2;
                  IP = -1
                  if (CURSL && .NOT.LST2SL) INFO = N + 2;
               } else {

                  // First eigenvalue of conjugate pair

                  IP = 1
               }
            }
            LST2SL = LASTSL
            LASTSL = CURSL
         } // 50

      }

      } // 60

      WORK( 1 ) = MAXWRK
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of DGGESX

      }
