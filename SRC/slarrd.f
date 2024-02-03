      SUBROUTINE SLARRD( RANGE, ORDER, N, VL, VU, IL, IU, GERS, RELTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT, M, W, WERR, WL, WU, IBLOCK, INDEXW, WORK, IWORK, INFO );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             ORDER, RANGE;
      int                IL, INFO, IU, M, N, NSPLIT;
      REAL                PIVMIN, RELTOL, VL, VU, WL, WU;
      // ..
      // .. Array Arguments ..
      int                IBLOCK( * ), INDEXW( * ), ISPLIT( * ), IWORK( * );
      REAL               D( * ), E( * ), E2( * ), GERS( * ), W( * ), WERR( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, HALF, FUDGE;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = ONE/TWO, FUDGE = TWO ;
      int       ALLRNG, VALRNG, INDRNG;
      const     ALLRNG = 1, VALRNG = 2, INDRNG = 3 ;
      // ..
      // .. Local Scalars ..
      bool               NCNVRG, TOOFEW;
      int                I, IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO, IM, IN, IOFF, IOUT, IRANGE, ITMAX, ITMP1, ITMP2, IW, IWOFF, J, JBLK, JDISC, JE, JEE, NB, NWL, NWU;
      REAL               ATOLI, EPS, GL, GU, RTOLI, TMP1, TMP2, TNORM, UFLOW, WKILL, WLU, WUL;

      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH;
      // EXTERNAL LSAME, ILAENV, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAEBZ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0;
      M = 0;

      // Quick return if possible

      if ( N <= 0 ) {
         RETURN;
      }

      // Decode RANGE

      if ( LSAME( RANGE, 'A' ) ) {
         IRANGE = ALLRNG;
      } else if ( LSAME( RANGE, 'V' ) ) {
         IRANGE = VALRNG;
      } else if ( LSAME( RANGE, 'I' ) ) {
         IRANGE = INDRNG;
      } else {
         IRANGE = 0;
      }

      // Check for Errors

      if ( IRANGE <= 0 ) {
         INFO = -1;
      } else if ( !(LSAME(ORDER,'B') || LSAME(ORDER,'E')) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( IRANGE == VALRNG ) {
         if ( VL >= VU ) INFO = -5       ELSE IF( IRANGE == INDRNG && ( IL < 1 || IL > MAX( 1, N ) ) ) {
         INFO = -6;
      } else if ( IRANGE == INDRNG && ( IU < MIN( N, IL ) || IU > N ) ) {
         INFO = -7;
      }

      if ( INFO != 0 ) {
         RETURN;
      }

      // Initialize error flags
      NCNVRG = false;
      TOOFEW = false;

      // Simplification:
      if (IRANGE == INDRNG && IL == 1 && IU == N) IRANGE = 1;

      // Get machine constants
      EPS = SLAMCH( 'P' );
      UFLOW = SLAMCH( 'U' );


      // Special Case when N=1
      // Treat case of 1x1 matrix for quick return
      if ( N == 1 ) {
         if ( (IRANGE == ALLRNG) || ((IRANGE == VALRNG) && (D(1) > VL) && (D(1) <= VU)) || ((IRANGE == INDRNG) && (IL == 1) && (IU == 1)) ) {
            M = 1;
            W(1) = D(1);
            // The computation error of the eigenvalue is zero
            WERR(1) = ZERO;
            IBLOCK( 1 ) = 1;
            INDEXW( 1 ) = 1;
         }
         RETURN;
      }

      // NB is the minimum vector length for vector bisection, or 0
      // if only scalar is to be done.
      NB = ILAENV( 1, 'SSTEBZ', ' ', N, -1, -1, -1 );
      if (NB <= 1) NB = 0;

      // Find global spectral radius
      GL = D(1);
      GU = D(1);
      for (I = 1; I <= N; I++) { // 5
         GL =  MIN( GL, GERS( 2*I - 1));
         GU = MAX( GU, GERS(2*I) );
      } // 5
      // Compute global Gerschgorin bounds and spectral diameter
      TNORM = MAX( ABS( GL ), ABS( GU ) );
      GL = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN;
      GU = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN;
      // [JAN/28/2009] remove the line below since SPDIAM variable not use
      // SPDIAM = GU - GL
      // Input arguments for SLAEBZ:
      // The relative tolerance.  An interval (a,b] lies within
      // "relative tolerance" if  b-a < RELTOL*max(|a|,|b|),
      RTOLI = RELTOL;
      // Set the absolute tolerance for interval convergence to zero to force
      // interval convergence based on relative size of the interval.
      // This is dangerous because intervals might not converge when RELTOL is
      // small. But at least a very small number should be selected so that for
      // strongly graded matrices, the code can get relatively accurate
      // eigenvalues.
      ATOLI = FUDGE*TWO*UFLOW + FUDGE*TWO*PIVMIN;

      if ( IRANGE == INDRNG ) {

         // RANGE='I': Compute an interval containing eigenvalues
         // IL through IU. The initial interval [GL,GU] from the global
         // Gerschgorin bounds GL and GU is refined by SLAEBZ.
         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2;
         WORK( N+1 ) = GL;
         WORK( N+2 ) = GL;
         WORK( N+3 ) = GU;
         WORK( N+4 ) = GU;
         WORK( N+5 ) = GL;
         WORK( N+6 ) = GU;
         IWORK( 1 ) = -1;
         IWORK( 2 ) = -1;
         IWORK( 3 ) = N + 1;
         IWORK( 4 ) = N + 1;
         IWORK( 5 ) = IL - 1;
         IWORK( 6 ) = IU;

         slaebz(3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E, E2, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT, IWORK, W, IBLOCK, IINFO );
         if ( IINFO != 0 ) {
            INFO = IINFO;
            RETURN;
         }
         // On exit, output intervals may not be ordered by ascending negcount
         if ( IWORK( 6 ) == IU ) {
            WL = WORK( N+1 );
            WLU = WORK( N+3 );
            NWL = IWORK( 1 );
            WU = WORK( N+4 );
            WUL = WORK( N+2 );
            NWU = IWORK( 4 );
         } else {
            WL = WORK( N+2 );
            WLU = WORK( N+4 );
            NWL = IWORK( 2 );
            WU = WORK( N+3 );
            WUL = WORK( N+1 );
            NWU = IWORK( 3 );
         }
         // On exit, the interval [WL, WLU] contains a value with negcount NWL,
         // and [WUL, WU] contains a value with negcount NWU.
         if ( NWL < 0 || NWL >= N || NWU < 1 || NWU > N ) {
            INFO = 4;
            RETURN;
         }

      } else if ( IRANGE == VALRNG ) {
         WL = VL;
         WU = VU;

      } else if ( IRANGE == ALLRNG ) {
         WL = GL;
         WU = GU;
      }



      // Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU.
      // NWL accumulates the number of eigenvalues <= WL,
      // NWU accumulates the number of eigenvalues <= WU
      M = 0;
      IEND = 0;
      INFO = 0;
      NWL = 0;
      NWU = 0;

      for (JBLK = 1; JBLK <= NSPLIT; JBLK++) { // 70
         IOFF = IEND;
         IBEGIN = IOFF + 1;
         IEND = ISPLIT( JBLK );
         IN = IEND - IOFF;

         if ( IN == 1 ) {
            // 1x1 block
            if ( WL >= D( IBEGIN )-PIVMIN ) NWL = NWL + 1             IF( WU >= D( IBEGIN )-PIVMIN ) NWU = NWU + 1             IF( IRANGE == ALLRNG || ( WL < D( IBEGIN )-PIVMIN && WU >= D( IBEGIN )-PIVMIN ) ) {
               M = M + 1;
               W( M ) = D( IBEGIN );
               WERR(M) = ZERO;
               // The gap for a single block doesn't matter for the later
               // algorithm and is assigned an arbitrary large value
               IBLOCK( M ) = JBLK;
               INDEXW( M ) = 1;
            }

         // Disabled 2x2 case because of a failure on the following matrix
         // RANGE = 'I', IL = IU = 4
           // Original Tridiagonal, d = [
            // -0.150102010615740e+00
            // -0.849897989384260e+00
            // -0.128208148052635e-15
             // 0.128257718286320e-15
           // ];
           // e = [
            // -0.357171383266986e+00
            // -0.180411241501588e-15
            // -0.175152352710251e-15
           // ];

          // ELSE IF( IN == 2 ) THEN
**           2x2 block
             // DISC = SQRT( (HALF*(D(IBEGIN)-D(IEND)))**2 + E(IBEGIN)**2 )
             // TMP1 = HALF*(D(IBEGIN)+D(IEND))
             // L1 = TMP1 - DISC
             // IF( WL >= L1-PIVMIN )
      // $         NWL = NWL + 1
             // IF( WU >= L1-PIVMIN )
      // $         NWU = NWU + 1
             // IF( IRANGE == ALLRNG || ( WL < L1-PIVMIN && WU.GE.
      // $          L1-PIVMIN ) ) THEN
                // M = M + 1
                // W( M ) = L1
**              The uncertainty of eigenvalues of a 2x2 matrix is very small
                // WERR( M ) = EPS * ABS( W( M ) ) * TWO
                // IBLOCK( M ) = JBLK
                // INDEXW( M ) = 1
             // ENDIF
             // L2 = TMP1 + DISC
             // IF( WL >= L2-PIVMIN )
      // $         NWL = NWL + 1
             // IF( WU >= L2-PIVMIN )
      // $         NWU = NWU + 1
             // IF( IRANGE == ALLRNG || ( WL < L2-PIVMIN && WU.GE.
      // $          L2-PIVMIN ) ) THEN
                // M = M + 1
                // W( M ) = L2
**              The uncertainty of eigenvalues of a 2x2 matrix is very small
                // WERR( M ) = EPS * ABS( W( M ) ) * TWO
                // IBLOCK( M ) = JBLK
                // INDEXW( M ) = 2
             // ENDIF
         } else {
            // General Case - block of size IN >= 2
            // Compute local Gerschgorin interval and use it as the initial
            // interval for SLAEBZ
            GU = D( IBEGIN );
            GL = D( IBEGIN );
            TMP1 = ZERO;

            for (J = IBEGIN; J <= IEND; J++) { // 40
               GL =  MIN( GL, GERS( 2*J - 1));
               GU = MAX( GU, GERS(2*J) );
            } // 40
            // [JAN/28/2009]
            // change SPDIAM by TNORM in lines 2 and 3 thereafter
            // line 1: remove computation of SPDIAM (not useful anymore)
            // SPDIAM = GU - GL
            // GL = GL - FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN
            // GU = GU + FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN
            GL = GL - FUDGE*TNORM*EPS*IN - FUDGE*PIVMIN;
            GU = GU + FUDGE*TNORM*EPS*IN + FUDGE*PIVMIN;

            if ( IRANGE > 1 ) {
               if ( GU < WL ) {
                  // the local block contains none of the wanted eigenvalues
                  NWL = NWL + IN;
                  NWU = NWU + IN;
                  GO TO 70;
               }
               // refine search interval if possible, only range (WL,WU] matters
               GL = MAX( GL, WL );
               GU = MIN( GU, WU );
               if (GL >= GU) GO TO 70;
            }

            // Find negcount of initial interval boundaries GL and GU
            WORK( N+1 ) = GL;
            WORK( N+IN+1 ) = GU;
            slaebz(1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D( IBEGIN ), E( IBEGIN ), E2( IBEGIN ), IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM, IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO );
            if ( IINFO != 0 ) {
               INFO = IINFO;
               RETURN;
            }

            NWL = NWL + IWORK( 1 );
            NWU = NWU + IWORK( IN+1 );
            IWOFF = M - IWORK( 1 );

            // Compute Eigenvalues
            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2;
            slaebz(2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D( IBEGIN ), E( IBEGIN ), E2( IBEGIN ), IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT, IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO );
            if ( IINFO != 0 ) {
               INFO = IINFO;
               RETURN;
            }

            // Copy eigenvalues into W and IBLOCK
            // Use -JBLK for block number for unconverged eigenvalues.
            // Loop over the number of output intervals from SLAEBZ
            for (J = 1; J <= IOUT; J++) { // 60
               // eigenvalue approximation is middle point of interval
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) );
               // semi length of error interval
               TMP2 = HALF*ABS( WORK( J+N )-WORK( J+IN+N ) );
               if ( J > IOUT-IINFO ) {
                  // Flag non-convergence.
                  NCNVRG = true;
                  IB = -JBLK;
               } else {
                  IB = JBLK;
               }
               for (JE = IWORK( J ) + 1 + IWOFF; JE <= IWORK( J+IN ) + IWOFF; JE++) { // 50
                  W( JE ) = TMP1;
                  WERR( JE ) = TMP2;
                  INDEXW( JE ) = JE - IWOFF;
                  IBLOCK( JE ) = IB;
               } // 50
            } // 60

            M = M + IM;
         }
      } // 70

      // If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
      // If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
      if ( IRANGE == INDRNG ) {
         IDISCL = IL - 1 - NWL;
         IDISCU = NWU - IU;

         if ( IDISCL > 0 ) {
            IM = 0;
            for (JE = 1; JE <= M; JE++) { // 80
               // Remove some of the smallest eigenvalues from the left so that
               // at the end IDISCL =0. Move all eigenvalues up to the left.
               if ( W( JE ) <= WLU && IDISCL > 0 ) {
                  IDISCL = IDISCL - 1;
               } else {
                  IM = IM + 1;
                  W( IM ) = W( JE );
                  WERR( IM ) = WERR( JE );
                  INDEXW( IM ) = INDEXW( JE );
                  IBLOCK( IM ) = IBLOCK( JE );
               }
            } // 80
            M = IM;
         }
         if ( IDISCU > 0 ) {
            // Remove some of the largest eigenvalues from the right so that
            // at the end IDISCU =0. Move all eigenvalues up to the left.
            IM=M+1;
            DO 81 JE = M, 1, -1;
               if ( W( JE ) >= WUL && IDISCU > 0 ) {
                  IDISCU = IDISCU - 1;
               } else {
                  IM = IM - 1;
                  W( IM ) = W( JE );
                  WERR( IM ) = WERR( JE );
                  INDEXW( IM ) = INDEXW( JE );
                  IBLOCK( IM ) = IBLOCK( JE );
               }
            } // 81
            JEE = 0;
            for (JE = IM; JE <= M; JE++) { // 82
               JEE = JEE + 1;
               W( JEE ) = W( JE );
               WERR( JEE ) = WERR( JE );
               INDEXW( JEE ) = INDEXW( JE );
               IBLOCK( JEE ) = IBLOCK( JE );
            } // 82
            M = M-IM+1;
         }

         if ( IDISCL > 0 || IDISCU > 0 ) {
            // Code to deal with effects of bad arithmetic. (If N(w) is
            // monotone non-decreasing, this should never happen.)
            // Some low eigenvalues to be discarded are not in (WL,WLU],
            // or high eigenvalues to be discarded are not in (WUL,WU]
            // so just kill off the smallest IDISCL/largest IDISCU
            // eigenvalues, by marking the corresponding IBLOCK = 0
            if ( IDISCL > 0 ) {
               WKILL = WU;
               for (JDISC = 1; JDISC <= IDISCL; JDISC++) { // 100
                  IW = 0;
                  for (JE = 1; JE <= M; JE++) { // 90
                     if ( IBLOCK( JE ) != 0 && ( W( JE ) < WKILL || IW == 0 ) ) {
                        IW = JE;
                        WKILL = W( JE );
                     }
                  } // 90
                  IBLOCK( IW ) = 0;
               } // 100
            }
            if ( IDISCU > 0 ) {
               WKILL = WL;
               for (JDISC = 1; JDISC <= IDISCU; JDISC++) { // 120
                  IW = 0;
                  for (JE = 1; JE <= M; JE++) { // 110
                     if ( IBLOCK( JE ) != 0 && ( W( JE ) >= WKILL || IW == 0 ) ) {
                        IW = JE;
                        WKILL = W( JE );
                     }
                  } // 110
                  IBLOCK( IW ) = 0;
               } // 120
            }
            // Now erase all eigenvalues with IBLOCK set to zero
            IM = 0;
            for (JE = 1; JE <= M; JE++) { // 130
               if ( IBLOCK( JE ) != 0 ) {
                  IM = IM + 1;
                  W( IM ) = W( JE );
                  WERR( IM ) = WERR( JE );
                  INDEXW( IM ) = INDEXW( JE );
                  IBLOCK( IM ) = IBLOCK( JE );
               }
            } // 130
            M = IM;
         }
         if ( IDISCL < 0 || IDISCU < 0 ) {
            TOOFEW = true;
         }
      }

      if (( IRANGE == ALLRNG && M != N ) || ( IRANGE == INDRNG && M != IU-IL+1 ) ) {
         TOOFEW = true;
      }

      // If ORDER='B', do nothing the eigenvalues are already sorted by
         // block.
      // If ORDER='E', sort the eigenvalues from smallest to largest

      if ( LSAME(ORDER,'E') && NSPLIT > 1 ) {
         for (JE = 1; JE <= M - 1; JE++) { // 150
            IE = 0;
            TMP1 = W( JE );
            for (J = JE + 1; J <= M; J++) { // 140
               if ( W( J ) < TMP1 ) {
                  IE = J;
                  TMP1 = W( J );
               }
            } // 140
            if ( IE != 0 ) {
               TMP2 = WERR( IE );
               ITMP1 = IBLOCK( IE );
               ITMP2 = INDEXW( IE );
               W( IE ) = W( JE );
               WERR( IE ) = WERR( JE );
               IBLOCK( IE ) = IBLOCK( JE );
               INDEXW( IE ) = INDEXW( JE );
               W( JE ) = TMP1;
               WERR( JE ) = TMP2;
               IBLOCK( JE ) = ITMP1;
               INDEXW( JE ) = ITMP2;
            }
         } // 150
      }

      INFO = 0;
      if (NCNVRG) INFO = INFO + 1       IF( TOOFEW ) INFO = INFO + 2;
      RETURN;

      // End of SLARRD

      }
