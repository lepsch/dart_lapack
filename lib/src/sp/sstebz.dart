      void sstebz(final int RANGE, final int ORDER, final int N, final int VL, final int VU, final int IL, final int IU, final int ABSTOL, final int D, final int E, final int M, final int NSPLIT, final int W, final int IBLOCK, final int ISPLIT, final Array<double> _WORK_, final Array<int> IWORK_, final Box<int> INFO,) {
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             ORDER, RANGE;
      int                IL, INFO, IU, M, N, NSPLIT;
      double               ABSTOL, VL, VU;
      int                IBLOCK( * ), ISPLIT( * ), IWORK( * );
      double               D( * ), E( * ), W( * ), WORK( * );
      // ..

      double               ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 1.0 / TWO ;
      double               FUDGE, RELFAC;
      const              FUDGE = 2.1, RELFAC = 2.0 ;
      bool               NCNVRG, TOOFEW;
      int                IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO, IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX, ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL, NWU;
      double               ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN, TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL;
      int                IDUMMA( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ILAENV, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAEBZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, MIN, SQRT

      INFO = 0;

      // Decode RANGE

      if ( lsame( RANGE, 'A' ) ) {
         IRANGE = 1;
      } else if ( lsame( RANGE, 'V' ) ) {
         IRANGE = 2;
      } else if ( lsame( RANGE, 'I' ) ) {
         IRANGE = 3;
      } else {
         IRANGE = 0;
      }

      // Decode ORDER

      if ( lsame( ORDER, 'B' ) ) {
         IORDER = 2;
      } else if ( lsame( ORDER, 'E' ) ) {
         IORDER = 1;
      } else {
         IORDER = 0;
      }

      // Check for Errors

      if ( IRANGE <= 0 ) {
         INFO = -1;
      } else if ( IORDER <= 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( IRANGE == 2 ) {
         if (VL >= VU) INFO = -5;
      } else if ( IRANGE == 3 && ( IL < 1 || IL > max( 1, N ) ) ) {
         INFO = -6;
      } else if ( IRANGE == 3 && ( IU < min( N, IL ) || IU > N ) ) {
         INFO = -7;
      }

      if ( INFO != 0 ) {
         xerbla('SSTEBZ', -INFO );
         return;
      }

      // Initialize error flags

      INFO = 0;
      NCNVRG = false;
      TOOFEW = false;

      // Quick return if possible

      M = 0;
      if (N == 0) return;

      // Simplifications:

      if (IRANGE == 3 && IL == 1 && IU == N) IRANGE = 1;

      // Get machine constants
      // NB is the minimum vector length for vector bisection, or 0
      // if only scalar is to be done.

      SAFEMN = SLAMCH( 'S' );
      ULP = SLAMCH( 'P' );
      RTOLI = ULP*RELFAC;
      NB = ilaenv( 1, 'SSTEBZ', ' ', N, -1, -1, -1 );
      if (NB <= 1) NB = 0;

      // Special Case when N=1

      if ( N == 1 ) {
         NSPLIT = 1;
         ISPLIT[1] = 1;
         if ( IRANGE == 2 && ( VL >= D( 1 ) || VU < D( 1 ) ) ) {
            M = 0;
         } else {
            W[1] = D( 1 );
            IBLOCK[1] = 1;
            M = 1;
         }
         return;
      }

      // Compute Splitting Points

      NSPLIT = 1;
      WORK[N] = ZERO;
      PIVMIN = ONE;

      for (J = 2; J <= N; J++) { // 10
         TMP1 = E( J-1 )**2;
         if ( ABS( D( J )*D( J-1 ) )*ULP**2+SAFEMN > TMP1 ) {
            ISPLIT[NSPLIT] = J - 1;
            NSPLIT = NSPLIT + 1;
            WORK[J-1] = ZERO;
         } else {
            WORK[J-1] = TMP1;
            PIVMIN = max( PIVMIN, TMP1 );
         }
      } // 10
      ISPLIT[NSPLIT] = N;
      PIVMIN = PIVMIN*SAFEMN;

      // Compute Interval and ATOLI

      if ( IRANGE == 3 ) {

         // RANGE='I': Compute the interval containing eigenvalues
         //            IL through IU.

         // Compute Gershgorin interval for entire (split) matrix
         // and use it as the initial interval

         GU = D( 1 );
         GL = D( 1 );
         TMP1 = ZERO;

         for (J = 1; J <= N - 1; J++) { // 20
            TMP2 = sqrt( WORK( J ) );
            GU = max( GU, D( J )+TMP1+TMP2 );
            GL = min( GL, D( J )-TMP1-TMP2 );
            TMP1 = TMP2;
         } // 20

         GU = max( GU, D( N )+TMP1 );
         GL = min( GL, D( N )-TMP1 );
         TNORM = max( ( GL ).abs(), ( GU ).abs() );
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN;
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN;

         // Compute Iteration parameters

         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2;
         if ( ABSTOL <= ZERO ) {
            ATOLI = ULP*TNORM;
         } else {
            ATOLI = ABSTOL;
         }

         WORK[N+1] = GL;
         WORK[N+2] = GL;
         WORK[N+3] = GU;
         WORK[N+4] = GU;
         WORK[N+5] = GL;
         WORK[N+6] = GU;
         IWORK[1] = -1;
         IWORK[2] = -1;
         IWORK[3] = N + 1;
         IWORK[4] = N + 1;
         IWORK[5] = IL - 1;
         IWORK[6] = IU;

         slaebz(3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E, WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT, IWORK, W, IBLOCK, IINFO );

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

         if ( NWL < 0 || NWL >= N || NWU < 1 || NWU > N ) {
            INFO = 4;
            return;
         }
      } else {

         // RANGE='A' or 'V' -- Set ATOLI

         TNORM = max( ( D( 1 ) ).abs()+( E( 1 ) ).abs(), ( D( N ) ).abs()+( E( N-1 ) ).abs() );

         for (J = 2; J <= N - 1; J++) { // 30
            TNORM = max( TNORM, ( D( J ) ).abs()+( E( J-1 ) ).abs()+ ( E( J ) ).abs() );
         } // 30

         if ( ABSTOL <= ZERO ) {
            ATOLI = ULP*TNORM;
         } else {
            ATOLI = ABSTOL;
         }

         if ( IRANGE == 2 ) {
            WL = VL;
            WU = VU;
         } else {
            WL = ZERO;
            WU = ZERO;
         }
      }

      // Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
      // NWL accumulates the number of eigenvalues <= WL,
      // NWU accumulates the number of eigenvalues <= WU

      M = 0;
      IEND = 0;
      INFO = 0;
      NWL = 0;
      NWU = 0;

      for (JB = 1; JB <= NSPLIT; JB++) { // 70
         IOFF = IEND;
         IBEGIN = IOFF + 1;
         IEND = ISPLIT( JB );
         IN = IEND - IOFF;

         if ( IN == 1 ) {

            // Special Case -- IN=1

            if ( IRANGE == 1 || WL >= D( IBEGIN )-PIVMIN ) NWL = NWL + 1;
            if( IRANGE == 1 || WU >= D( IBEGIN )-PIVMIN ) NWU = NWU + 1;
            IF( IRANGE == 1 || ( WL < D( IBEGIN )-PIVMIN && WU >= D( IBEGIN )-PIVMIN ) ) {
               M = M + 1;
               W[M] = D( IBEGIN );
               IBLOCK[M] = JB;
            }
         } else {

            // General Case -- IN > 1

            // Compute Gershgorin Interval
            // and use it as the initial interval

            GU = D( IBEGIN );
            GL = D( IBEGIN );
            TMP1 = ZERO;

            for (J = IBEGIN; J <= IEND - 1; J++) { // 40
               TMP2 = ( E( J ) ).abs();
               GU = max( GU, D( J )+TMP1+TMP2 );
               GL = min( GL, D( J )-TMP1-TMP2 );
               TMP1 = TMP2;
            } // 40

            GU = max( GU, D( IEND )+TMP1 );
            GL = min( GL, D( IEND )-TMP1 );
            BNORM = max( ( GL ).abs(), ( GU ).abs() );
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN;
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN;

            // Compute ATOLI for the current submatrix

            if ( ABSTOL <= ZERO ) {
               ATOLI = ULP*max( ( GL ).abs(), ( GU ).abs() );
            } else {
               ATOLI = ABSTOL;
            }

            if ( IRANGE > 1 ) {
               if ( GU < WL ) {
                  NWL = NWL + IN;
                  NWU = NWU + IN;
                  GO TO 70;
               }
               GL = max( GL, WL );
               GU = min( GU, WU );
               if (GL >= GU) GO TO 70;
            }

            // Set Up Initial Interval

            WORK[N+1] = GL;
            WORK[N+IN+1] = GU;
            slaebz(1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM, IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO );

            NWL = NWL + IWORK( 1 );
            NWU = NWU + IWORK( IN+1 );
            IWOFF = M - IWORK( 1 );

            // Compute Eigenvalues

            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2;
            slaebz(2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT, IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO );

            // Copy Eigenvalues Into W and IBLOCK
            // Use -JB for block number for unconverged eigenvalues.

            for (J = 1; J <= IOUT; J++) { // 60
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) );

               // Flag non-convergence.

               if ( J > IOUT-IINFO ) {
                  NCNVRG = true;
                  IB = -JB;
               } else {
                  IB = JB;
               }
               for (JE = IWORK( J ) + 1 + IWOFF; JE <= IWORK( J+IN ) + IWOFF; JE++) { // 50
                  W[JE] = TMP1;
                  IBLOCK[JE] = IB;
               } // 50
            } // 60

            M = M + IM;
         }
      } // 70

      // If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
      // If NWL+1 < IL or NWU > IU, discard extra eigenvalues.

      if ( IRANGE == 3 ) {
         IM = 0;
         IDISCL = IL - 1 - NWL;
         IDISCU = NWU - IU;

         if ( IDISCL > 0 || IDISCU > 0 ) {
            for (JE = 1; JE <= M; JE++) { // 80
               if ( W( JE ) <= WLU && IDISCL > 0 ) {
                  IDISCL = IDISCL - 1;
               } else if ( W( JE ) >= WUL && IDISCU > 0 ) {
                  IDISCU = IDISCU - 1;
               } else {
                  IM = IM + 1;
                  W[IM] = W( JE );
                  IBLOCK[IM] = IBLOCK( JE );
               }
            } // 80
            M = IM;
         }
         if ( IDISCL > 0 || IDISCU > 0 ) {

            // Code to deal with effects of bad arithmetic:
            // Some low eigenvalues to be discarded are not in (WL,WLU],
            // or high eigenvalues to be discarded are not in (WUL,WU]
            // so just kill off the smallest IDISCL/largest IDISCU
            // eigenvalues, by simply finding the smallest/largest
            // eigenvalue(s).

            // (If N(w) is monotone non-decreasing, this should never
            //     happen.)

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
                  IBLOCK[IW] = 0;
               } // 100
            }
            if ( IDISCU > 0 ) {

               WKILL = WL;
               for (JDISC = 1; JDISC <= IDISCU; JDISC++) { // 120
                  IW = 0;
                  for (JE = 1; JE <= M; JE++) { // 110
                     if ( IBLOCK( JE ) != 0 && ( W( JE ) > WKILL || IW == 0 ) ) {
                        IW = JE;
                        WKILL = W( JE );
                     }
                  } // 110
                  IBLOCK[IW] = 0;
               } // 120
            }
            IM = 0;
            for (JE = 1; JE <= M; JE++) { // 130
               if ( IBLOCK( JE ) != 0 ) {
                  IM = IM + 1;
                  W[IM] = W( JE );
                  IBLOCK[IM] = IBLOCK( JE );
               }
            } // 130
            M = IM;
         }
         if ( IDISCL < 0 || IDISCU < 0 ) {
            TOOFEW = true;
         }
      }

      // If ORDER='B', do nothing -- the eigenvalues are already sorted
      //    by block.
      // If ORDER='E', sort the eigenvalues from smallest to largest

      if ( IORDER == 1 && NSPLIT > 1 ) {
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
               ITMP1 = IBLOCK( IE );
               W[IE] = W( JE );
               IBLOCK[IE] = IBLOCK( JE );
               W[JE] = TMP1;
               IBLOCK[JE] = ITMP1;
            }
         } // 150
      }

      INFO = 0;
      if (NCNVRG) INFO = INFO + 1;
      IF( TOOFEW ) INFO = INFO + 2;
      }
