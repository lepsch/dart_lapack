      void slarrk(final int N, final int IW, final int GL, final int GU, final int D, final int E2, final int PIVMIN, final int RELTOL, final int W, final int WERR, final Box<int> INFO,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int       INFO, IW, N;
      double                PIVMIN, RELTOL, GL, GU, W, WERR;
      double               D( * ), E2( * );
      // ..

      double               FUDGE, HALF, TWO, ZERO;
      const              HALF = 0.5, TWO = 2.0, FUDGE = TWO, ZERO = 0.0 ;
      int       I, IT, ITMAX, NEGCNT;
      double               ATOLI, EPS, LEFT, MID, RIGHT, RTOLI, TMP1, TMP2, TNORM;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX

      // Quick return if possible

      if ( N <= 0 ) {
         INFO = 0;
         return;
      }

      // Get machine constants
      EPS = SLAMCH( 'P' );

      TNORM = max( ( GL ).abs(), ( GU ).abs() );
      RTOLI = RELTOL;
      ATOLI = FUDGE*TWO*PIVMIN;
       ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2;

      INFO = -1;

      LEFT = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN;
      RIGHT = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN;
      IT = 0;

      } // 10

      // Check if interval converged or maximum number of iterations reached

      TMP1 = ( RIGHT - LEFT ).abs();
      TMP2 = max( (RIGHT).abs(), (LEFT).abs() );
      if ( TMP1 < max( ATOLI, PIVMIN, RTOLI*TMP2 ) ) {
         INFO = 0;
         GOTO 30;
      }
      if (IT > ITMAX) GOTO 30;


      // Count number of negative pivots for mid-point

      IT = IT + 1;
      MID = HALF * (LEFT + RIGHT);
      NEGCNT = 0;
      TMP1 = D( 1 ) - MID;
      if( ( TMP1 ).abs() < PIVMIN ) TMP1 = -PIVMIN;
      IF( TMP1 <= ZERO ) NEGCNT = NEGCNT + 1;

      for (I = 2; I <= N; I++) { // 20
         TMP1 = D( I ) - E2( I-1 ) / TMP1 - MID;
         if( ( TMP1 ).abs() < PIVMIN ) TMP1 = -PIVMIN;
         IF( TMP1 <= ZERO ) NEGCNT = NEGCNT + 1;
      } // 20

      if (NEGCNT >= IW) {
         RIGHT = MID;
      } else {
         LEFT = MID;
      }
      GOTO 10;

      } // 30

      // Converged or maximum number of iterations reached

      W = HALF * (LEFT + RIGHT);
      WERR = HALF * ( RIGHT - LEFT ).abs();

      }
