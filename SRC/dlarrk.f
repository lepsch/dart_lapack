      void dlarrk(N, IW, GL, GU, D, E2, PIVMIN, RELTOL, W, WERR, INFO) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int       INFO, IW, N;
      double              PIVMIN, RELTOL, GL, GU, W, WERR;
      // ..
      // .. Array Arguments ..
      double             D( * ), E2( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             FUDGE, HALF, TWO, ZERO;
      const              HALF = 0.5, TWO = 2.0, FUDGE = TWO, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int       I, IT, ITMAX, NEGCNT;
      double             ATOLI, EPS, LEFT, MID, RIGHT, RTOLI, TMP1, TMP2, TNORM;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 0 ) {
         INFO = 0;
         return;
      }

      // Get machine constants
      EPS = DLAMCH( 'P' );

      TNORM = max( ABS( GL ), ABS( GU ) );
      RTOLI = RELTOL;
      ATOLI = FUDGE*TWO*PIVMIN;
       ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2;

      INFO = -1;

      LEFT = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN;
      RIGHT = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN;
      IT = 0;

      } // 10

      // Check if interval converged or maximum number of iterations reached

      TMP1 = ABS( RIGHT - LEFT );
      TMP2 = max( ABS(RIGHT), ABS(LEFT) );
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
      if( ABS( TMP1 ) < PIVMIN ) TMP1 = -PIVMIN;
      IF( TMP1 <= ZERO ) NEGCNT = NEGCNT + 1;

      for (I = 2; I <= N; I++) { // 20
         TMP1 = D( I ) - E2( I-1 ) / TMP1 - MID;
         if( ABS( TMP1 ) < PIVMIN ) TMP1 = -PIVMIN;
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
      WERR = HALF * ABS( RIGHT - LEFT );

      return;
      }
