      void sgebak(JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOB, SIDE;
      int                IHI, ILO, INFO, LDV, M, N;
      double               V( LDV, * ), SCALE( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      bool               LEFTV, RIGHTV;
      int                I, II, K;
      double               S;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Decode and Test the input parameters

      RIGHTV = lsame( SIDE, 'R' );
      LEFTV = lsame( SIDE, 'L' );

      INFO = 0;
      if ( !lsame( JOB, 'N' ) && !lsame( JOB, 'P' ) && !lsame( JOB, 'S' ) && !lsame( JOB, 'B' ) ) {
         INFO = -1;
      } else if ( !RIGHTV && !LEFTV ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ILO < 1 || ILO > max( 1, N ) ) {
         INFO = -4;
      } else if ( IHI < min( ILO, N ) || IHI > N ) {
         INFO = -5;
      } else if ( M < 0 ) {
         INFO = -7;
      } else if ( LDV < max( 1, N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('SGEBAK', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;
      if( M == 0 ) return;
      IF( lsame( JOB, 'N' ) ) return;

      if (ILO == IHI) GO TO 30;

      // Backward balance

      if ( lsame( JOB, 'S' ) || lsame( JOB, 'B' ) ) {

         if ( RIGHTV ) {
            for (I = ILO; I <= IHI; I++) { // 10
               S = SCALE( I );
               sscal(M, S, V( I, 1 ), LDV );
            } // 10
         }

         if ( LEFTV ) {
            for (I = ILO; I <= IHI; I++) { // 20
               S = ONE / SCALE( I );
               sscal(M, S, V( I, 1 ), LDV );
            } // 20
         }

      }

      // Backward permutation

      // For  I = ILO-1 step -1 until 1,
               // IHI+1 step 1 until N do --

      } // 30
      if ( lsame( JOB, 'P' ) || lsame( JOB, 'B' ) ) {
         if ( RIGHTV ) {
            for (II = 1; II <= N; II++) { // 40
               I = II;
               if (I >= ILO && I <= IHI) GO TO 40;
               IF( I < ILO ) I = ILO - II;
               K = INT( SCALE( I ) );
               if (K == I) GO TO 40;
               sswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 40
         }

         if ( LEFTV ) {
            for (II = 1; II <= N; II++) { // 50
               I = II;
               if (I >= ILO && I <= IHI) GO TO 50;
               IF( I < ILO ) I = ILO - II;
               K = INT( SCALE( I ) );
               if (K == I) GO TO 50;
               sswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 50
         }
      }

      }
