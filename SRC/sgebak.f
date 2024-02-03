      SUBROUTINE SGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB, SIDE;
      int                IHI, ILO, INFO, LDV, M, N;
      // ..
      // .. Array Arguments ..
      REAL               V( LDV, * ), SCALE( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LEFTV, RIGHTV;
      int                I, II, K;
      REAL               S
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and Test the input parameters

      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )

      INFO = 0
      if ( .NOT.LSAME( JOB, 'N' ) && .NOT.LSAME( JOB, 'P' ) && .NOT.LSAME( JOB, 'S' ) && .NOT.LSAME( JOB, 'B' ) ) {
         INFO = -1
      } else if ( .NOT.RIGHTV && .NOT.LEFTV ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( ILO < 1 || ILO > MAX( 1, N ) ) {
         INFO = -4
      } else if ( IHI < MIN( ILO, N ) || IHI > N ) {
         INFO = -5
      } else if ( M < 0 ) {
         INFO = -7
      } else if ( LDV < MAX( 1, N ) ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('SGEBAK', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN       IF( M == 0 ) RETURN       IF( LSAME( JOB, 'N' ) ) RETURN;

      if (ILO == IHI) GO TO 30;

      // Backward balance

      if ( LSAME( JOB, 'S' ) || LSAME( JOB, 'B' ) ) {

         if ( RIGHTV ) {
            for (I = ILO; I <= IHI; I++) { // 10
               S = SCALE( I )
               sscal(M, S, V( I, 1 ), LDV );
            } // 10
         }

         if ( LEFTV ) {
            for (I = ILO; I <= IHI; I++) { // 20
               S = ONE / SCALE( I )
               sscal(M, S, V( I, 1 ), LDV );
            } // 20
         }

      }

      // Backward permutation

      // For  I = ILO-1 step -1 until 1,
               // IHI+1 step 1 until N do --

      } // 30
      if ( LSAME( JOB, 'P' ) || LSAME( JOB, 'B' ) ) {
         if ( RIGHTV ) {
            for (II = 1; II <= N; II++) { // 40
               I = II
               if (I.GE.ILO && I.LE.IHI) GO TO 40                IF( I < ILO ) I = ILO - II;
               K = INT( SCALE( I ) )
               if (K == I) GO TO 40;
               sswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 40
         }

         if ( LEFTV ) {
            for (II = 1; II <= N; II++) { // 50
               I = II
               if (I.GE.ILO && I.LE.IHI) GO TO 50                IF( I < ILO ) I = ILO - II;
               K = INT( SCALE( I ) )
               if (K == I) GO TO 50;
               sswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 50
         }
      }

      RETURN

      // End of SGEBAK

      }
