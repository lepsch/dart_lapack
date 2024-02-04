      void cggbak(JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V, LDV, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOB, SIDE;
      int                IHI, ILO, INFO, LDV, M, N;
      // ..
      // .. Array Arguments ..
      double               LSCALE( * ), RSCALE( * );
      Complex            V( LDV, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LEFTV, RIGHTV;
      int                I, K;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSSCAL, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      RIGHTV = lsame( SIDE, 'R' );
      LEFTV = lsame( SIDE, 'L' );

      INFO = 0;
      if ( !lsame( JOB, 'N' ) && !lsame( JOB, 'P' ) && !lsame( JOB, 'S' ) && !lsame( JOB, 'B' ) ) {
         INFO = -1;
      } else if ( !RIGHTV && !LEFTV ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ILO < 1 ) {
         INFO = -4;
      } else if ( N == 0 && IHI == 0 && ILO != 1 ) {
         INFO = -4;
      } else if ( N > 0 && ( IHI < ILO || IHI > max( 1, N ) ) ) {
         INFO = -5;
      } else if ( N == 0 && ILO == 1 && IHI != 0 ) {
         INFO = -5;
      } else if ( M < 0 ) {
         INFO = -8;
      } else if ( LDV < max( 1, N ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('CGGBAK', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;
      if( M == 0 ) return;
      IF( lsame( JOB, 'N' ) ) return;

      if (ILO == IHI) GO TO 30;

      // Backward balance

      if ( lsame( JOB, 'S' ) || lsame( JOB, 'B' ) ) {

         // Backward transformation on right eigenvectors

         if ( RIGHTV ) {
            for (I = ILO; I <= IHI; I++) { // 10
               csscal(M, RSCALE( I ), V( I, 1 ), LDV );
            } // 10
         }

         // Backward transformation on left eigenvectors

         if ( LEFTV ) {
            for (I = ILO; I <= IHI; I++) { // 20
               csscal(M, LSCALE( I ), V( I, 1 ), LDV );
            } // 20
         }
      }

      // Backward permutation

      } // 30
      if ( lsame( JOB, 'P' ) || lsame( JOB, 'B' ) ) {

         // Backward permutation on right eigenvectors

         if ( RIGHTV ) {
            if (ILO == 1) GO TO 50;
            for (I = ILO - 1; I >= 1; I--) { // 40
               K = INT( RSCALE( I ) );
               if (K == I) GO TO 40;
               cswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 40

            } // 50
            if (IHI == N) GO TO 70;
            for (I = IHI + 1; I <= N; I++) { // 60
               K = INT( RSCALE( I ) );
               if (K == I) GO TO 60;
               cswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 60
         }

         // Backward permutation on left eigenvectors

         } // 70
         if ( LEFTV ) {
            if (ILO == 1) GO TO 90;
            for (I = ILO - 1; I >= 1; I--) { // 80
               K = INT( LSCALE( I ) );
               if (K == I) GO TO 80;
               cswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 80

            } // 90
            if (IHI == N) GO TO 110;
            for (I = IHI + 1; I <= N; I++) { // 100
               K = INT( LSCALE( I ) );
               if (K == I) GO TO 100;
               cswap(M, V( I, 1 ), LDV, V( K, 1 ), LDV );
            } // 100
         }
      }

      } // 110

      return;
      }
