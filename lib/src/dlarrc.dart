      void dlarrc(JOBT, N, VL, VU, D, E, PIVMIN, EIGCNT, LCNT, RCNT, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBT;
      int                EIGCNT, INFO, LCNT, N, RCNT;
      double             PIVMIN, VL, VU;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      bool               MATT;
      double             LPIVOT, RPIVOT, SL, SU, TMP, TMP2;

      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Executable Statements ..

      INFO = 0;
      LCNT = 0;
      RCNT = 0;
      EIGCNT = 0;

      // Quick return if possible

      if ( N <= 0 ) {
         return;
      }

      MATT = lsame( JOBT, 'T' );


      if (MATT) {
         // Sturm sequence count on T
         LPIVOT = D( 1 ) - VL;
         RPIVOT = D( 1 ) - VU;
         if ( LPIVOT <= ZERO ) {
            LCNT = LCNT + 1;
         }
         if ( RPIVOT <= ZERO ) {
            RCNT = RCNT + 1;
         }
         for (I = 1; I <= N-1; I++) { // 10
            TMP = E(I)**2;
            LPIVOT = ( D( I+1 )-VL ) - TMP/LPIVOT;
            RPIVOT = ( D( I+1 )-VU ) - TMP/RPIVOT;
            if ( LPIVOT <= ZERO ) {
               LCNT = LCNT + 1;
            }
            if ( RPIVOT <= ZERO ) {
               RCNT = RCNT + 1;
            }
         } // 10
      } else {
         // Sturm sequence count on L D L^T
         SL = -VL;
         SU = -VU;
         for (I = 1; I <= N - 1; I++) { // 20
            LPIVOT = D( I ) + SL;
            RPIVOT = D( I ) + SU;
            if ( LPIVOT <= ZERO ) {
               LCNT = LCNT + 1;
            }
            if ( RPIVOT <= ZERO ) {
               RCNT = RCNT + 1;
            }
            TMP = E(I) * D(I) * E(I);

            TMP2 = TMP / LPIVOT;
            if ( TMP2 == ZERO ) {
               SL =  TMP - VL;
            } else {
               SL = SL*TMP2 - VL;
            }

            TMP2 = TMP / RPIVOT;
            if ( TMP2 == ZERO ) {
               SU =  TMP - VU;
            } else {
               SU = SU*TMP2 - VU;
            }
         } // 20
         LPIVOT = D( N ) + SL;
         RPIVOT = D( N ) + SU;
         if ( LPIVOT <= ZERO ) {
            LCNT = LCNT + 1;
         }
         if ( RPIVOT <= ZERO ) {
            RCNT = RCNT + 1;
         }
      }
      EIGCNT = RCNT - LCNT;

      return;
      }
