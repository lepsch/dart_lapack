      SUBROUTINE SLARRC( JOBT, N, VL, VU, D, E, PIVMIN, EIGCNT, LCNT, RCNT, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBT;
      int                EIGCNT, INFO, LCNT, N, RCNT;
      REAL               PIVMIN, VL, VU
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      bool               MATT;
      REAL               LPIVOT, RPIVOT, SL, SU, TMP, TMP2

      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      INFO = 0
      LCNT = 0
      RCNT = 0
      EIGCNT = 0

      // Quick return if possible

      if ( N.LE.0 ) {
         RETURN
      }

      MATT = LSAME( JOBT, 'T' )


      if (MATT) {
         // Sturm sequence count on T
         LPIVOT = D( 1 ) - VL
         RPIVOT = D( 1 ) - VU
         if ( LPIVOT.LE.ZERO ) {
            LCNT = LCNT + 1
         ENDIF
         if ( RPIVOT.LE.ZERO ) {
            RCNT = RCNT + 1
         ENDIF
         DO 10 I = 1, N-1
            TMP = E(I)**2
            LPIVOT = ( D( I+1 )-VL ) - TMP/LPIVOT
            RPIVOT = ( D( I+1 )-VU ) - TMP/RPIVOT
            if ( LPIVOT.LE.ZERO ) {
               LCNT = LCNT + 1
            ENDIF
            if ( RPIVOT.LE.ZERO ) {
               RCNT = RCNT + 1
            ENDIF
 10      CONTINUE
      } else {
         // Sturm sequence count on L D L^T
         SL = -VL
         SU = -VU
         DO 20 I = 1, N - 1
            LPIVOT = D( I ) + SL
            RPIVOT = D( I ) + SU
            if ( LPIVOT.LE.ZERO ) {
               LCNT = LCNT + 1
            ENDIF
            if ( RPIVOT.LE.ZERO ) {
               RCNT = RCNT + 1
            ENDIF
            TMP = E(I) * D(I) * E(I)

            TMP2 = TMP / LPIVOT
            if ( TMP2.EQ.ZERO ) {
               SL =  TMP - VL
            } else {
               SL = SL*TMP2 - VL
            }

            TMP2 = TMP / RPIVOT
            if ( TMP2.EQ.ZERO ) {
               SU =  TMP - VU
            } else {
               SU = SU*TMP2 - VU
            }
 20      CONTINUE
         LPIVOT = D( N ) + SL
         RPIVOT = D( N ) + SU
         if ( LPIVOT.LE.ZERO ) {
            LCNT = LCNT + 1
         ENDIF
         if ( RPIVOT.LE.ZERO ) {
            RCNT = RCNT + 1
         ENDIF
      ENDIF
      EIGCNT = RCNT - LCNT

      RETURN

      // End of SLARRC

      }
