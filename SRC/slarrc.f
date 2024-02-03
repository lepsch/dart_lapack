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

      IF( N.LE.0 ) THEN
         RETURN
      END IF

      MATT = LSAME( JOBT, 'T' )


      IF (MATT) THEN
         // Sturm sequence count on T
         LPIVOT = D( 1 ) - VL
         RPIVOT = D( 1 ) - VU
         IF( LPIVOT.LE.ZERO ) THEN
            LCNT = LCNT + 1
         ENDIF
         IF( RPIVOT.LE.ZERO ) THEN
            RCNT = RCNT + 1
         ENDIF
         DO 10 I = 1, N-1
            TMP = E(I)**2
            LPIVOT = ( D( I+1 )-VL ) - TMP/LPIVOT
            RPIVOT = ( D( I+1 )-VU ) - TMP/RPIVOT
            IF( LPIVOT.LE.ZERO ) THEN
               LCNT = LCNT + 1
            ENDIF
            IF( RPIVOT.LE.ZERO ) THEN
               RCNT = RCNT + 1
            ENDIF
 10      CONTINUE
      ELSE
         // Sturm sequence count on L D L^T
         SL = -VL
         SU = -VU
         DO 20 I = 1, N - 1
            LPIVOT = D( I ) + SL
            RPIVOT = D( I ) + SU
            IF( LPIVOT.LE.ZERO ) THEN
               LCNT = LCNT + 1
            ENDIF
            IF( RPIVOT.LE.ZERO ) THEN
               RCNT = RCNT + 1
            ENDIF
            TMP = E(I) * D(I) * E(I)

            TMP2 = TMP / LPIVOT
            IF( TMP2.EQ.ZERO ) THEN
               SL =  TMP - VL
            ELSE
               SL = SL*TMP2 - VL
            END IF

            TMP2 = TMP / RPIVOT
            IF( TMP2.EQ.ZERO ) THEN
               SU =  TMP - VU
            ELSE
               SU = SU*TMP2 - VU
            END IF
 20      CONTINUE
         LPIVOT = D( N ) + SL
         RPIVOT = D( N ) + SU
         IF( LPIVOT.LE.ZERO ) THEN
            LCNT = LCNT + 1
         ENDIF
         IF( RPIVOT.LE.ZERO ) THEN
            RCNT = RCNT + 1
         ENDIF
      ENDIF
      EIGCNT = RCNT - LCNT

      RETURN

      // End of SLARRC

      }
