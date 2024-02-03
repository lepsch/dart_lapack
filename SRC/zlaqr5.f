      SUBROUTINE ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, WV, LDWV, NH, WH, LDWH )
      IMPLICIT NONE

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, LDWH, LDWV, LDZ, N, NH, NSHFTS, NV;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
      // ..

*  ================================================================
      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0d0, 0.0d0 ), ONE = ( 1.0d0, 0.0d0 ) ;
      double             RZERO, RONE;
      const              RZERO = 0.0d0, RONE = 1.0d0 ;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         ALPHA, BETA, CDUM, REFSUM, T1, T2, T3
      double             H11, H12, H21, H22, SAFMAX, SAFMIN, SCL, SMLNUM, TST1, TST2, ULP;
      int                I2, I4, INCOL, J, JBOT, JCOL, JLEN, JROW, JTOP, K, K1, KDU, KMS, KRCOL, M, M22, MBOT, MTOP, NBMPS, NDCOL, NS, NU;
      bool               ACCUM, BMP22;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..

      // INTRINSIC ABS, DBLE, DCONJG, DIMAG, MAX, MIN, MOD
      // ..
      // .. Local Arrays ..
      COMPLEX*16         VT( 3 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZLACPY, ZLAQR1, ZLARFG, ZLASET, ZTRMM
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // ==== If there are no shifts, then there is nothing to do. ====

      if (NSHFTS.LT.2) RETURN;

      // ==== If the active block is empty or 1-by-1, then there
      // .    is nothing to do. ====

      if (KTOP.GE.KBOT) RETURN;

      // ==== NSHFTS is supposed to be even, but if it is odd,
      // .    then simply reduce it by one.  ====

      NS = NSHFTS - MOD( NSHFTS, 2 )

      // ==== Machine constants for deflation ====

      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )

      // ==== Use accumulated reflections to update far-from-diagonal
      // .    entries ? ====

      ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )

      // ==== clear trash ====

      if (KTOP+2.LE.KBOT) H( KTOP+2, KTOP ) = ZERO;

      // ==== NBMPS = number of 2-shift bulges in the chain ====

      NBMPS = NS / 2

      // ==== KDU = width of slab ====

      KDU = 4*NBMPS

      // ==== Create and chase chains of NBMPS bulges ====

      DO 180 INCOL = KTOP - 2*NBMPS + 1, KBOT - 2, 2*NBMPS

         // JTOP = Index from which updates from the right start.

         if ( ACCUM ) {
            JTOP = MAX( KTOP, INCOL )
         } else if ( WANTT ) {
            JTOP = 1
         } else {
            JTOP = KTOP
         }

         NDCOL = INCOL + KDU
         if (ACCUM) CALL ZLASET( 'ALL', KDU, KDU, ZERO, ONE, U, LDU );

         // ==== Near-the-diagonal bulge chase.  The following loop
         // .    performs the near-the-diagonal part of a small bulge
         // .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
         // .    chunk extends from column INCOL to column NDCOL
         // .    (including both column INCOL and column NDCOL). The
         // .    following loop chases a 2*NBMPS+1 column long chain of
         // .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
         // .    may be less than KTOP and and NDCOL may be greater than
         // .    KBOT indicating phantom columns from which to chase
         // .    bulges before they are actually introduced or to which
         // .    to chase bulges beyond column KBOT.)  ====

         DO 145 KRCOL = INCOL, MIN( INCOL+2*NBMPS-1, KBOT-2 )

            // ==== Bulges number MTOP to MBOT are active double implicit
            // .    shift bulges.  There may or may not also be small
            // .    2-by-2 bulge, if there is room.  The inactive bulges
            // .    (if any) must wait until the active bulges have moved
            // .    down the diagonal to make room.  The phantom matrix
            // .    paradigm described above helps keep track.  ====

            MTOP = MAX( 1, ( KTOP-KRCOL ) / 2+1 )
            MBOT = MIN( NBMPS, ( KBOT-KRCOL-1 ) / 2 )
            M22 = MBOT + 1
            BMP22 = ( MBOT.LT.NBMPS ) .AND. ( KRCOL+2*( M22-1 ) ).EQ. ( KBOT-2 )

            // ==== Generate reflections to chase the chain right
            // .    one column.  (The minimum value of K is KTOP-1.) ====

            if ( BMP22 ) {

               // ==== Special case: 2-by-2 reflection at bottom treated
               // .    separately ====

               K = KRCOL + 2*( M22-1 )
               if ( K.EQ.KTOP-1 ) {
                  zlaqr1(2, H( K+1, K+1 ), LDH, S( 2*M22-1 ), S( 2*M22 ), V( 1, M22 ) );
                  BETA = V( 1, M22 )
                  zlarfg(2, BETA, V( 2, M22 ), 1, V( 1, M22 ) );
               } else {
                  BETA = H( K+1, K )
                  V( 2, M22 ) = H( K+2, K )
                  zlarfg(2, BETA, V( 2, M22 ), 1, V( 1, M22 ) );
                  H( K+1, K ) = BETA
                  H( K+2, K ) = ZERO
               }


               // ==== Perform update from right within
               // .    computational window. ====

               T1 = V( 1, M22 )
               T2 = T1*DCONJG( V( 2, M22 ) )
               DO 30 J = JTOP, MIN( KBOT, K+3 )
                  REFSUM = H( J, K+1 ) + V( 2, M22 )*H( J, K+2 )
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM*T1
                  H( J, K+2 ) = H( J, K+2 ) - REFSUM*T2
               } // 30

               // ==== Perform update from left within
               // .    computational window. ====

               if ( ACCUM ) {
                  JBOT = MIN( NDCOL, KBOT )
               } else if ( WANTT ) {
                  JBOT = N
               } else {
                  JBOT = KBOT
               }
               T1 = DCONJG( V( 1, M22 ) )
               T2 = T1*V( 2, M22 )
               for (J = K+1; J <= JBOT; J++) { // 40
                  REFSUM = H( K+1, J ) + DCONJG( V( 2, M22 ) )*H( K+2, J )
                  H( K+1, J ) = H( K+1, J ) - REFSUM*T1
                  H( K+2, J ) = H( K+2, J ) - REFSUM*T2
               } // 40

               // ==== The following convergence test requires that
               // .    the tradition small-compared-to-nearby-diagonals
               // .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
               // .    criteria both be satisfied.  The latter improves
               // .    accuracy in some examples. Falling back on an
               // .    alternate convergence criterion when TST1 or TST2
               // .    is zero (as done here) is traditional but probably
               // .    unnecessary. ====

               if ( K.GE.KTOP ) {
                  if ( H( K+1, K ).NE.ZERO ) {
                     TST1 = CABS1( H( K, K ) ) + CABS1( H( K+1, K+1 ) )
                     if ( TST1.EQ.RZERO ) {
                        if (K.GE.KTOP+1) TST1 = TST1 + CABS1( H( K, K-1 ) )                         IF( K.GE.KTOP+2 ) TST1 = TST1 + CABS1( H( K, K-2 ) )                         IF( K.GE.KTOP+3 ) TST1 = TST1 + CABS1( H( K, K-3 ) )                         IF( K.LE.KBOT-2 ) TST1 = TST1 + CABS1( H( K+2, K+1 ) )                         IF( K.LE.KBOT-3 ) TST1 = TST1 + CABS1( H( K+3, K+1 ) )                         IF( K.LE.KBOT-4 ) TST1 = TST1 + CABS1( H( K+4, K+1 ) );
                     }
                     IF( CABS1( H( K+1, K ) ) .LE.MAX( SMLNUM, ULP*TST1 ) ) THEN                         H12 = MAX( CABS1( H( K+1, K ) ), CABS1( H( K, K+1 ) ) )                         H21 = MIN( CABS1( H( K+1, K ) ), CABS1( H( K, K+1 ) ) )                         H11 = MAX( CABS1( H( K+1, K+1 ) ), CABS1( H( K, K )-H( K+1, K+1 ) ) )                         H22 = MIN( CABS1( H( K+1, K+1 ) ), CABS1( H( K, K )-H( K+1, K+1 ) ) )
                        SCL = H11 + H12
                        TST2 = H22*( H11 / SCL )

                        IF( TST2.EQ.RZERO .OR. H21*( H12 / SCL ).LE. MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
                     }
                  }
               }

               // ==== Accumulate orthogonal transformations. ====

               if ( ACCUM ) {
                  KMS = K - INCOL
                  DO 50 J = MAX( 1, KTOP-INCOL ), KDU
                     REFSUM = V( 1, M22 )*( U( J, KMS+1 )+ V( 2, M22 )*U( J, KMS+2 ) )
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                     U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*DCONJG( V( 2, M22 ) )
                     } // 50
               } else if ( WANTZ ) {
                  for (J = ILOZ; J <= IHIZ; J++) { // 60
                     REFSUM = V( 1, M22 )*( Z( J, K+1 )+V( 2, M22 )* Z( J, K+2 ) )
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                     Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*DCONJG( V( 2, M22 ) )
                  } // 60
               }
            }

            // ==== Normal case: Chain of 3-by-3 reflections ====

            DO 80 M = MBOT, MTOP, -1
               K = KRCOL + 2*( M-1 )
               if ( K.EQ.KTOP-1 ) {
                  zlaqr1(3, H( KTOP, KTOP ), LDH, S( 2*M-1 ), S( 2*M ), V( 1, M ) );
                  ALPHA = V( 1, M )
                  zlarfg(3, ALPHA, V( 2, M ), 1, V( 1, M ) );
               } else {

                  // ==== Perform delayed transformation of row below
                  // .    Mth bulge. Exploit fact that first two elements
                  // .    of row are actually zero. ====

                  T1 = V( 1, M )
                  T2 = T1*DCONJG( V( 2, M ) )
                  T3 = T1*DCONJG( V( 3, M ) )
                  REFSUM = V( 3, M )*H( K+3, K+2 )
                  H( K+3, K   ) = -REFSUM*T1
                  H( K+3, K+1 ) = -REFSUM*T2
                  H( K+3, K+2 ) = H( K+3, K+2 ) - REFSUM*T3

                  // ==== Calculate reflection to move
                  // .    Mth bulge one step. ====

                  BETA      = H( K+1, K )
                  V( 2, M ) = H( K+2, K )
                  V( 3, M ) = H( K+3, K )
                  zlarfg(3, BETA, V( 2, M ), 1, V( 1, M ) );

                  // ==== A Bulge may collapse because of vigilant
                  // .    deflation or destructive underflow.  In the
                  // .    underflow case, try the two-small-subdiagonals
                  // .    trick to try to reinflate the bulge.  ====

                  if ( H( K+3, K ).NE.ZERO .OR. H( K+3, K+1 ).NE. ZERO .OR. H( K+3, K+2 ).EQ.ZERO ) {

                     // ==== Typical case: not collapsed (yet). ====

                     H( K+1, K ) = BETA
                     H( K+2, K ) = ZERO
                     H( K+3, K ) = ZERO
                  } else {

                     // ==== Atypical case: collapsed.  Attempt to
                     // .    reintroduce ignoring H(K+1,K) and H(K+2,K).
                     // .    If the fill resulting from the new
                     // .    reflector is too large, then abandon it.
                     // .    Otherwise, use the new one. ====

                     zlaqr1(3, H( K+1, K+1 ), LDH, S( 2*M-1 ), S( 2*M ), VT );
                     ALPHA = VT( 1 )
                     zlarfg(3, ALPHA, VT( 2 ), 1, VT( 1 ) );
                     T1 = DCONJG( VT( 1 ) )
                     T2 = T1*VT( 2 )
                     T3 = T1*VT( 3 )
                     REFSUM = H( K+1, K )+DCONJG( VT( 2 ) )*H( K+2, K )

                     if ( CABS1( H( K+2, K )-REFSUM*T2 )+ CABS1( REFSUM*T3 ).GT.ULP* ( CABS1( H( K, K ) )+CABS1( H( K+1, K+1 ) )+CABS1( H( K+2, K+2 ) ) ) ) {

                        // ==== Starting a new bulge here would
                        // .    create non-negligible fill.  Use
                        // .    the old one with trepidation. ====

                        H( K+1, K ) = BETA
                        H( K+2, K ) = ZERO
                        H( K+3, K ) = ZERO
                     } else {

                        // ==== Starting a new bulge here would
                        // .    create only negligible fill.
                        // .    Replace the old reflector with
                        // .    the new one. ====

                        H( K+1, K ) = H( K+1, K ) - REFSUM*T1
                        H( K+2, K ) = ZERO
                        H( K+3, K ) = ZERO
                        V( 1, M ) = VT( 1 )
                        V( 2, M ) = VT( 2 )
                        V( 3, M ) = VT( 3 )
                     }
                  }
               }

               // ====  Apply reflection from the right and
               // .     the first column of update from the left.
               // .     These updates are required for the vigilant
               // .     deflation check. We still delay most of the
               // .     updates from the left for efficiency. ====

               T1 = V( 1, M )
               T2 = T1*DCONJG( V( 2, M ) )
               T3 = T1*DCONJG( V( 3, M ) )
               DO 70 J = JTOP, MIN( KBOT, K+3 )
                  REFSUM = H( J, K+1 ) + V( 2, M )*H( J, K+2 ) + V( 3, M )*H( J, K+3 )
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM*T1
                  H( J, K+2 ) = H( J, K+2 ) - REFSUM*T2
                  H( J, K+3 ) = H( J, K+3 ) - REFSUM*T3
               } // 70

               // ==== Perform update from left for subsequent
               // .    column. ====

               T1 = DCONJG( V( 1, M ) )
               T2 = T1*V( 2, M )
               T3 = T1*V( 3, M )
               REFSUM = H( K+1, K+1 ) + DCONJG( V( 2, M ) )*H( K+2, K+1 ) + DCONJG( V( 3, M ) )*H( K+3, K+1 )
               H( K+1, K+1 ) = H( K+1, K+1 ) - REFSUM*T1
               H( K+2, K+1 ) = H( K+2, K+1 ) - REFSUM*T2
               H( K+3, K+1 ) = H( K+3, K+1 ) - REFSUM*T3

               // ==== The following convergence test requires that
               // .    the tradition small-compared-to-nearby-diagonals
               // .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
               // .    criteria both be satisfied.  The latter improves
               // .    accuracy in some examples. Falling back on an
               // .    alternate convergence criterion when TST1 or TST2
               // .    is zero (as done here) is traditional but probably
               // .    unnecessary. ====

               if (K.LT.KTOP) CYCLE;
               if ( H( K+1, K ).NE.ZERO ) {
                  TST1 = CABS1( H( K, K ) ) + CABS1( H( K+1, K+1 ) )
                  if ( TST1.EQ.RZERO ) {
                     if (K.GE.KTOP+1) TST1 = TST1 + CABS1( H( K, K-1 ) )                      IF( K.GE.KTOP+2 ) TST1 = TST1 + CABS1( H( K, K-2 ) )                      IF( K.GE.KTOP+3 ) TST1 = TST1 + CABS1( H( K, K-3 ) )                      IF( K.LE.KBOT-2 ) TST1 = TST1 + CABS1( H( K+2, K+1 ) )                      IF( K.LE.KBOT-3 ) TST1 = TST1 + CABS1( H( K+3, K+1 ) )                      IF( K.LE.KBOT-4 ) TST1 = TST1 + CABS1( H( K+4, K+1 ) );
                  }
                  IF( CABS1( H( K+1, K ) ).LE.MAX( SMLNUM, ULP*TST1 ) ) THEN                      H12 = MAX( CABS1( H( K+1, K ) ), CABS1( H( K, K+1 ) ) )                      H21 = MIN( CABS1( H( K+1, K ) ), CABS1( H( K, K+1 ) ) )                      H11 = MAX( CABS1( H( K+1, K+1 ) ), CABS1( H( K, K )-H( K+1, K+1 ) ) )                      H22 = MIN( CABS1( H( K+1, K+1 ) ), CABS1( H( K, K )-H( K+1, K+1 ) ) )
                     SCL = H11 + H12
                     TST2 = H22*( H11 / SCL )

                     IF( TST2.EQ.RZERO .OR. H21*( H12 / SCL ).LE. MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
                  }
               }
            } // 80

            // ==== Multiply H by reflections from the left ====

            if ( ACCUM ) {
               JBOT = MIN( NDCOL, KBOT )
            } else if ( WANTT ) {
               JBOT = N
            } else {
               JBOT = KBOT
            }

            DO 100 M = MBOT, MTOP, -1
               K = KRCOL + 2*( M-1 )
               T1 = DCONJG( V( 1, M ) )
               T2 = T1*V( 2, M )
               T3 = T1*V( 3, M )
               DO 90 J = MAX( KTOP, KRCOL + 2*M ), JBOT
                  REFSUM = H( K+1, J ) + DCONJG( V( 2, M ) )*H( K+2, J ) + DCONJG( V( 3, M ) )*H( K+3, J )
                  H( K+1, J ) = H( K+1, J ) - REFSUM*T1
                  H( K+2, J ) = H( K+2, J ) - REFSUM*T2
                  H( K+3, J ) = H( K+3, J ) - REFSUM*T3
               } // 90
            } // 100

            // ==== Accumulate orthogonal transformations. ====

            if ( ACCUM ) {

               // ==== Accumulate U. (If needed, update Z later
               // .    with an efficient matrix-matrix
               // .    multiply.) ====

               DO 120 M = MBOT, MTOP, -1
                  K = KRCOL + 2*( M-1 )
                  KMS = K - INCOL
                  I2 = MAX( 1, KTOP-INCOL )
                  I2 = MAX( I2, KMS-(KRCOL-INCOL)+1 )
                  I4 = MIN( KDU, KRCOL + 2*( MBOT-1 ) - INCOL + 5 )
                  T1 = V( 1, M )
                  T2 = T1*DCONJG( V( 2, M ) )
                  T3 = T1*DCONJG( V( 3, M ) )
                  for (J = I2; J <= I4; J++) { // 110
                     REFSUM = U( J, KMS+1 ) + V( 2, M )*U( J, KMS+2 ) + V( 3, M )*U( J, KMS+3 )
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM*T1
                     U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*T2
                     U( J, KMS+3 ) = U( J, KMS+3 ) - REFSUM*T3
                  } // 110
               } // 120
            } else if ( WANTZ ) {

               // ==== U is not accumulated, so update Z
               // .    now by multiplying by reflections
               // .    from the right. ====

               DO 140 M = MBOT, MTOP, -1
                  K = KRCOL + 2*( M-1 )
                  T1 = V( 1, M )
                  T2 = T1*DCONJG( V( 2, M ) )
                  T3 = T1*DCONJG( V( 3, M ) )
                  for (J = ILOZ; J <= IHIZ; J++) { // 130
                     REFSUM = Z( J, K+1 ) + V( 2, M )*Z( J, K+2 ) + V( 3, M )*Z( J, K+3 )
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM*T1
                     Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*T2
                     Z( J, K+3 ) = Z( J, K+3 ) - REFSUM*T3
                  } // 130
               } // 140
            }

            // ==== End of near-the-diagonal bulge chase. ====

         } // 145

         // ==== Use U (if accumulated) to update far-from-diagonal
         // .    entries in H.  If required, use U to update Z as
         // .    well. ====

         if ( ACCUM ) {
            if ( WANTT ) {
               JTOP = 1
               JBOT = N
            } else {
               JTOP = KTOP
               JBOT = KBOT
            }
            K1 = MAX( 1, KTOP-INCOL )
            NU = ( KDU-MAX( 0, NDCOL-KBOT ) ) - K1 + 1

            // ==== Horizontal Multiply ====

            DO 150 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
               JLEN = MIN( NH, JBOT-JCOL+1 )
               zgemm('C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ), LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH, LDWH );
               zlacpy('ALL', NU, JLEN, WH, LDWH, H( INCOL+K1, JCOL ), LDH );
            } // 150

            // ==== Vertical multiply ====

            DO 160 JROW = JTOP, MAX( KTOP, INCOL ) - 1, NV
               JLEN = MIN( NV, MAX( KTOP, INCOL )-JROW )
               zgemm('N', 'N', JLEN, NU, NU, ONE, H( JROW, INCOL+K1 ), LDH, U( K1, K1 ), LDU, ZERO, WV, LDWV );
               zlacpy('ALL', JLEN, NU, WV, LDWV, H( JROW, INCOL+K1 ), LDH );
            } // 160

            // ==== Z multiply (also vertical) ====

            if ( WANTZ ) {
               DO 170 JROW = ILOZ, IHIZ, NV
                  JLEN = MIN( NV, IHIZ-JROW+1 )
                  zgemm('N', 'N', JLEN, NU, NU, ONE, Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ), LDU, ZERO, WV, LDWV );
                  zlacpy('ALL', JLEN, NU, WV, LDWV, Z( JROW, INCOL+K1 ), LDZ );
               } // 170
            }
         }
      } // 180

      // ==== End of ZLAQR5 ====

      }
