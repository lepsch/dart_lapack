      void slaqr5(WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV, WV, LDWV, NH, WH, LDWH ) {
      // IMPLICIT NONE

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, LDWH, LDWV, LDZ, N, NH, NSHFTS, NV;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      REAL               H( LDH, * ), SI( * ), SR( * ), U( LDU, * ), V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * );
      // ..

// ================================================================
      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      REAL               ALPHA, BETA, H11, H12, H21, H22, REFSUM, SAFMAX, SAFMIN, SCL, SMLNUM, SWAP, T1, T2, T3, TST1, TST2, ULP;
      int                I, I2, I4, INCOL, J, JBOT, JCOL, JLEN, JROW, JTOP, K, K1, KDU, KMS, KRCOL, M, M22, MBOT, MTOP, NBMPS, NDCOL, NS, NU;
      bool               ACCUM, BMP22;
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..

      // INTRINSIC ABS, MAX, MIN, MOD, REAL
      // ..
      // .. Local Arrays ..
      REAL               VT( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLAQR1, SLARFG, SLASET, STRMM
      // ..
      // .. Executable Statements ..

      // ==== If there are no shifts, then there is nothing to do. ====

      if (NSHFTS < 2) return;

      // ==== If the active block is empty or 1-by-1, then there
      // .    is nothing to do. ====

      if (KTOP >= KBOT) return;

      // ==== Shuffle shifts into pairs of real shifts and pairs
      // .    of complex conjugate shifts assuming complex
      // .    conjugate shifts are already adjacent to one
      // .    another. ====

      DO 10 I = 1, NSHFTS - 2, 2;
         if ( SI( I ) != -SI( I+1 ) ) {

            SWAP = SR( I );
            SR( I ) = SR( I+1 );
            SR( I+1 ) = SR( I+2 );
            SR( I+2 ) = SWAP;

            SWAP = SI( I );
            SI( I ) = SI( I+1 );
            SI( I+1 ) = SI( I+2 );
            SI( I+2 ) = SWAP;
         }
      } // 10

      // ==== NSHFTS is supposed to be even, but if it is odd,
      // .    then simply reduce it by one.  The shuffle above
      // .    ensures that the dropped shift is real and that
      // .    the remaining shifts are paired. ====

      NS = NSHFTS - MOD( NSHFTS, 2 );

      // ==== Machine constants for deflation ====

      SAFMIN = SLAMCH( 'SAFE MINIMUM' );
      SAFMAX = ONE / SAFMIN;
      ULP = SLAMCH( 'PRECISION' );
      SMLNUM = SAFMIN*( REAL( N ) / ULP );

      // ==== Use accumulated reflections to update far-from-diagonal
      // .    entries ? ====

      ACCUM = ( KACC22 == 1 ) || ( KACC22 == 2 );

      // ==== clear trash ====

      if (KTOP+2 <= KBOT) H( KTOP+2, KTOP ) = ZERO;

      // ==== NBMPS = number of 2-shift bulges in the chain ====

      NBMPS = NS / 2;

      // ==== KDU = width of slab ====

      KDU = 4*NBMPS;

      // ==== Create and chase chains of NBMPS bulges ====

      DO 180 INCOL = KTOP - 2*NBMPS + 1, KBOT - 2, 2*NBMPS;

         // JTOP = Index from which updates from the right start.

         if ( ACCUM ) {
            JTOP = max( KTOP, INCOL );
         } else if ( WANTT ) {
            JTOP = 1;
         } else {
            JTOP = KTOP;
         }

         NDCOL = INCOL + KDU;
         if (ACCUM) slaset( 'ALL', KDU, KDU, ZERO, ONE, U, LDU );

         // ==== Near-the-diagonal bulge chase.  The following loop
         // .    performs the near-the-diagonal part of a small bulge
         // .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
         // .    chunk extends from column INCOL to column NDCOL
         // .    (including both column INCOL and column NDCOL). The
         // .    following loop chases a 2*NBMPS+1 column long chain of
         // .    NBMPS bulges 2*NBMPS-1 columns to the right.  (INCOL
         // .    may be less than KTOP and and NDCOL may be greater than
         // .    KBOT indicating phantom columns from which to chase
         // .    bulges before they are actually introduced or to which
         // .    to chase bulges beyond column KBOT.)  ====

         DO 145 KRCOL = INCOL, min( INCOL+2*NBMPS-1, KBOT-2 );

            // ==== Bulges number MTOP to MBOT are active double implicit
            // .    shift bulges.  There may or may not also be small
            // .    2-by-2 bulge, if there is room.  The inactive bulges
            // .    (if any) must wait until the active bulges have moved
            // .    down the diagonal to make room.  The phantom matrix
            // .    paradigm described above helps keep track.  ====

            MTOP = max( 1, ( KTOP-KRCOL ) / 2+1 );
            MBOT = min( NBMPS, ( KBOT-KRCOL-1 ) / 2 );
            M22 = MBOT + 1;
            BMP22 = ( MBOT < NBMPS ) && ( KRCOL+2*( M22-1 ) ) == ( KBOT-2 );

            // ==== Generate reflections to chase the chain right
            // .    one column.  (The minimum value of K is KTOP-1.) ====

            if ( BMP22 ) {

               // ==== Special case: 2-by-2 reflection at bottom treated
               // .    separately ====

               K = KRCOL + 2*( M22-1 );
               if ( K == KTOP-1 ) {
                  slaqr1(2, H( K+1, K+1 ), LDH, SR( 2*M22-1 ), SI( 2*M22-1 ), SR( 2*M22 ), SI( 2*M22 ), V( 1, M22 ) );
                  BETA = V( 1, M22 );
                  slarfg(2, BETA, V( 2, M22 ), 1, V( 1, M22 ) );
               } else {
                  BETA = H( K+1, K );
                  V( 2, M22 ) = H( K+2, K );
                  slarfg(2, BETA, V( 2, M22 ), 1, V( 1, M22 ) );
                  H( K+1, K ) = BETA;
                  H( K+2, K ) = ZERO;
               }


               // ==== Perform update from right within
               // .    computational window. ====

               T1 = V( 1, M22 );
               T2 = T1*V( 2, M22 );
               DO 30 J = JTOP, min( KBOT, K+3 );
                  REFSUM = H( J, K+1 ) + V( 2, M22 )*H( J, K+2 );
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM*T1;
                  H( J, K+2 ) = H( J, K+2 ) - REFSUM*T2;
               } // 30

               // ==== Perform update from left within
               // .    computational window. ====

               if ( ACCUM ) {
                  JBOT = min( NDCOL, KBOT );
               } else if ( WANTT ) {
                  JBOT = N;
               } else {
                  JBOT = KBOT;
               }
               T1 = V( 1, M22 );
               T2 = T1*V( 2, M22 );
               for (J = K+1; J <= JBOT; J++) { // 40
                  REFSUM = H( K+1, J ) + V( 2, M22 )*H( K+2, J );
                  H( K+1, J ) = H( K+1, J ) - REFSUM*T1;
                  H( K+2, J ) = H( K+2, J ) - REFSUM*T2;
               } // 40

               // ==== The following convergence test requires that
               // .    the tradition small-compared-to-nearby-diagonals
               // .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
               // .    criteria both be satisfied.  The latter improves
               // .    accuracy in some examples. Falling back on an
               // .    alternate convergence criterion when TST1 or TST2
               // .    is zero (as done here) is traditional but probably
               // .    unnecessary. ====

               if ( K >= KTOP ) {
                  if ( H( K+1, K ) != ZERO ) {
                     TST1 = ABS( H( K, K ) ) + ABS( H( K+1, K+1 ) );
                     if ( TST1 == ZERO ) {
                        if (K >= KTOP+1) TST1 = TST1 + ABS( H( K, K-1 ) );
                        if( K >= KTOP+2 ) TST1 = TST1 + ABS( H( K, K-2 ) );
                        if( K >= KTOP+3 ) TST1 = TST1 + ABS( H( K, K-3 ) );
                        if( K <= KBOT-2 ) TST1 = TST1 + ABS( H( K+2, K+1 ) );
                        if( K <= KBOT-3 ) TST1 = TST1 + ABS( H( K+3, K+1 ) );
                        IF( K <= KBOT-4 ) TST1 = TST1 + ABS( H( K+4, K+1 ) );
                     }
                     if( ABS( H( K+1, K ) ) <= max( SMLNUM, ULP*TST1 ) ) {
                        H12 = max( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )                         H21 = min( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )                         H11 = max( ABS( H( K+1, K+1 ) ), ABS( H( K, K )-H( K+1, K+1 ) ) )                         H22 = min( ABS( H( K+1, K+1 ) ), ABS( H( K, K )-H( K+1, K+1 ) ) );
                        SCL = H11 + H12;
                        TST2 = H22*( H11 / SCL );

                        if ( TST2 == ZERO || H21*( H12 / SCL ) <= max( SMLNUM, ULP*TST2 ) ) {
                           H( K+1, K ) = ZERO;
                        }
                     }
                  }
               }

               // ==== Accumulate orthogonal transformations. ====

               if ( ACCUM ) {
                  KMS = K - INCOL;
                  T1 = V( 1, M22 );
                  T2 = T1*V( 2, M22 );
                  DO 50 J = max( 1, KTOP-INCOL ), KDU;
                     REFSUM = U( J, KMS+1 ) + V( 2, M22 )*U( J, KMS+2 );
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM*T1;
                     U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*T2;
                     } // 50
               } else if ( WANTZ ) {
                  T1 = V( 1, M22 );
                  T2 = T1*V( 2, M22 );
                  for (J = ILOZ; J <= IHIZ; J++) { // 60
                     REFSUM = Z( J, K+1 )+V( 2, M22 )*Z( J, K+2 );
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM*T1;
                     Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*T2;
                  } // 60
               }
            }

            // ==== Normal case: Chain of 3-by-3 reflections ====

            DO 80 M = MBOT, MTOP, -1;
               K = KRCOL + 2*( M-1 );
               if ( K == KTOP-1 ) {
                  slaqr1(3, H( KTOP, KTOP ), LDH, SR( 2*M-1 ), SI( 2*M-1 ), SR( 2*M ), SI( 2*M ), V( 1, M ) );
                  ALPHA = V( 1, M );
                  slarfg(3, ALPHA, V( 2, M ), 1, V( 1, M ) );
               } else {

                  // ==== Perform delayed transformation of row below
                  // .    Mth bulge. Exploit fact that first two elements
                  // .    of row are actually zero. ====

                  T1 = V( 1, M );
                  T2 = T1*V( 2, M );
                  T3 = T1*V( 3, M );
                  REFSUM = V( 3, M )*H( K+3, K+2 );
                  H( K+3, K   ) = -REFSUM*T1;
                  H( K+3, K+1 ) = -REFSUM*T2;
                  H( K+3, K+2 ) = H( K+3, K+2 ) - REFSUM*T3;

                  // ==== Calculate reflection to move
                  // .    Mth bulge one step. ====

                  BETA      = H( K+1, K );
                  V( 2, M ) = H( K+2, K );
                  V( 3, M ) = H( K+3, K );
                  slarfg(3, BETA, V( 2, M ), 1, V( 1, M ) );

                  // ==== A Bulge may collapse because of vigilant
                  // .    deflation or destructive underflow.  In the
                  // .    underflow case, try the two-small-subdiagonals
                  // .    trick to try to reinflate the bulge.  ====

                  if ( H( K+3, K ) != ZERO || H( K+3, K+1 ) != ZERO || H( K+3, K+2 ) == ZERO ) {

                     // ==== Typical case: not collapsed (yet). ====

                     H( K+1, K ) = BETA;
                     H( K+2, K ) = ZERO;
                     H( K+3, K ) = ZERO;
                  } else {

                     // ==== Atypical case: collapsed.  Attempt to
                     // .    reintroduce ignoring H(K+1,K) and H(K+2,K).
                     // .    If the fill resulting from the new
                     // .    reflector is too large, then abandon it.
                     // .    Otherwise, use the new one. ====

                     slaqr1(3, H( K+1, K+1 ), LDH, SR( 2*M-1 ), SI( 2*M-1 ), SR( 2*M ), SI( 2*M ), VT );
                     ALPHA = VT( 1 );
                     slarfg(3, ALPHA, VT( 2 ), 1, VT( 1 ) );
                     T1 = VT( 1 );
                     T2 = T1*VT( 2 );
                     T3 = T2*VT( 3 );
                     REFSUM = H( K+1, K )+VT( 2 )*H( K+2, K );

                     if ( ABS( H( K+2, K )-REFSUM*T2 )+ ABS( REFSUM*T3 ) > ULP* ( ABS( H( K, K ) )+ABS( H( K+1, K+1 ) )+ABS( H( K+2, K+2 ) ) ) ) {

                        // ==== Starting a new bulge here would
                        // .    create non-negligible fill.  Use
                        // .    the old one with trepidation. ====

                        H( K+1, K ) = BETA;
                        H( K+2, K ) = ZERO;
                        H( K+3, K ) = ZERO;
                     } else {

                        // ==== Starting a new bulge here would
                        // .    create only negligible fill.
                        // .    Replace the old reflector with
                        // .    the new one. ====

                        H( K+1, K ) = H( K+1, K ) - REFSUM*T1;
                        H( K+2, K ) = ZERO;
                        H( K+3, K ) = ZERO;
                        V( 1, M ) = VT( 1 );
                        V( 2, M ) = VT( 2 );
                        V( 3, M ) = VT( 3 );
                     }
                  }
               }

               // ====  Apply reflection from the right and
               // .     the first column of update from the left.
               // .     These updates are required for the vigilant
               // .     deflation check. We still delay most of the
               // .     updates from the left for efficiency. ====

               T1 = V( 1, M );
               T2 = T1*V( 2, M );
               T3 = T1*V( 3, M );
               DO 70 J = JTOP, min( KBOT, K+3 );
                  REFSUM = H( J, K+1 ) + V( 2, M )*H( J, K+2 ) + V( 3, M )*H( J, K+3 );
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM*T1;
                  H( J, K+2 ) = H( J, K+2 ) - REFSUM*T2;
                  H( J, K+3 ) = H( J, K+3 ) - REFSUM*T3;
               } // 70

               // ==== Perform update from left for subsequent
               // .    column. ====

               REFSUM = H( K+1, K+1 ) + V( 2, M )*H( K+2, K+1 ) + V( 3, M )*H( K+3, K+1 );
               H( K+1, K+1 ) = H( K+1, K+1 ) - REFSUM*T1;
               H( K+2, K+1 ) = H( K+2, K+1 ) - REFSUM*T2;
               H( K+3, K+1 ) = H( K+3, K+1 ) - REFSUM*T3;

               // ==== The following convergence test requires that
               // .    the tradition small-compared-to-nearby-diagonals
               // .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
               // .    criteria both be satisfied.  The latter improves
               // .    accuracy in some examples. Falling back on an
               // .    alternate convergence criterion when TST1 or TST2
               // .    is zero (as done here) is traditional but probably
               // .    unnecessary. ====

               if (K < KTOP) CYCLE;
               if ( H( K+1, K ) != ZERO ) {
                  TST1 = ABS( H( K, K ) ) + ABS( H( K+1, K+1 ) );
                  if ( TST1 == ZERO ) {
                     if (K >= KTOP+1) TST1 = TST1 + ABS( H( K, K-1 ) );
                     if( K >= KTOP+2 ) TST1 = TST1 + ABS( H( K, K-2 ) );
                     if( K >= KTOP+3 ) TST1 = TST1 + ABS( H( K, K-3 ) );
                     if( K <= KBOT-2 ) TST1 = TST1 + ABS( H( K+2, K+1 ) );
                     if( K <= KBOT-3 ) TST1 = TST1 + ABS( H( K+3, K+1 ) );
                     IF( K <= KBOT-4 ) TST1 = TST1 + ABS( H( K+4, K+1 ) );
                  }
                  if ( ABS( H( K+1, K ) ) <= max( SMLNUM, ULP*TST1 ) ) {
                     H12 = max( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) );
                     H21 = min( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) );
                     H11 = max( ABS( H( K+1, K+1 ) ), ABS( H( K, K )-H( K+1, K+1 ) ) )                      H22 = min( ABS( H( K+1, K+1 ) ), ABS( H( K, K )-H( K+1, K+1 ) ) );
                     SCL = H11 + H12;
                     TST2 = H22*( H11 / SCL );

                     if ( TST2 == ZERO || H21*( H12 / SCL ) <= max( SMLNUM, ULP*TST2 ) ) {
                        H( K+1, K ) = ZERO;
                     }
                  }
               }
            } // 80

            // ==== Multiply H by reflections from the left ====

            if ( ACCUM ) {
               JBOT = min( NDCOL, KBOT );
            } else if ( WANTT ) {
               JBOT = N;
            } else {
               JBOT = KBOT;
            }

            DO 100 M = MBOT, MTOP, -1;
               K = KRCOL + 2*( M-1 );
               T1 = V( 1, M );
               T2 = T1*V( 2, M );
               T3 = T1*V( 3, M );
               DO 90 J = max( KTOP, KRCOL + 2*M ), JBOT;
                  REFSUM = H( K+1, J ) + V( 2, M )*H( K+2, J ) + V( 3, M )*H( K+3, J );
                  H( K+1, J ) = H( K+1, J ) - REFSUM*T1;
                  H( K+2, J ) = H( K+2, J ) - REFSUM*T2;
                  H( K+3, J ) = H( K+3, J ) - REFSUM*T3;
               } // 90
            } // 100

            // ==== Accumulate orthogonal transformations. ====

            if ( ACCUM ) {

               // ==== Accumulate U. (If needed, update Z later
               // .    with an efficient matrix-matrix
               // .    multiply.) ====

               DO 120 M = MBOT, MTOP, -1;
                  K = KRCOL + 2*( M-1 );
                  KMS = K - INCOL;
                  I2 = max( 1, KTOP-INCOL );
                  I2 = max( I2, KMS-(KRCOL-INCOL)+1 );
                  I4 = min( KDU, KRCOL + 2*( MBOT-1 ) - INCOL + 5 );
                  T1 = V( 1, M );
                  T2 = T1*V( 2, M );
                  T3 = T1*V( 3, M );
                  for (J = I2; J <= I4; J++) { // 110
                     REFSUM = U( J, KMS+1 ) + V( 2, M )*U( J, KMS+2 ) + V( 3, M )*U( J, KMS+3 );
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM*T1;
                     U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*T2;
                     U( J, KMS+3 ) = U( J, KMS+3 ) - REFSUM*T3;
                  } // 110
               } // 120
            } else if ( WANTZ ) {

               // ==== U is not accumulated, so update Z
               // .    now by multiplying by reflections
               // .    from the right. ====

               DO 140 M = MBOT, MTOP, -1;
                  K = KRCOL + 2*( M-1 );
                  T1 = V( 1, M );
                  T2 = T1*V( 2, M );
                  T3 = T1*V( 3, M );
                  for (J = ILOZ; J <= IHIZ; J++) { // 130
                     REFSUM = Z( J, K+1 ) + V( 2, M )*Z( J, K+2 ) + V( 3, M )*Z( J, K+3 );
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM*T1;
                     Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*T2;
                     Z( J, K+3 ) = Z( J, K+3 ) - REFSUM*T3;
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
               JTOP = 1;
               JBOT = N;
            } else {
               JTOP = KTOP;
               JBOT = KBOT;
            }
            K1 = max( 1, KTOP-INCOL );
            NU = ( KDU-max( 0, NDCOL-KBOT ) ) - K1 + 1;

            // ==== Horizontal Multiply ====

            DO 150 JCOL = min( NDCOL, KBOT ) + 1, JBOT, NH;
               JLEN = min( NH, JBOT-JCOL+1 );
               sgemm('C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ), LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH, LDWH );
               slacpy('ALL', NU, JLEN, WH, LDWH, H( INCOL+K1, JCOL ), LDH );
            } // 150

            // ==== Vertical multiply ====

            DO 160 JROW = JTOP, max( KTOP, INCOL ) - 1, NV;
               JLEN = min( NV, max( KTOP, INCOL )-JROW );
               sgemm('N', 'N', JLEN, NU, NU, ONE, H( JROW, INCOL+K1 ), LDH, U( K1, K1 ), LDU, ZERO, WV, LDWV );
               slacpy('ALL', JLEN, NU, WV, LDWV, H( JROW, INCOL+K1 ), LDH );
            } // 160

            // ==== Z multiply (also vertical) ====

            if ( WANTZ ) {
               DO 170 JROW = ILOZ, IHIZ, NV;
                  JLEN = min( NV, IHIZ-JROW+1 );
                  sgemm('N', 'N', JLEN, NU, NU, ONE, Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ), LDU, ZERO, WV, LDWV );
                  slacpy('ALL', JLEN, NU, WV, LDWV, Z( JROW, INCOL+K1 ), LDZ );
               } // 170
            }
         }
      } // 180

      // ==== End of SLAQR5 ====

      }
