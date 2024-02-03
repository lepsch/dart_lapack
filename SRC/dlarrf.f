      SUBROUTINE DLARRF( N, D, L, LD, CLSTRT, CLEND, W, WGAP, WERR, SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA, DPLUS, LPLUS, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                CLSTRT, CLEND, INFO, N;
      double             CLGAPL, CLGAPR, PIVMIN, SIGMA, SPDIAM;
      // ..
      // .. Array Arguments ..
      double             D( * ), DPLUS( * ), L( * ), LD( * ), LPLUS( * ), W( * ), WGAP( * ), WERR( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             FOUR, MAXGROWTH1, MAXGROWTH2, ONE, QUART, TWO;
      const              ONE = 1.0D0, TWO = 2.0D0, FOUR = 4.0D0, QUART = 0.25D0, MAXGROWTH1 = 8.D0, MAXGROWTH2 = 8.D0 ;
      // ..
      // .. Local Scalars ..
      bool      DORRR1, FORCER, NOFAIL, SAWNAN1, SAWNAN2, TRYRRR1;
      int                I, INDX, KTRY, KTRYMAX, SLEFT, SRIGHT, SHIFT;
      const              KTRYMAX = 1, SLEFT = 1, SRIGHT = 2 ;
      double             AVGAP, BESTSHIFT, CLWDTH, EPS, FACT, FAIL, FAIL2, GROWTHBOUND, LDELTA, LDMAX, LSIGMA, MAX1, MAX2, MINGAP, OLDP, PROD, RDELTA, RDMAX, RRR1, RRR2, RSIGMA, S, SMLGROWTH, TMP, ZNM2;
      // ..
      // .. External Functions ..
      bool    DISNAN;
      double             DLAMCH;
      // EXTERNAL DISNAN, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Quick return if possible

      if ( N.LE.0 ) {
         RETURN
      }

      FACT = DBLE(2**KTRYMAX)
      EPS = DLAMCH( 'Precision' )
      SHIFT = 0
      FORCER = false;


      // Note that we cannot guarantee that for any of the shifts tried,
      // the factorization has a small or even moderate element growth.
      // There could be Ritz values at both ends of the cluster and despite
      // backing off, there are examples where all factorizations tried
      // (in IEEE mode, allowing zero pivots & infinities) have INFINITE
      // element growth.
      // For this reason, we should use PIVMIN in this subroutine so that at
      // least the L D L^T factorization exists. It can be checked afterwards
      // whether the element growth caused bad residuals/orthogonality.

      // Decide whether the code should accept the best among all
      // representations despite large element growth or signal INFO=1
      // Setting NOFAIL to false for quick fix for bug 113
      NOFAIL = false;


      // Compute the average gap length of the cluster
      CLWDTH = ABS(W(CLEND)-W(CLSTRT)) + WERR(CLEND) + WERR(CLSTRT)
      AVGAP = CLWDTH / DBLE(CLEND-CLSTRT)
      MINGAP = MIN(CLGAPL, CLGAPR)
      // Initial values for shifts to both ends of cluster
      LSIGMA = MIN(W( CLSTRT ),W( CLEND )) - WERR( CLSTRT )
      RSIGMA = MAX(W( CLSTRT ),W( CLEND )) + WERR( CLEND )

      // Use a small fudge to make sure that we really shift to the outside
      LSIGMA = LSIGMA - ABS(LSIGMA)* FOUR * EPS
      RSIGMA = RSIGMA + ABS(RSIGMA)* FOUR * EPS

      // Compute upper bounds for how much to back off the initial shifts
      LDMAX = QUART * MINGAP + TWO * PIVMIN
      RDMAX = QUART * MINGAP + TWO * PIVMIN

      LDELTA = MAX(AVGAP,WGAP( CLSTRT ))/FACT
      RDELTA = MAX(AVGAP,WGAP( CLEND-1 ))/FACT

      // Initialize the record of the best representation found

      S = DLAMCH( 'S' )
      SMLGROWTH = ONE / S
      FAIL = DBLE(N-1)*MINGAP/(SPDIAM*EPS)
      FAIL2 = DBLE(N-1)*MINGAP/(SPDIAM*SQRT(EPS))
      BESTSHIFT = LSIGMA

      // while (KTRY <= KTRYMAX)
      KTRY = 0
      GROWTHBOUND = MAXGROWTH1*SPDIAM

      } // 5
      SAWNAN1 = false;
      SAWNAN2 = false;
      // Ensure that we do not back off too much of the initial shifts
      LDELTA = MIN(LDMAX,LDELTA)
      RDELTA = MIN(RDMAX,RDELTA)

      // Compute the element growth when shifting to both ends of the cluster
      // accept the shift if there is no element growth at one of the two ends

      // Left end
      S = -LSIGMA
      DPLUS( 1 ) = D( 1 ) + S
      if (ABS(DPLUS(1)).LT.PIVMIN) {
         DPLUS(1) = -PIVMIN
         // Need to set SAWNAN1 because refined RRR test should not be used
         // in this case
         SAWNAN1 = true;
      }
      MAX1 = ABS( DPLUS( 1 ) )
      for (I = 1; I <= N - 1; I++) { // 6
         LPLUS( I ) = LD( I ) / DPLUS( I )
         S = S*LPLUS( I )*L( I ) - LSIGMA
         DPLUS( I+1 ) = D( I+1 ) + S
         if (ABS(DPLUS(I+1)).LT.PIVMIN) {
            DPLUS(I+1) = -PIVMIN
            // Need to set SAWNAN1 because refined RRR test should not be used
            // in this case
            SAWNAN1 = true;
         }
         MAX1 = MAX( MAX1,ABS(DPLUS(I+1)) )
      } // 6
      SAWNAN1 = SAWNAN1 .OR.  DISNAN( MAX1 )
       if ( FORCER .OR. (MAX1.LE.GROWTHBOUND && .NOT.SAWNAN1 ) ) {
         SIGMA = LSIGMA
         SHIFT = SLEFT
         GOTO 100
      }

      // Right end
      S = -RSIGMA
      WORK( 1 ) = D( 1 ) + S
      if (ABS(WORK(1)).LT.PIVMIN) {
         WORK(1) = -PIVMIN
         // Need to set SAWNAN2 because refined RRR test should not be used
         // in this case
         SAWNAN2 = true;
      }
      MAX2 = ABS( WORK( 1 ) )
      for (I = 1; I <= N - 1; I++) { // 7
         WORK( N+I ) = LD( I ) / WORK( I )
         S = S*WORK( N+I )*L( I ) - RSIGMA
         WORK( I+1 ) = D( I+1 ) + S
         if (ABS(WORK(I+1)).LT.PIVMIN) {
            WORK(I+1) = -PIVMIN
            // Need to set SAWNAN2 because refined RRR test should not be used
            // in this case
            SAWNAN2 = true;
         }
         MAX2 = MAX( MAX2,ABS(WORK(I+1)) )
      } // 7
      SAWNAN2 = SAWNAN2 .OR.  DISNAN( MAX2 )
       if ( FORCER .OR. (MAX2.LE.GROWTHBOUND && .NOT.SAWNAN2 ) ) {
         SIGMA = RSIGMA
         SHIFT = SRIGHT
         GOTO 100
      }
      // If we are at this point, both shifts led to too much element growth

      // Record the better of the two shifts (provided it didn't lead to NaN)
      if (SAWNAN1 && SAWNAN2) {
         // both MAX1 and MAX2 are NaN
         GOTO 50
      } else {
         if ( .NOT.SAWNAN1 ) {
            INDX = 1
            if (MAX1.LE.SMLGROWTH) {
               SMLGROWTH = MAX1
               BESTSHIFT = LSIGMA
            }
         }
         if ( .NOT.SAWNAN2 ) {
            if (SAWNAN1 .OR. MAX2.LE.MAX1) INDX = 2;
            if (MAX2.LE.SMLGROWTH) {
               SMLGROWTH = MAX2
               BESTSHIFT = RSIGMA
            }
         }
      }

      // If we are here, both the left and the right shift led to
      // element growth. If the element growth is moderate, then
      // we may still accept the representation, if it passes a
      // refined test for RRR. This test supposes that no NaN occurred.
      // Moreover, we use the refined RRR test only for isolated clusters.
      if ((CLWDTH.LT.MINGAP/DBLE(128)) && (MIN(MAX1,MAX2).LT.FAIL2) && (.NOT.SAWNAN1) && (.NOT.SAWNAN2)) {
         DORRR1 = true;
      } else {
         DORRR1 = false;
      }
      TRYRRR1 = true;
      if ( TRYRRR1 && DORRR1 ) {
      if (INDX == 1) {
         TMP = ABS( DPLUS( N ) )
         ZNM2 = ONE
         PROD = ONE
         OLDP = ONE
         DO 15 I = N-1, 1, -1
            if ( PROD .LE. EPS ) {
               PROD = ((DPLUS(I+1)*WORK(N+I+1))/(DPLUS(I)*WORK(N+I)))*OLDP
            } else {
               PROD = PROD*ABS(WORK(N+I))
            }
            OLDP = PROD
            ZNM2 = ZNM2 + PROD**2
            TMP = MAX( TMP, ABS( DPLUS( I ) * PROD ))
         } // 15
         RRR1 = TMP/( SPDIAM * SQRT( ZNM2 ) )
         if (RRR1.LE.MAXGROWTH2) {
            SIGMA = LSIGMA
            SHIFT = SLEFT
            GOTO 100
         }
      } else if (INDX == 2) {
         TMP = ABS( WORK( N ) )
         ZNM2 = ONE
         PROD = ONE
         OLDP = ONE
         DO 16 I = N-1, 1, -1
            if ( PROD .LE. EPS ) {
               PROD = ((WORK(I+1)*LPLUS(I+1))/(WORK(I)*LPLUS(I)))*OLDP
            } else {
               PROD = PROD*ABS(LPLUS(I))
            }
            OLDP = PROD
            ZNM2 = ZNM2 + PROD**2
            TMP = MAX( TMP, ABS( WORK( I ) * PROD ))
         } // 16
         RRR2 = TMP/( SPDIAM * SQRT( ZNM2 ) )
         if (RRR2.LE.MAXGROWTH2) {
            SIGMA = RSIGMA
            SHIFT = SRIGHT
            GOTO 100
         }
      }
      }

      } // 50

      if (KTRY.LT.KTRYMAX) {
         // If we are here, both shifts failed also the RRR test.
         // Back off to the outside
         LSIGMA = MAX( LSIGMA - LDELTA, LSIGMA - LDMAX)          RSIGMA = MIN( RSIGMA + RDELTA, RSIGMA + RDMAX )
         LDELTA = TWO * LDELTA
         RDELTA = TWO * RDELTA
         KTRY = KTRY + 1
         GOTO 5
      } else {
         // None of the representations investigated satisfied our
         // criteria. Take the best one we found.
         if ((SMLGROWTH.LT.FAIL).OR.NOFAIL) {
            LSIGMA = BESTSHIFT
            RSIGMA = BESTSHIFT
            FORCER = true;
            GOTO 5
         } else {
            INFO = 1
            RETURN
         }
      }

      } // 100
      if (SHIFT == SLEFT) {
      } else if (SHIFT == SRIGHT) {
         // store new L and D back into DPLUS, LPLUS
         dcopy(N, WORK, 1, DPLUS, 1 );
         dcopy(N-1, WORK(N+1), 1, LPLUS, 1 );
      }

      RETURN

      // End of DLARRF

      }
