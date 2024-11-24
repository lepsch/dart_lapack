// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slaqr3(final int WANTT, final int WANTZ, final int N, final int KTOP, final int KBOT, final int NW, final Matrix<double> H_, final int LDH, final int ILOZ, final int IHIZ, final Matrix<double> Z_, final int LDZ, final int NS, final int ND, final int SR, final int SI, final Matrix<double> V_, final int LDV, final int NH, final Matrix<double> T_, final int LDT, final int NV, final Matrix<double> WV_, final int LDWV, final Array<double> WORK_, final int LWORK,) {
  final H = H_.dim();
  final Z = Z_.dim();
  final V = V_.dim();
  final T = T_.dim();
  final WV = WV_.dim();
  final WORK = WORK_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, LDZ, LWORK, N, ND, NH, NS, NV, NW;
      bool               WANTT, WANTZ;
      double               H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), V( LDV, * ), WORK( * ), WV( LDWV, * ), Z( LDZ, * );
      // ..

// ================================================================
      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP;
      int                I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3, LWKOPT, NMIN;
      bool               BULGE, SORTED;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SROUNDUP_LWORK;
      //- int                ILAENV;
      // EXTERNAL SLAMCH, SROUNDUP_LWORK, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEHRD, SGEMM, SLACPY, SLAHQR, SLANV2, SLAQR4, SLARF, SLARFG, SLASET, SORMHR, STREXC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, MAX, MIN, REAL, SQRT

      // ==== Estimate optimal workspace. ====

      JW = min( NW, KBOT-KTOP+1 );
      if ( JW <= 2 ) {
         LWKOPT = 1;
      } else {

         // ==== Workspace query call to SGEHRD ====

         sgehrd(JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO );
         LWK1 = INT( WORK( 1 ) );

         // ==== Workspace query call to SORMHR ====

         sormhr('R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, WORK, -1, INFO );
         LWK2 = INT( WORK( 1 ) );

         // ==== Workspace query call to SLAQR4 ====

         slaqr4( true , true , JW, 1, JW, T, LDT, SR, SI, 1, JW, V, LDV, WORK, -1, INFQR );
         LWK3 = INT( WORK( 1 ) );

         // ==== Optimal workspace ====

         LWKOPT = max( JW+max( LWK1, LWK2 ), LWK3 );
      }

      // ==== Quick return in case of workspace query. ====

      if ( LWORK == -1 ) {
         WORK[1] = SROUNDUP_LWORK( LWKOPT );
         return;
      }

      // ==== Nothing to do ...
      // ... for an empty active block ... ====
      NS = 0;
      ND = 0;
      WORK[1] = ONE;
      if (KTOP > KBOT) return;
      // ... nor for an empty deflation window. ====
      if (NW < 1) return;

      // ==== Machine constants ====

      SAFMIN = SLAMCH( 'SAFE MINIMUM' );
      SAFMAX = ONE / SAFMIN;
      ULP = SLAMCH( 'PRECISION' );
      SMLNUM = SAFMIN*( double( N ) / ULP );

      // ==== Setup deflation window ====

      JW = min( NW, KBOT-KTOP+1 );
      KWTOP = KBOT - JW + 1;
      if ( KWTOP == KTOP ) {
         S = ZERO;
      } else {
         S = H( KWTOP, KWTOP-1 );
      }

      if ( KBOT == KWTOP ) {

         // ==== 1-by-1 deflation window: not much to do ====

         SR[KWTOP] = H( KWTOP, KWTOP );
         SI[KWTOP] = ZERO;
         NS = 1;
         ND = 0;
         if ( ( S ).abs() <= max( SMLNUM, ULP*( H( KWTOP, KWTOP ) ).abs() ) ) {
            NS = 0;
            ND = 1;
            if (KWTOP > KTOP) H( KWTOP, KWTOP-1 ) = ZERO;
         }
         WORK[1] = ONE;
         return;
      }

      // ==== Convert to spike-triangular form.  (In case of a
      // .    rare QR failure, this routine continues to do
      // .    aggressive early deflation using that part of
      // .    the deflation window that converged using INFQR
      // .    here and there to keep track.) ====

      slacpy('U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT );
      scopy(JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 );

      slaset('A', JW, JW, ZERO, ONE, V, LDV );
      NMIN = ilaenv( 12, 'SLAQR3', 'SV', JW, 1, JW, LWORK );
      if ( JW > NMIN ) {
         slaqr4( true , true , JW, 1, JW, T, LDT, SR( KWTOP ), SI( KWTOP ), 1, JW, V, LDV, WORK, LWORK, INFQR );
      } else {
         slahqr( true , true , JW, 1, JW, T, LDT, SR( KWTOP ), SI( KWTOP ), 1, JW, V, LDV, INFQR );
      }

      // ==== STREXC needs a clean margin near the diagonal ====

      for (J = 1; J <= JW - 3; J++) { // 10
         T[J+2][J] = ZERO;
         T[J+3][J] = ZERO;
      } // 10
      if (JW > 2) T( JW, JW-2 ) = ZERO;

      // ==== Deflation detection loop ====

      NS = JW;
      ILST = INFQR + 1;
      } // 20
      if ( ILST <= NS ) {
         if ( NS == 1 ) {
            BULGE = false;
         } else {
            BULGE = T( NS, NS-1 ) != ZERO;
         }

         // ==== Small spike tip test for deflation ====

         if ( !BULGE ) {

            // ==== Real eigenvalue ====

            FOO = ( T( NS, NS ) ).abs();
            if (FOO == ZERO) FOO = ( S ).abs();
            if ( ( S*V( 1, NS ) ).abs() <= max( SMLNUM, ULP*FOO ) ) {

               // ==== Deflatable ====

               NS = NS - 1;
            } else {

               // ==== Undeflatable.   Move it up out of the way.
               // .    (STREXC can not fail in this case.) ====

               IFST = NS;
               strexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               ILST = ILST + 1;
            }
         } else {

            // ==== Complex conjugate pair ====

            FOO = ( T( NS, NS ) ).abs() + sqrt( ( T( NS, NS-1 ) ).abs() )* sqrt( ( T( NS-1, NS ) ).abs() )             IF( FOO == ZERO ) FOO = ( S ).abs()             IF( max( ( S*V( 1, NS ) ).abs(), ( S*V( 1, NS-1 ) ).abs() ) <= max( SMLNUM, ULP*FOO ) ) THEN;

               // ==== Deflatable ====

               NS = NS - 2;
            } else {

               // ==== Undeflatable. Move them up out of the way.
               // .    Fortunately, STREXC does the right thing with
               // .    ILST in case of a rare exchange failure. ====

               IFST = NS;
               strexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               ILST = ILST + 2;
            }
         }

         // ==== End deflation detection loop ====

         GO TO 20;
      }

         // ==== Return to Hessenberg form ====

      if (NS == 0) S = ZERO;

      if ( NS < JW ) {

         // ==== sorting diagonal blocks of T improves accuracy for
         // .    graded matrices.  Bubble sort deals well with
         // .    exchange failures. ====

         SORTED = false;
         I = NS + 1;
         } // 30
         if (SORTED) GO TO 50;
         SORTED = true;

         KEND = I - 1;
         I = INFQR + 1;
         if ( I == NS ) {
            K = I + 1;
         } else if ( T( I+1, I ) == ZERO ) {
            K = I + 1;
         } else {
            K = I + 2;
         }
         } // 40
         if ( K <= KEND ) {
            if ( K == I+1 ) {
               EVI = ( T( I, I ) ).abs();
            } else {
               EVI = ( T( I, I ) ).abs() + sqrt( ( T( I+1, I ) ).abs() )* sqrt( ( T( I, I+1 ) ).abs() );
            }

            if ( K == KEND ) {
               EVK = ( T( K, K ) ).abs();
            } else if ( T( K+1, K ) == ZERO ) {
               EVK = ( T( K, K ) ).abs();
            } else {
               EVK = ( T( K, K ) ).abs() + sqrt( ( T( K+1, K ) ).abs() )* sqrt( ( T( K, K+1 ) ).abs() );
            }

            if ( EVI >= EVK ) {
               I = K;
            } else {
               SORTED = false;
               IFST = I;
               ILST = K;
               strexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               if ( INFO == 0 ) {
                  I = ILST;
               } else {
                  I = K;
               }
            }
            if ( I == KEND ) {
               K = I + 1;
            } else if ( T( I+1, I ) == ZERO ) {
               K = I + 1;
            } else {
               K = I + 2;
            }
            GO TO 40;
         }
         GO TO 30;
         } // 50
      }

      // ==== Restore shift/eigenvalue array from T ====

      I = JW;
      } // 60
      if ( I >= INFQR+1 ) {
         if ( I == INFQR+1 ) {
            SR[KWTOP+I-1] = T( I, I );
            SI[KWTOP+I-1] = ZERO;
            I = I - 1;
         } else if ( T( I, I-1 ) == ZERO ) {
            SR[KWTOP+I-1] = T( I, I );
            SI[KWTOP+I-1] = ZERO;
            I = I - 1;
         } else {
            AA = T( I-1, I-1 );
            CC = T( I, I-1 );
            BB = T( I-1, I );
            DD = T( I, I );
            slanv2(AA, BB, CC, DD, SR( KWTOP+I-2 ), SI( KWTOP+I-2 ), SR( KWTOP+I-1 ), SI( KWTOP+I-1 ), CS, SN );
            I = I - 2;
         }
         GO TO 60;
      }

      if ( NS < JW || S == ZERO ) {
         if ( NS > 1 && S != ZERO ) {

            // ==== Reflect spike back into lower triangle ====

            scopy(NS, V, LDV, WORK, 1 );
            BETA = WORK( 1 );
            slarfg(NS, BETA, WORK( 2 ), 1, TAU );
            WORK[1] = ONE;

            slaset('L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT );

            slarf('L', NS, JW, WORK, 1, TAU, T, LDT, WORK( JW+1 ) );
            slarf('R', NS, NS, WORK, 1, TAU, T, LDT, WORK( JW+1 ) );
            slarf('R', JW, NS, WORK, 1, TAU, V, LDV, WORK( JW+1 ) );

            sgehrd(JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), LWORK-JW, INFO );
         }

         // ==== Copy updated reduced window into place ====

         if (KWTOP > 1) H( KWTOP, KWTOP-1 ) = S*V( 1, 1 );
         slacpy('U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH );
         scopy(JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), LDH+1 );

         // ==== Accumulate orthogonal matrix in order update
         // .    H and Z, if requested.  ====

         if (NS > 1 && S != ZERO) sormhr( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, WORK( JW+1 ), LWORK-JW, INFO );

         // ==== Update vertical slab in H ====

         if ( WANTT ) {
            LTOP = 1;
         } else {
            LTOP = KTOP;
         }
         for (KROW = LTOP; NV < 0 ? KROW >= KWTOP - 1 : KROW <= KWTOP - 1; KROW += NV) { // 70
            KLN = min( NV, KWTOP-KROW );
            sgemm('N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), LDH, V, LDV, ZERO, WV, LDWV );
            slacpy('A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH );
         } // 70

         // ==== Update horizontal slab in H ====

         if ( WANTT ) {
            for (KCOL = KBOT + 1; NH < 0 ? KCOL >= N : KCOL <= N; KCOL += NH) { // 80
               KLN = min( NH, N-KCOL+1 );
               sgemm('C', 'N', JW, KLN, JW, ONE, V, LDV, H( KWTOP, KCOL ), LDH, ZERO, T, LDT );
               slacpy('A', JW, KLN, T, LDT, H( KWTOP, KCOL ), LDH );
            } // 80
         }

         // ==== Update vertical slab in Z ====

         if ( WANTZ ) {
            for (KROW = ILOZ; NV < 0 ? KROW >= IHIZ : KROW <= IHIZ; KROW += NV) { // 90
               KLN = min( NV, IHIZ-KROW+1 );
               sgemm('N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), LDZ, V, LDV, ZERO, WV, LDWV );
               slacpy('A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), LDZ );
            } // 90
         }
      }

      // ==== Return the number of deflations ... ====

      ND = JW - NS;

      // ==== ... and the number of shifts. (Subtracting
      // .    INFQR from the spike length takes care
      // .    of the case of a rare QR failure while
      // .    calculating eigenvalues of the deflation
      // .    window.)  ====

      NS = NS - INFQR;

       // ==== Return optimal workspace. ====

      WORK[1] = SROUNDUP_LWORK( LWKOPT );

      // ==== End of SLAQR3 ====

      }
