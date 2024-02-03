      SUBROUTINE SLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, LDZ, LWORK, N, ND, NH, NS, NV, NW;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      REAL               H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), V( LDV, * ), WORK( * ), WV( LDWV, * ), Z( LDZ, * )
      // ..

*  ================================================================
      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0e0, ONE = 1.0e0 ;
      // ..
      // .. Local Scalars ..
      REAL               AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP;
      int                I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWKOPT;
      bool               BULGE, SORTED;
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SROUNDUP_LWORK
      // EXTERNAL SLAMCH, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEHRD, SGEMM, SLACPY, SLAHQR, SLANV2, SLARF, SLARFG, SLASET, SORMHR, STREXC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // ==== Estimate optimal workspace. ====

      JW = MIN( NW, KBOT-KTOP+1 )
      if ( JW.LE.2 ) {
         LWKOPT = 1
      } else {

         // ==== Workspace query call to SGEHRD ====

         sgehrd(JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO );
         LWK1 = INT( WORK( 1 ) )

         // ==== Workspace query call to SORMHR ====

         sormhr('R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, WORK, -1, INFO );
         LWK2 = INT( WORK( 1 ) )

         // ==== Optimal workspace ====

         LWKOPT = JW + MAX( LWK1, LWK2 )
      }

      // ==== Quick return in case of workspace query. ====

      if ( LWORK.EQ.-1 ) {
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         RETURN
      }

      // ==== Nothing to do ...
      // ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      if (KTOP.GT.KBOT) RETURN;
      // ... nor for an empty deflation window. ====
      if (NW.LT.1) RETURN;

      // ==== Machine constants ====

      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( N ) / ULP )

      // ==== Setup deflation window ====

      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      if ( KWTOP.EQ.KTOP ) {
         S = ZERO
      } else {
         S = H( KWTOP, KWTOP-1 )
      }

      if ( KBOT.EQ.KWTOP ) {

         // ==== 1-by-1 deflation window: not much to do ====

         SR( KWTOP ) = H( KWTOP, KWTOP )
         SI( KWTOP ) = ZERO
         NS = 1
         ND = 0
         if ( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( H( KWTOP, KWTOP ) ) ) ) {
            NS = 0
            ND = 1
            if (KWTOP.GT.KTOP) H( KWTOP, KWTOP-1 ) = ZERO;
         }
         WORK( 1 ) = ONE
         RETURN
      }

      // ==== Convert to spike-triangular form.  (In case of a
      // .    rare QR failure, this routine continues to do
      // .    aggressive early deflation using that part of
      // .    the deflation window that converged using INFQR
      // .    here and there to keep track.) ====

      slacpy('U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT );
      scopy(JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 );

      slaset('A', JW, JW, ZERO, ONE, V, LDV );
      slahqr(.true., .true., JW, 1, JW, T, LDT, SR( KWTOP ), SI( KWTOP ), 1, JW, V, LDV, INFQR );

      // ==== STREXC needs a clean margin near the diagonal ====

      for (J = 1; J <= JW - 3; J++) { // 10
         T( J+2, J ) = ZERO
         T( J+3, J ) = ZERO
      } // 10
      if (JW.GT.2) T( JW, JW-2 ) = ZERO;

      // ==== Deflation detection loop ====

      NS = JW
      ILST = INFQR + 1
      } // 20
      if ( ILST.LE.NS ) {
         if ( NS.EQ.1 ) {
            BULGE = .FALSE.
         } else {
            BULGE = T( NS, NS-1 ).NE.ZERO
         }

         // ==== Small spike tip test for deflation ====

         if ( .NOT.BULGE ) {

            // ==== Real eigenvalue ====

            FOO = ABS( T( NS, NS ) )
            if (FOO.EQ.ZERO) FOO = ABS( S );
            if ( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) {

               // ==== Deflatable ====

               NS = NS - 1
            } else {

               // ==== Undeflatable.   Move it up out of the way.
               // .    (STREXC can not fail in this case.) ====

               IFST = NS
               strexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               ILST = ILST + 1
            }
         } else {

            // ==== Complex conjugate pair ====

            FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )* SQRT( ABS( T( NS-1, NS ) ) )             IF( FOO.EQ.ZERO ) FOO = ABS( S )             IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE. MAX( SMLNUM, ULP*FOO ) ) THEN

               // ==== Deflatable ====

               NS = NS - 2
            } else {

               // ==== Undeflatable. Move them up out of the way.
               // .    Fortunately, STREXC does the right thing with
               // .    ILST in case of a rare exchange failure. ====

               IFST = NS
               strexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               ILST = ILST + 2
            }
         }

         // ==== End deflation detection loop ====

         GO TO 20
      }

         // ==== Return to Hessenberg form ====

      if (NS.EQ.0) S = ZERO;

      if ( NS.LT.JW ) {

         // ==== sorting diagonal blocks of T improves accuracy for
         // .    graded matrices.  Bubble sort deals well with
         // .    exchange failures. ====

         SORTED = .false.
         I = NS + 1
         } // 30
         if (SORTED) GO TO 50;
         SORTED = .true.

         KEND = I - 1
         I = INFQR + 1
         if ( I.EQ.NS ) {
            K = I + 1
         } else if ( T( I+1, I ).EQ.ZERO ) {
            K = I + 1
         } else {
            K = I + 2
         }
         } // 40
         if ( K.LE.KEND ) {
            if ( K.EQ.I+1 ) {
               EVI = ABS( T( I, I ) )
            } else {
               EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )* SQRT( ABS( T( I, I+1 ) ) )
            }

            if ( K.EQ.KEND ) {
               EVK = ABS( T( K, K ) )
            } else if ( T( K+1, K ).EQ.ZERO ) {
               EVK = ABS( T( K, K ) )
            } else {
               EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )* SQRT( ABS( T( K, K+1 ) ) )
            }

            if ( EVI.GE.EVK ) {
               I = K
            } else {
               SORTED = .false.
               IFST = I
               ILST = K
               strexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               if ( INFO.EQ.0 ) {
                  I = ILST
               } else {
                  I = K
               }
            }
            if ( I.EQ.KEND ) {
               K = I + 1
            } else if ( T( I+1, I ).EQ.ZERO ) {
               K = I + 1
            } else {
               K = I + 2
            }
            GO TO 40
         }
         GO TO 30
         } // 50
      }

      // ==== Restore shift/eigenvalue array from T ====

      I = JW
      } // 60
      if ( I.GE.INFQR+1 ) {
         if ( I.EQ.INFQR+1 ) {
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         } else if ( T( I, I-1 ).EQ.ZERO ) {
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         } else {
            AA = T( I-1, I-1 )
            CC = T( I, I-1 )
            BB = T( I-1, I )
            DD = T( I, I )
            slanv2(AA, BB, CC, DD, SR( KWTOP+I-2 ), SI( KWTOP+I-2 ), SR( KWTOP+I-1 ), SI( KWTOP+I-1 ), CS, SN );
            I = I - 2
         }
         GO TO 60
      }

      if ( NS.LT.JW .OR. S.EQ.ZERO ) {
         if ( NS.GT.1 .AND. S.NE.ZERO ) {

            // ==== Reflect spike back into lower triangle ====

            scopy(NS, V, LDV, WORK, 1 );
            BETA = WORK( 1 )
            slarfg(NS, BETA, WORK( 2 ), 1, TAU );
            WORK( 1 ) = ONE

            slaset('L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT );

            slarf('L', NS, JW, WORK, 1, TAU, T, LDT, WORK( JW+1 ) );
            slarf('R', NS, NS, WORK, 1, TAU, T, LDT, WORK( JW+1 ) );
            slarf('R', JW, NS, WORK, 1, TAU, V, LDV, WORK( JW+1 ) );

            sgehrd(JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), LWORK-JW, INFO );
         }

         // ==== Copy updated reduced window into place ====

         if (KWTOP.GT.1) H( KWTOP, KWTOP-1 ) = S*V( 1, 1 );
         slacpy('U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH );
         scopy(JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), LDH+1 );

         // ==== Accumulate orthogonal matrix in order update
         // .    H and Z, if requested.  ====

         if (NS.GT.1 .AND. S.NE.ZERO) CALL SORMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, WORK( JW+1 ), LWORK-JW, INFO );

         // ==== Update vertical slab in H ====

         if ( WANTT ) {
            LTOP = 1
         } else {
            LTOP = KTOP
         }
         DO 70 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            sgemm('N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), LDH, V, LDV, ZERO, WV, LDWV );
            slacpy('A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH );
         } // 70

         // ==== Update horizontal slab in H ====

         if ( WANTT ) {
            DO 80 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               sgemm('C', 'N', JW, KLN, JW, ONE, V, LDV, H( KWTOP, KCOL ), LDH, ZERO, T, LDT );
               slacpy('A', JW, KLN, T, LDT, H( KWTOP, KCOL ), LDH );
            } // 80
         }

         // ==== Update vertical slab in Z ====

         if ( WANTZ ) {
            DO 90 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               sgemm('N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), LDZ, V, LDV, ZERO, WV, LDWV );
               slacpy('A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), LDZ );
            } // 90
         }
      }

      // ==== Return the number of deflations ... ====

      ND = JW - NS

      // ==== ... and the number of shifts. (Subtracting
      // .    INFQR from the spike length takes care
      // .    of the case of a rare QR failure while
      // .    calculating eigenvalues of the deflation
      // .    window.)  ====

      NS = NS - INFQR

       // ==== Return optimal workspace. ====

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      // ==== End of SLAQR2 ====

      }
