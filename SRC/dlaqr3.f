      SUBROUTINE DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, LDZ, LWORK, N, ND, NH, NS, NV, NW;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      double             H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), V( LDV, * ), WORK( * ), WV( LDWV, * ), Z( LDZ, * );
      // ..

*  ================================================================
      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0d0, ONE = 1.0d0 ;
      // ..
      // .. Local Scalars ..
      double             AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP;
      int                I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3, LWKOPT, NMIN;
      bool               BULGE, SORTED;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      int                ILAENV;
      // EXTERNAL DLAMCH, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEHRD, DGEMM, DLACPY, DLAHQR, DLANV2, DLAQR4, DLARF, DLARFG, DLASET, DORMHR, DTREXC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // ==== Estimate optimal workspace. ====

      JW = MIN( NW, KBOT-KTOP+1 )
      if ( JW.LE.2 ) {
         LWKOPT = 1
      } else {

         // ==== Workspace query call to DGEHRD ====

         dgehrd(JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO );
         LWK1 = INT( WORK( 1 ) )

         // ==== Workspace query call to DORMHR ====

         dormhr('R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, WORK, -1, INFO );
         LWK2 = INT( WORK( 1 ) )

         // ==== Workspace query call to DLAQR4 ====

         dlaqr4( true , true , JW, 1, JW, T, LDT, SR, SI, 1, JW, V, LDV, WORK, -1, INFQR );
         LWK3 = INT( WORK( 1 ) )

         // ==== Optimal workspace ====

         LWKOPT = MAX( JW+MAX( LWK1, LWK2 ), LWK3 )
      }

      // ==== Quick return in case of workspace query. ====

      if ( LWORK == -1 ) {
         WORK( 1 ) = DBLE( LWKOPT )
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

      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )

      // ==== Setup deflation window ====

      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      if ( KWTOP == KTOP ) {
         S = ZERO
      } else {
         S = H( KWTOP, KWTOP-1 )
      }

      if ( KBOT == KWTOP ) {

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

      dlacpy('U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT );
      dcopy(JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 );

      dlaset('A', JW, JW, ZERO, ONE, V, LDV );
      NMIN = ILAENV( 12, 'DLAQR3', 'SV', JW, 1, JW, LWORK )
      if ( JW.GT.NMIN ) {
         dlaqr4( true , true , JW, 1, JW, T, LDT, SR( KWTOP ), SI( KWTOP ), 1, JW, V, LDV, WORK, LWORK, INFQR );
      } else {
         dlahqr( true , true , JW, 1, JW, T, LDT, SR( KWTOP ), SI( KWTOP ), 1, JW, V, LDV, INFQR );
      }

      // ==== DTREXC needs a clean margin near the diagonal ====

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
         if ( NS == 1 ) {
            BULGE = false;
         } else {
            BULGE = T( NS, NS-1 ) != ZERO
         }

         // ==== Small spike tip test for deflation ====

         if ( .NOT. BULGE ) {

            // ==== Real eigenvalue ====

            FOO = ABS( T( NS, NS ) )
            if (FOO == ZERO) FOO = ABS( S );
            if ( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) {

               // ==== Deflatable ====

               NS = NS - 1
            } else {

               // ==== Undeflatable.   Move it up out of the way.
               // .    (DTREXC can not fail in this case.) ====

               IFST = NS
               dtrexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               ILST = ILST + 1
            }
         } else {

            // ==== Complex conjugate pair ====

            FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )* SQRT( ABS( T( NS-1, NS ) ) )             IF( FOO == ZERO ) FOO = ABS( S )             IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE. MAX( SMLNUM, ULP*FOO ) ) THEN

               // ==== Deflatable ====

               NS = NS - 2
            } else {

               // ==== Undeflatable. Move them up out of the way.
               // .    Fortunately, DTREXC does the right thing with
               // .    ILST in case of a rare exchange failure. ====

               IFST = NS
               dtrexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               ILST = ILST + 2
            }
         }

         // ==== End deflation detection loop ====

         GO TO 20
      }

         // ==== Return to Hessenberg form ====

      if (NS == 0) S = ZERO;

      if ( NS.LT.JW ) {

         // ==== sorting diagonal blocks of T improves accuracy for
         // .    graded matrices.  Bubble sort deals well with
         // .    exchange failures. ====

         SORTED = false;
         I = NS + 1
         } // 30
         if (SORTED) GO TO 50;
         SORTED = true;

         KEND = I - 1
         I = INFQR + 1
         if ( I == NS ) {
            K = I + 1
         } else if ( T( I+1, I ) == ZERO ) {
            K = I + 1
         } else {
            K = I + 2
         }
         } // 40
         if ( K.LE.KEND ) {
            if ( K == I+1 ) {
               EVI = ABS( T( I, I ) )
            } else {
               EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )* SQRT( ABS( T( I, I+1 ) ) )
            }

            if ( K == KEND ) {
               EVK = ABS( T( K, K ) )
            } else if ( T( K+1, K ) == ZERO ) {
               EVK = ABS( T( K, K ) )
            } else {
               EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )* SQRT( ABS( T( K, K+1 ) ) )
            }

            if ( EVI.GE.EVK ) {
               I = K
            } else {
               SORTED = false;
               IFST = I
               ILST = K
               dtrexc('V', JW, T, LDT, V, LDV, IFST, ILST, WORK, INFO );
               if ( INFO == 0 ) {
                  I = ILST
               } else {
                  I = K
               }
            }
            if ( I == KEND ) {
               K = I + 1
            } else if ( T( I+1, I ) == ZERO ) {
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
         if ( I == INFQR+1 ) {
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         } else if ( T( I, I-1 ) == ZERO ) {
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         } else {
            AA = T( I-1, I-1 )
            CC = T( I, I-1 )
            BB = T( I-1, I )
            DD = T( I, I )
            dlanv2(AA, BB, CC, DD, SR( KWTOP+I-2 ), SI( KWTOP+I-2 ), SR( KWTOP+I-1 ), SI( KWTOP+I-1 ), CS, SN );
            I = I - 2
         }
         GO TO 60
      }

      if ( NS.LT.JW .OR. S == ZERO ) {
         if ( NS.GT.1 .AND. S != ZERO ) {

            // ==== Reflect spike back into lower triangle ====

            dcopy(NS, V, LDV, WORK, 1 );
            BETA = WORK( 1 )
            dlarfg(NS, BETA, WORK( 2 ), 1, TAU );
            WORK( 1 ) = ONE

            dlaset('L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT );

            dlarf('L', NS, JW, WORK, 1, TAU, T, LDT, WORK( JW+1 ) );
            dlarf('R', NS, NS, WORK, 1, TAU, T, LDT, WORK( JW+1 ) );
            dlarf('R', JW, NS, WORK, 1, TAU, V, LDV, WORK( JW+1 ) );

            dgehrd(JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), LWORK-JW, INFO );
         }

         // ==== Copy updated reduced window into place ====

         if (KWTOP.GT.1) H( KWTOP, KWTOP-1 ) = S*V( 1, 1 );
         dlacpy('U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH );
         dcopy(JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), LDH+1 );

         // ==== Accumulate orthogonal matrix in order update
         // .    H and Z, if requested.  ====

         if (NS.GT.1 .AND. S != ZERO) CALL DORMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, WORK( JW+1 ), LWORK-JW, INFO );

         // ==== Update vertical slab in H ====

         if ( WANTT ) {
            LTOP = 1
         } else {
            LTOP = KTOP
         }
         DO 70 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            dgemm('N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), LDH, V, LDV, ZERO, WV, LDWV );
            dlacpy('A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH );
         } // 70

         // ==== Update horizontal slab in H ====

         if ( WANTT ) {
            DO 80 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               dgemm('C', 'N', JW, KLN, JW, ONE, V, LDV, H( KWTOP, KCOL ), LDH, ZERO, T, LDT );
               dlacpy('A', JW, KLN, T, LDT, H( KWTOP, KCOL ), LDH );
            } // 80
         }

         // ==== Update vertical slab in Z ====

         if ( WANTZ ) {
            DO 90 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               dgemm('N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), LDZ, V, LDV, ZERO, WV, LDWV );
               dlacpy('A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), LDZ );
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

      WORK( 1 ) = DBLE( LWKOPT )

      // ==== End of DLAQR3 ====

      }
