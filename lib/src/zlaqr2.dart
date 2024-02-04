      void zlaqr2(WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, LDZ, LWORK, N, ND, NH, NS, NV, NW;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      Complex         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), WORK( * ), WV( LDWV, * ), Z( LDZ, * );
      // ..

// ================================================================

      // .. Parameters ..
      Complex         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      double             RZERO, RONE;
      const              RZERO = 0.0, RONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      Complex         BETA, CDUM, S, TAU;
      double             FOO, SAFMAX, SAFMIN, SMLNUM, ULP;
      int                I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN, KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWKOPT;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZGEHRD, ZGEMM, ZLACPY, ZLAHQR, ZLARF, ZLARFG, ZLASET, ZTREXC, ZUNMHR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, INT, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[CDUM] = ( CDUM.toDouble() ).abs() + ( DIMAG( CDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // ==== Estimate optimal workspace. ====

      JW = min( NW, KBOT-KTOP+1 );
      if ( JW <= 2 ) {
         LWKOPT = 1;
      } else {

         // ==== Workspace query call to ZGEHRD ====

         zgehrd(JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO );
         LWK1 = INT( WORK( 1 ) );

         // ==== Workspace query call to ZUNMHR ====

         zunmhr('R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, WORK, -1, INFO );
         LWK2 = INT( WORK( 1 ) );

         // ==== Optimal workspace ====

         LWKOPT = JW + max( LWK1, LWK2 );
      }

      // ==== Quick return in case of workspace query. ====

      if ( LWORK == -1 ) {
         WORK[1] = DCMPLX( LWKOPT, 0 );
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

      SAFMIN = DLAMCH( 'SAFE MINIMUM' );
      SAFMAX = RONE / SAFMIN;
      ULP = DLAMCH( 'PRECISION' );
      SMLNUM = SAFMIN*( N.toDouble() / ULP );

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

         SH[KWTOP] = H( KWTOP, KWTOP );
         NS = 1;
         ND = 0;
         if ( CABS1( S ) <= max( SMLNUM, ULP*CABS1( H( KWTOP, KWTOP ) ) ) ) {
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

      zlacpy('U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT );
      zcopy(JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 );

      zlaset('A', JW, JW, ZERO, ONE, V, LDV );
      zlahqr( true , true , JW, 1, JW, T, LDT, SH( KWTOP ), 1, JW, V, LDV, INFQR );

      // ==== Deflation detection loop ====

      NS = JW;
      ILST = INFQR + 1;
      for (KNT = INFQR + 1; KNT <= JW; KNT++) { // 10

         // ==== Small spike tip deflation test ====

         FOO = CABS1( T( NS, NS ) );
         if ( FOO == RZERO ) FOO = CABS1( S );
         IF( CABS1( S )*CABS1( V( 1, NS ) ) <= max( SMLNUM, ULP*FOO ) ) {

            // ==== One more converged eigenvalue ====

            NS = NS - 1;
         } else {

            // ==== One undeflatable eigenvalue.  Move it up out of the
            // .    way.   (ZTREXC can not fail in this case.) ====

            IFST = NS;
            ztrexc('V', JW, T, LDT, V, LDV, IFST, ILST, INFO );
            ILST = ILST + 1;
         }
      } // 10

         // ==== Return to Hessenberg form ====

      if (NS == 0) S = ZERO;

      if ( NS < JW ) {

         // ==== sorting the diagonal of T improves accuracy for
         // .    graded matrices.  ====

         for (I = INFQR + 1; I <= NS; I++) { // 30
            IFST = I;
            for (J = I + 1; J <= NS; J++) { // 20
               if( CABS1( T( J, J ) ) > CABS1( T( IFST, IFST ) ) ) IFST = J;
            } // 20
            ILST = I;
            if (IFST != ILST) ztrexc( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO );
         } // 30
      }

      // ==== Restore shift/eigenvalue array from T ====

      for (I = INFQR + 1; I <= JW; I++) { // 40
         SH[KWTOP+I-1] = T( I, I );
      } // 40


      if ( NS < JW || S == ZERO ) {
         if ( NS > 1 && S != ZERO ) {

            // ==== Reflect spike back into lower triangle ====

            zcopy(NS, V, LDV, WORK, 1 );
            for (I = 1; I <= NS; I++) { // 50
               WORK[I] = DCONJG( WORK( I ) );
            } // 50
            BETA = WORK( 1 );
            zlarfg(NS, BETA, WORK( 2 ), 1, TAU );
            WORK[1] = ONE;

            zlaset('L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT );

            zlarf('L', NS, JW, WORK, 1, DCONJG( TAU ), T, LDT, WORK( JW+1 ) );
            zlarf('R', NS, NS, WORK, 1, TAU, T, LDT, WORK( JW+1 ) );
            zlarf('R', JW, NS, WORK, 1, TAU, V, LDV, WORK( JW+1 ) );

            zgehrd(JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), LWORK-JW, INFO );
         }

         // ==== Copy updated reduced window into place ====

         if (KWTOP > 1) H( KWTOP, KWTOP-1 ) = S*DCONJG( V( 1, 1 ) );
         zlacpy('U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH );
         zcopy(JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), LDH+1 );

         // ==== Accumulate orthogonal matrix in order update
         // .    H and Z, if requested.  ====

         if (NS > 1 && S != ZERO) zunmhr( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, WORK( JW+1 ), LWORK-JW, INFO );

         // ==== Update vertical slab in H ====

         if ( WANTT ) {
            LTOP = 1;
         } else {
            LTOP = KTOP;
         }
         for (KROW = LTOP; NV < 0 ? KROW >= KWTOP - 1 : KROW <= KWTOP - 1; KROW += NV) { // 60
            KLN = min( NV, KWTOP-KROW );
            zgemm('N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), LDH, V, LDV, ZERO, WV, LDWV );
            zlacpy('A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH );
         } // 60

         // ==== Update horizontal slab in H ====

         if ( WANTT ) {
            for (KCOL = KBOT + 1; NH < 0 ? KCOL >= N : KCOL <= N; KCOL += NH) { // 70
               KLN = min( NH, N-KCOL+1 );
               zgemm('C', 'N', JW, KLN, JW, ONE, V, LDV, H( KWTOP, KCOL ), LDH, ZERO, T, LDT );
               zlacpy('A', JW, KLN, T, LDT, H( KWTOP, KCOL ), LDH );
            } // 70
         }

         // ==== Update vertical slab in Z ====

         if ( WANTZ ) {
            for (KROW = ILOZ; NV < 0 ? KROW >= IHIZ : KROW <= IHIZ; KROW += NV) { // 80
               KLN = min( NV, IHIZ-KROW+1 );
               zgemm('N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), LDZ, V, LDV, ZERO, WV, LDWV );
               zlacpy('A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), LDZ );
            } // 80
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

      WORK[1] = DCMPLX( LWKOPT, 0 );

      // ==== End of ZLAQR2 ====

      }