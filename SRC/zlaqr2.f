      SUBROUTINE ZLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, LDZ, LWORK, N, ND, NH, NS, NV, NW;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), WORK( * ), WV( LDWV, * ), Z( LDZ, * )
      // ..

*  ================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0d0, 0.0d0 ), ONE = ( 1.0d0, 0.0d0 ) ;
      double             RZERO, RONE;
      const              RZERO = 0.0d0, RONE = 1.0d0 ;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         BETA, CDUM, S, TAU
      double             FOO, SAFMAX, SAFMIN, SMLNUM, ULP;
      int                I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN, KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWKOPT;
      // ..
      // .. External Functions ..
      double             DLAMCH;
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
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..

      // ==== Estimate optimal workspace. ====

      JW = MIN( NW, KBOT-KTOP+1 )
      if ( JW.LE.2 ) {
         LWKOPT = 1
      } else {

         // ==== Workspace query call to ZGEHRD ====

         zgehrd(JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO );
         LWK1 = INT( WORK( 1 ) )

         // ==== Workspace query call to ZUNMHR ====

         zunmhr('R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, WORK, -1, INFO );
         LWK2 = INT( WORK( 1 ) )

         // ==== Optimal workspace ====

         LWKOPT = JW + MAX( LWK1, LWK2 )
      }

      // ==== Quick return in case of workspace query. ====

      if ( LWORK.EQ.-1 ) {
         WORK( 1 ) = DCMPLX( LWKOPT, 0 )
         RETURN
      }

      // ==== Nothing to do ...
      // ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      IF( KTOP.GT.KBOT ) RETURN
      // ... nor for an empty deflation window. ====
      IF( NW.LT.1 ) RETURN

      // ==== Machine constants ====

      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )

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

         SH( KWTOP ) = H( KWTOP, KWTOP )
         NS = 1
         ND = 0
         if ( CABS1( S ).LE.MAX( SMLNUM, ULP*CABS1( H( KWTOP, KWTOP ) ) ) ) {
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP ) H( KWTOP, KWTOP-1 ) = ZERO
         }
         WORK( 1 ) = ONE
         RETURN
      }

      // ==== Convert to spike-triangular form.  (In case of a
      // .    rare QR failure, this routine continues to do
      // .    aggressive early deflation using that part of
      // .    the deflation window that converged using INFQR
      // .    here and there to keep track.) ====

      zlacpy('U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT );
      zcopy(JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 );

      zlaset('A', JW, JW, ZERO, ONE, V, LDV );
      zlahqr(.true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1, JW, V, LDV, INFQR );

      // ==== Deflation detection loop ====

      NS = JW
      ILST = INFQR + 1
      DO 10 KNT = INFQR + 1, JW

         // ==== Small spike tip deflation test ====

         FOO = CABS1( T( NS, NS ) )
         if ( FOO.EQ.RZERO ) FOO = CABS1( S )          IF( CABS1( S )*CABS1( V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) {

            // ==== One more converged eigenvalue ====

            NS = NS - 1
         } else {

            // ==== One undeflatable eigenvalue.  Move it up out of the
            // .    way.   (ZTREXC can not fail in this case.) ====

            IFST = NS
            ztrexc('V', JW, T, LDT, V, LDV, IFST, ILST, INFO );
            ILST = ILST + 1
         }
   10 CONTINUE

         // ==== Return to Hessenberg form ====

      IF( NS.EQ.0 ) S = ZERO

      if ( NS.LT.JW ) {

         // ==== sorting the diagonal of T improves accuracy for
         // .    graded matrices.  ====

         DO 30 I = INFQR + 1, NS
            IFST = I
            DO 20 J = I + 1, NS
               IF( CABS1( T( J, J ) ).GT.CABS1( T( IFST, IFST ) ) ) IFST = J
   20       CONTINUE
            ILST = I
            IF( IFST.NE.ILST ) CALL ZTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
   30    CONTINUE
      }

      // ==== Restore shift/eigenvalue array from T ====

      DO 40 I = INFQR + 1, JW
         SH( KWTOP+I-1 ) = T( I, I )
   40 CONTINUE


      if ( NS.LT.JW .OR. S.EQ.ZERO ) {
         if ( NS.GT.1 .AND. S.NE.ZERO ) {

            // ==== Reflect spike back into lower triangle ====

            zcopy(NS, V, LDV, WORK, 1 );
            DO 50 I = 1, NS
               WORK( I ) = DCONJG( WORK( I ) )
   50       CONTINUE
            BETA = WORK( 1 )
            zlarfg(NS, BETA, WORK( 2 ), 1, TAU );
            WORK( 1 ) = ONE

            zlaset('L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT );

            zlarf('L', NS, JW, WORK, 1, DCONJG( TAU ), T, LDT, WORK( JW+1 ) )             CALL ZLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT, WORK( JW+1 ) )             CALL ZLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV, WORK( JW+1 ) );

            zgehrd(JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), LWORK-JW, INFO );
         }

         // ==== Copy updated reduced window into place ====

         IF( KWTOP.GT.1 ) H( KWTOP, KWTOP-1 ) = S*DCONJG( V( 1, 1 ) )
         zlacpy('U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH );
         zcopy(JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), LDH+1 );

         // ==== Accumulate orthogonal matrix in order update
         // .    H and Z, if requested.  ====

         IF( NS.GT.1 .AND. S.NE.ZERO ) CALL ZUNMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, WORK( JW+1 ), LWORK-JW, INFO )

         // ==== Update vertical slab in H ====

         if ( WANTT ) {
            LTOP = 1
         } else {
            LTOP = KTOP
         }
         DO 60 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            zgemm('N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), LDH, V, LDV, ZERO, WV, LDWV );
            zlacpy('A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH );
   60    CONTINUE

         // ==== Update horizontal slab in H ====

         if ( WANTT ) {
            DO 70 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               zgemm('C', 'N', JW, KLN, JW, ONE, V, LDV, H( KWTOP, KCOL ), LDH, ZERO, T, LDT )                CALL ZLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ), LDH );
   70       CONTINUE
         }

         // ==== Update vertical slab in Z ====

         if ( WANTZ ) {
            DO 80 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               zgemm('N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), LDZ, V, LDV, ZERO, WV, LDWV )                CALL ZLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), LDZ );
   80       CONTINUE
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

      WORK( 1 ) = DCMPLX( LWKOPT, 0 )

      // ==== End of ZLAQR2 ====

      }
