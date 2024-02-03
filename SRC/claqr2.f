      SUBROUTINE CLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT, NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, LDZ, LWORK, N, ND, NH, NS, NV, NW;
      bool               WANTT, WANTZ;
*     ..
*     .. Array Arguments ..
      COMPLEX            H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ), WORK( * ), WV( LDWV, * ), Z( LDZ, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0e0, 0.0e0 ), ONE = ( 1.0e0, 0.0e0 ) )
      REAL               RZERO, RONE
      PARAMETER          ( RZERO = 0.0e0, RONE = 1.0e0 )
*     ..
*     .. Local Scalars ..
      COMPLEX            BETA, CDUM, S, TAU
      REAL               FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      int                I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN, KNT, KROW, KWTOP, LTOP, LWK1, LWK2, LWKOPT;
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEHRD, CGEMM, CLACPY, CLAHQR, CLARF, CLARFG, CLASET, CTREXC, CUNMHR
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, INT, MAX, MIN, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      JW = MIN( NW, KBOT-KTOP+1 )
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
*
*        ==== Workspace query call to CGEHRD ====
*
         CALL CGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
*
*        ==== Workspace query call to CUNMHR ====
*
         CALL CUNMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
*
*        ==== Optimal workspace ====
*
         LWKOPT = JW + MAX( LWK1, LWK2 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( LWORK.EQ.-1 ) THEN
         WORK( 1 ) = CMPLX( LWKOPT, 0 )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      IF( KTOP.GT.KBOT ) RETURN
*     ... nor for an empty deflation window. ====
      IF( NW.LT.1 ) RETURN
*
*     ==== Machine constants ====
*
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = RONE / SAFMIN
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( N ) / ULP )
*
*     ==== Setup deflation window ====
*
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         S = H( KWTOP, KWTOP-1 )
      END IF
*
      IF( KBOT.EQ.KWTOP ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         SH( KWTOP ) = H( KWTOP, KWTOP )
         NS = 1
         ND = 0
         IF( CABS1( S ).LE.MAX( SMLNUM, ULP*CABS1( H( KWTOP, KWTOP ) ) ) ) THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP ) H( KWTOP, KWTOP-1 ) = ZERO
         END IF
         WORK( 1 ) = ONE
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL CLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL CCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
*
      CALL CLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      CALL CLAHQR( .true., .true., JW, 1, JW, T, LDT, SH( KWTOP ), 1, JW, V, LDV, INFQR )
*
*     ==== Deflation detection loop ====
*
      NS = JW
      ILST = INFQR + 1
      DO 10 KNT = INFQR + 1, JW
*
*        ==== Small spike tip deflation test ====
*
         FOO = CABS1( T( NS, NS ) )
         IF( FOO.EQ.RZERO ) FOO = CABS1( S )          IF( CABS1( S )*CABS1( V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
*
*           ==== One more converged eigenvalue ====
*
            NS = NS - 1
         ELSE
*
*           ==== One undeflatable eigenvalue.  Move it up out of the
*           .    way.   (CTREXC can not fail in this case.) ====
*
            IFST = NS
            CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
            ILST = ILST + 1
         END IF
   10 CONTINUE
*
*        ==== Return to Hessenberg form ====
*
      IF( NS.EQ.0 ) S = ZERO
*
      IF( NS.LT.JW ) THEN
*
*        ==== sorting the diagonal of T improves accuracy for
*        .    graded matrices.  ====
*
         DO 30 I = INFQR + 1, NS
            IFST = I
            DO 20 J = I + 1, NS
               IF( CABS1( T( J, J ) ).GT.CABS1( T( IFST, IFST ) ) ) IFST = J
   20       CONTINUE
            ILST = I
            IF( IFST.NE.ILST ) CALL CTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, INFO )
   30    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      DO 40 I = INFQR + 1, JW
         SH( KWTOP+I-1 ) = T( I, I )
   40 CONTINUE
*
*
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL CCOPY( NS, V, LDV, WORK, 1 )
            DO 50 I = 1, NS
               WORK( I ) = CONJG( WORK( I ) )
   50       CONTINUE
            BETA = WORK( 1 )
            CALL CLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
*
            CALL CLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
*
            CALL CLARF( 'L', NS, JW, WORK, 1, CONJG( TAU ), T, LDT, WORK( JW+1 ) )             CALL CLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT, WORK( JW+1 ) )             CALL CLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV, WORK( JW+1 ) )
*
            CALL CGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), LWORK-JW, INFO )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( KWTOP.GT.1 ) H( KWTOP, KWTOP-1 ) = S*CONJG( V( 1, 1 ) )
         CALL CLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL CCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), LDH+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( NS.GT.1 .AND. S.NE.ZERO ) CALL CUNMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, WORK( JW+1 ), LWORK-JW, INFO )
*
*        ==== Update vertical slab in H ====
*
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         DO 60 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), LDH, V, LDV, ZERO, WV, LDWV )
            CALL CLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   60    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( WANTT ) THEN
            DO 70 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL CGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV, H( KWTOP, KCOL ), LDH, ZERO, T, LDT )                CALL CLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ), LDH )
   70       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( WANTZ ) THEN
            DO 80 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL CGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), LDZ, V, LDV, ZERO, WV, LDWV )                CALL CLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), LDZ )
   80       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      ND = JW - NS
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      NS = NS - INFQR
*
*      ==== Return optimal workspace. ====
*
      WORK( 1 ) = CMPLX( LWKOPT, 0 )
*
*     ==== End of CLAQR2 ====
*
      END
