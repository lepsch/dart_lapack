      SUBROUTINE DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      bool               WANTT, WANTZ;
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ), Z( LDZ, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    DLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      int                NTINY
      PARAMETER          ( NTINY = 15 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      int                KEXNW
      PARAMETER          ( KEXNW = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      int                KEXSH
      PARAMETER          ( KEXSH = 6 )
*
*     ==== The constants WILK1 and WILK2 are used to form the
*     .    exceptional shifts. ====
      DOUBLE PRECISION   WILK1, WILK2
      PARAMETER          ( WILK1 = 0.75d0, WILK2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      int                I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS, LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS, NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD
      bool               SORTED;
      String             JBCMPZ*2;
*     ..
*     .. External Functions ..
      int                ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLAHQR, DLANV2, DLAQR3, DLAQR4, DLAQR5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      INFO = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF
*
      IF( N.LE.NTINY ) THEN
*
*        ==== Tiny matrices must use DLAHQR. ====
*
         LWKOPT = 1
         IF( LWORK.NE.-1 ) CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILOZ, IHIZ, Z, LDZ, INFO )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         INFO = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( WANTT ) THEN
            JBCMPZ( 1: 1 ) = 'S'
         ELSE
            JBCMPZ( 1: 1 ) = 'E'
         END IF
         IF( WANTZ ) THEN
            JBCMPZ( 2: 2 ) = 'V'
         ELSE
            JBCMPZ( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 15, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.4.) ====
*
         NWR = ILAENV( 13, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 15, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         NSR = ILAENV( 15, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N-3 ) / 6, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to DLAQR3 ====
*
         CALL DLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ, IHIZ, Z, LDZ, LS, LD, WR, WI, H, LDH, N, H, LDH, N, H, LDH, WORK, -1 )
*
*        ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
*
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = DBLE( LWKOPT )
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         NMIN = ILAENV( 12, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
*
*        ==== Nibble crossover point ====
*
         NIBBLE = ILAENV( 14, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         KACC22 = ILAENV( 16, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
         NW = NWMAX
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         NSMAX = MIN( ( N-3 ) / 6, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         NDFL = 1
*
*        ==== ITMAX = iteration limit ====
*
         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
*
*        ==== Last row and column in the active block ====
*
         KBOT = IHI
*
*        ==== Main Loop ====
*
         DO 80 IT = 1, ITMAX
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( KBOT.LT.ILO ) GO TO 90
*
*           ==== Locate active block ====
*
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO ) GO TO 20
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            NH = KBOT - KTOP + 1
            NWUPBD = MIN( NH, NWMAX )
            IF( NDFL.LT.KEXNW ) THEN
               NW = MIN( NWUPBD, NWR )
            ELSE
               NW = MIN( NWUPBD, 2*NW )
            END IF
            IF( NW.LT.NWMAX ) THEN
               IF( NW.GE.NH-1 ) THEN
                  NW = NH
               ELSE
                  KWTOP = KBOT - NW + 1
                  IF( ABS( H( KWTOP, KWTOP-1 ) ).GT. ABS( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               END IF
            END IF
            IF( NDFL.LT.KEXNW ) THEN
               NDEC = -1
            ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
               NDEC = NDEC + 1
               IF( NW-NDEC.LT.2 ) NDEC = 0
               NW = NW - NDEC
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1
*
*           ==== Aggressive early deflation ====
*
            CALL DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, LS, LD, WR, WI, H( KV, 1 ), LDH, NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, WORK, LWORK )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            KBOT = KBOT - LD
*
*           ==== KS points to the shifts. ====
*
            KS = KBOT - LS + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT- KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if DLAQR3
*              .    did not provide that many shifts. ====
*
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    DLAQR3 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, MAX( KS+1, KTOP+2 ), -2
                     SS = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
                     AA = WILK1*SS + H( I, I )
                     BB = SS
                     CC = WILK2*SS
                     DD = AA
                     CALL DLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ), WR( I ), WI( I ), CS, SN )
   30             CONTINUE
                  IF( KS.EQ.KTOP ) THEN
                     WR( KS+1 ) = H( KS+1, KS+1 )
                     WI( KS+1 ) = ZERO
                     WR( KS ) = WR( KS+1 )
                     WI( KS ) = WI( KS+1 )
                  END IF
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use DLAQR4 or
*                 .    DLAHQR on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL DLACPY( 'A', NS, NS, H( KS, KS ), LDH, H( KT, 1 ), LDH )
                     IF( NS.GT.NMIN ) THEN
                        CALL DLAQR4( .false., .false., NS, 1, NS, H( KT, 1 ), LDH, WR( KS ), WI( KS ), 1, 1, ZDUM, 1, WORK, LWORK, INF )
                     ELSE
                        CALL DLAHQR( .false., .false., NS, 1, NS, H( KT, 1 ), LDH, WR( KS ), WI( KS ), 1, 1, ZDUM, 1, INF )
                     END IF
                     KS = KS + INF
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  ====
*
                     IF( KS.GE.KBOT ) THEN
                        AA = H( KBOT-1, KBOT-1 )
                        CC = H( KBOT, KBOT-1 )
                        BB = H( KBOT-1, KBOT )
                        DD = H( KBOT, KBOT )
                        CALL DLANV2( AA, BB, CC, DD, WR( KBOT-1 ), WI( KBOT-1 ), WR( KBOT ), WI( KBOT ), CS, SN )
                        KS = KBOT - 1
                     END IF
                  END IF
*
                  IF( KBOT-KS+1.GT.NS ) THEN
*
*                    ==== Sort the shifts (Helps a little)
*                    .    Bubble sort keeps complex conjugate
*                    .    pairs together. ====
*
                     SORTED = .false.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED ) GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( ABS( WR( I ) )+ABS( WI( I ) ).LT. ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
                              SORTED = .false.
*
                              SWAP = WR( I )
                              WR( I ) = WR( I+1 )
                              WR( I+1 ) = SWAP
*
                              SWAP = WI( I )
                              WI( I ) = WI( I+1 )
                              WI( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
*
*                 ==== Shuffle shifts into pairs of real shifts
*                 .    and pairs of complex conjugate shifts
*                 .    assuming complex conjugate shifts are
*                 .    already adjacent to one another. (Yes,
*                 .    they are.)  ====
*
                  DO 70 I = KBOT, KS + 2, -2
                     IF( WI( I ).NE.-WI( I-1 ) ) THEN
*
                        SWAP = WR( I )
                        WR( I ) = WR( I-1 )
                        WR( I-1 ) = WR( I-2 )
                        WR( I-2 ) = SWAP
*
                        SWAP = WI( I )
                        WI( I ) = WI( I-1 )
                        WI( I-1 ) = WI( I-2 )
                        WI( I-2 ) = SWAP
                     END IF
   70             CONTINUE
               END IF
*
*              ==== If there are only two shifts and both are
*              .    real, then use only one.  ====
*
               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( WI( KBOT ).EQ.ZERO ) THEN
                     IF( ABS( WR( KBOT )-H( KBOT, KBOT ) ).LT. ABS( WR( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                        WR( KBOT-1 ) = WR( KBOT )
                     ELSE
                        WR( KBOT ) = WR( KBOT-1 )
                     END IF
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnitude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               NS = MIN( NS, KBOT-KS+1 )
               NS = NS - MOD( NS, 2 )
               KS = KBOT - NS + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               KDU = 2*NS
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS, WR( KS ), WI( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK, 3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF
*
*           ==== End of main loop ====
   80    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         INFO = KBOT
   90    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      WORK( 1 ) = DBLE( LWKOPT )
*
*     ==== End of DLAQR0 ====
*
      END
