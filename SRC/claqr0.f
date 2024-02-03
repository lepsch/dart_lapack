      SUBROUTINE CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  ================================================================
      // .. Parameters ..

      // ==== Matrices of order NTINY or smaller must be processed by
      // .    CLAHQR because of insufficient subdiagonal scratch space.
      // .    (This is a hard limit.) ====
      int                NTINY;
      const              NTINY = 15 ;

      // ==== Exceptional deflation windows:  try to cure rare
      // .    slow convergence by varying the size of the
      // .    deflation window after KEXNW iterations. ====
      int                KEXNW;
      const              KEXNW = 5 ;

      // ==== Exceptional shifts: try to cure rare slow convergence
      // .    with ad-hoc exceptional shifts every KEXSH iterations.
      // .    ====
      int                KEXSH;
      const              KEXSH = 6 ;

      // ==== The constant WILK1 is used to form the exceptional
      // .    shifts. ====
      REAL               WILK1
      const              WILK1 = 0.75e0 ;
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0e0, 0.0e0 ), ONE = ( 1.0e0, 0.0e0 ) ;
      REAL               TWO
      const              TWO = 2.0e0 ;
      // ..
      // .. Local Scalars ..
      COMPLEX            AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      REAL               S
      int                I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS, LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS, NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD;
      bool               SORTED;
      String             JBCMPZ*2;
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Local Arrays ..
      COMPLEX            ZDUM( 1, 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACPY, CLAHQR, CLAQR3, CLAQR4, CLAQR5
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, INT, MAX, MIN, MOD, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( REAL( CDUM ) ) + ABS( AIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..
      INFO = 0

      // ==== Quick return for N = 0: nothing to do. ====

      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF

      IF( N.LE.NTINY ) THEN

         // ==== Tiny matrices must use CLAHQR. ====

         LWKOPT = 1
         IF( LWORK.NE.-1 ) CALL CLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ, INFO )
      ELSE

         // ==== Use small bulge multi-shift QR with aggressive early
         // .    deflation on larger-than-tiny matrices. ====

         // ==== Hope for the best. ====

         INFO = 0

         // ==== Set up job flags for ILAENV. ====

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

         // ==== NWR = recommended deflation window size.  At this
         // .    point,  N .GT. NTINY = 15, so there is enough
         // .    subdiagonal workspace for NWR.GE.2 as required.
         // .    (In fact, there is enough subdiagonal space for
         // .    NWR.GE.4.) ====

         NWR = ILAENV( 13, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )

         // ==== NSR = recommended number of simultaneous shifts.
         // .    At this point N .GT. NTINY = 15, so there is at
         // .    enough subdiagonal workspace for NSR to be even
         // .    and greater than or equal to two as required. ====

         NSR = ILAENV( 15, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N-3 ) / 6, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )

         // ==== Estimate optimal workspace ====

         // ==== Workspace query call to CLAQR3 ====

         CALL CLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ, IHIZ, Z, LDZ, LS, LD, W, H, LDH, N, H, LDH, N, H, LDH, WORK, -1 )

         // ==== Optimal workspace = MAX(CLAQR5, CLAQR3) ====

         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )

         // ==== Quick return in case of workspace query. ====

         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = CMPLX( LWKOPT, 0 )
            RETURN
         END IF

         // ==== CLAHQR/CLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )

         // ==== Nibble crossover point ====

         NIBBLE = ILAENV( 14, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )

         // ==== Accumulate reflections during ttswp?  Use block
         // .    2-by-2 structure during matrix-matrix multiply? ====

         KACC22 = ILAENV( 16, 'CLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )

         // ==== NWMAX = the largest possible deflation window for
         // .    which there is sufficient workspace. ====

         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
         NW = NWMAX

         // ==== NSMAX = the Largest number of simultaneous shifts
         // .    for which there is sufficient workspace. ====

         NSMAX = MIN( ( N-3 ) / 6, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )

         // ==== NDFL: an iteration count restarted at deflation. ====

         NDFL = 1

         // ==== ITMAX = iteration limit ====

         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )

         // ==== Last row and column in the active block ====

         KBOT = IHI

         // ==== Main Loop ====

         DO 70 IT = 1, ITMAX

            // ==== Done when KBOT falls below ILO ====

            IF( KBOT.LT.ILO ) GO TO 80

            // ==== Locate active block ====

            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO ) GO TO 20
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K

            // ==== Select deflation window size:
            // .    Typical Case:
            // .      If possible and advisable, nibble the entire
            // .      active block.  If not, use size MIN(NWR,NWMAX)
            // .      or MIN(NWR+1,NWMAX) depending upon which has
            // .      the smaller corresponding subdiagonal entry
            // .      (a heuristic).
            // .
            // .    Exceptional Case:
            // .      If there have been no deflations in KEXNW or
            // .      more iterations, then vary the deflation window
            // .      size.   At first, because, larger windows are,
            // .      in general, more powerful than smaller ones,
            // .      rapidly increase the window to the maximum possible.
            // .      Then, gradually reduce the window size. ====

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
                  IF( CABS1( H( KWTOP, KWTOP-1 ) ).GT. CABS1( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               END IF
            END IF
            IF( NDFL.LT.KEXNW ) THEN
               NDEC = -1
            ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
               NDEC = NDEC + 1
               IF( NW-NDEC.LT.2 ) NDEC = 0
               NW = NW - NDEC
            END IF

            // ==== Aggressive early deflation:
            // .    split workspace under the subdiagonal into
            // .      - an nw-by-nw work array V in the lower
            // .        left-hand-corner,
            // .      - an NW-by-at-least-NW-but-more-is-better
            // .        (NW-by-NHO) horizontal work array along
            // .        the bottom edge,
            // .      - an at-least-NW-but-more-is-better (NHV-by-NW)
            // .        vertical work array along the left-hand-edge.
            // .        ====

            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1

            // ==== Aggressive early deflation ====

            CALL CLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, LS, LD, W, H( KV, 1 ), LDH, NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, WORK, LWORK )

            // ==== Adjust KBOT accounting for new deflations. ====

            KBOT = KBOT - LD

            // ==== KS points to the shifts. ====

            KS = KBOT - LS + 1

            // ==== Skip an expensive QR sweep if there is a (partly
            // .    heuristic) reason to expect that many eigenvalues
            // .    will deflate without it.  Here, the QR sweep is
            // .    skipped if many eigenvalues have just been deflated
            // .    or if the remaining active block is small.

            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT- KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN

               // ==== NS = nominal number of simultaneous shifts.
               // .    This may be lowered (slightly) if CLAQR3
               // .    did not provide that many shifts. ====

               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )

               // ==== If there have been no deflations
               // .    in a multiple of KEXSH iterations,
               // .    then try exceptional shifts.
               // .    Otherwise use shifts provided by
               // .    CLAQR3 above or from the eigenvalues
               // .    of a trailing principal submatrix. ====

               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, KS + 1, -2
                     W( I ) = H( I, I ) + WILK1*CABS1( H( I, I-1 ) )
                     W( I-1 ) = W( I )
   30             CONTINUE
               ELSE

                  // ==== Got NS/2 or fewer shifts? Use CLAQR4 or
                  // .    CLAHQR on a trailing principal submatrix to
                  // .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
                  // .    there is enough space below the subdiagonal
                  // .    to fit an NS-by-NS scratch array.) ====

                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL CLACPY( 'A', NS, NS, H( KS, KS ), LDH, H( KT, 1 ), LDH )
                     IF( NS.GT.NMIN ) THEN
                        CALL CLAQR4( .false., .false., NS, 1, NS, H( KT, 1 ), LDH, W( KS ), 1, 1, ZDUM, 1, WORK, LWORK, INF )
                     ELSE
                        CALL CLAHQR( .false., .false., NS, 1, NS, H( KT, 1 ), LDH, W( KS ), 1, 1, ZDUM, 1, INF )
                     END IF
                     KS = KS + INF

                     // ==== In case of a rare QR failure use
                     // .    eigenvalues of the trailing 2-by-2
                     // .    principal submatrix.  Scale to avoid
                     // .    overflows, underflows and subnormals.
                     // .    (The scale factor S can not be zero,
                     // .    because H(KBOT,KBOT-1) is nonzero.) ====

                     IF( KS.GE.KBOT ) THEN
                        S = CABS1( H( KBOT-1, KBOT-1 ) ) + CABS1( H( KBOT, KBOT-1 ) ) + CABS1( H( KBOT-1, KBOT ) ) + CABS1( H( KBOT, KBOT ) )
                        AA = H( KBOT-1, KBOT-1 ) / S
                        CC = H( KBOT, KBOT-1 ) / S
                        BB = H( KBOT-1, KBOT ) / S
                        DD = H( KBOT, KBOT ) / S
                        TR2 = ( AA+DD ) / TWO
                        DET = ( AA-TR2 )*( DD-TR2 ) - BB*CC
                        RTDISC = SQRT( -DET )
                        W( KBOT-1 ) = ( TR2+RTDISC )*S
                        W( KBOT ) = ( TR2-RTDISC )*S

                        KS = KBOT - 1
                     END IF
                  END IF

                  IF( KBOT-KS+1.GT.NS ) THEN

                     // ==== Sort the shifts (Helps a little) ====

                     SORTED = .false.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED ) GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( CABS1( W( I ) ).LT.CABS1( W( I+1 ) ) ) THEN
                              SORTED = .false.
                              SWAP = W( I )
                              W( I ) = W( I+1 )
                              W( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
               END IF

               // ==== If there are only two shifts, then use
               // .    only one.  ====

               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( CABS1( W( KBOT )-H( KBOT, KBOT ) ).LT. CABS1( W( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                     W( KBOT-1 ) = W( KBOT )
                  ELSE
                     W( KBOT ) = W( KBOT-1 )
                  END IF
               END IF

               // ==== Use up to NS of the the smallest magnitude
               // .    shifts.  If there aren't NS shifts available,
               // .    then use them all, possibly dropping one to
               // .    make the number of shifts even. ====

               NS = MIN( NS, KBOT-KS+1 )
               NS = NS - MOD( NS, 2 )
               KS = KBOT - NS + 1

               // ==== Small-bulge multi-shift QR sweep:
               // .    split workspace under the subdiagonal into
               // .    - a KDU-by-KDU work array U in the lower
               // .      left-hand-corner,
               // .    - a KDU-by-at-least-KDU-but-more-is-better
               // .      (KDU-by-NHo) horizontal work array WH along
               // .      the bottom edge,
               // .    - and an at-least-KDU-but-more-is-better-by-KDU
               // .      (NVE-by-KDU) vertical work WV arrow along
               // .      the left-hand-edge. ====

               KDU = 2*NS
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1

               // ==== Small-bulge multi-shift QR sweep ====

               CALL CLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS, W( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK, 3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
            END IF

            // ==== Note progress (or the lack of it). ====

            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF

            // ==== End of main loop ====
   70    CONTINUE

         // ==== Iteration limit exceeded.  Set INFO to show where
         // .    the problem occurred and exit. ====

         INFO = KBOT
   80    CONTINUE
      END IF

      // ==== Return the optimal value of LWORK. ====

      WORK( 1 ) = CMPLX( LWKOPT, 0 )

      // ==== End of CLAQR0 ====

      }
