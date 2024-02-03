      SUBROUTINE ZLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N;
      bool               WANTT, WANTZ;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  ================================================================

      // .. Parameters ..

      // ==== Matrices of order NTINY or smaller must be processed by
      // .    ZLAHQR because of insufficient subdiagonal scratch space.
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
      double             WILK1;
      const              WILK1 = 0.75d0 ;
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0d0, 0.0d0 ), ONE = ( 1.0d0, 0.0d0 ) ;
      double             TWO;
      const              TWO = 2.0d0 ;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      double             S;
      int                I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS, LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS, NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD;
      bool               SORTED;
      String             JBCMPZ*2;
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Local Arrays ..
      COMPLEX*16         ZDUM( 1, 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACPY, ZLAHQR, ZLAQR2, ZLAQR5
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG, INT, MAX, MIN, MOD, SQRT
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
      // ..
      // .. Executable Statements ..
      INFO = 0

      // ==== Quick return for N = 0: nothing to do. ====

      if ( N == 0 ) {
         WORK( 1 ) = ONE
         RETURN
      }

      if ( N.LE.NTINY ) {

         // ==== Tiny matrices must use ZLAHQR. ====

         LWKOPT = 1
         if (LWORK != -1) CALL ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ, IHIZ, Z, LDZ, INFO );
      } else {

         // ==== Use small bulge multi-shift QR with aggressive early
         // .    deflation on larger-than-tiny matrices. ====

         // ==== Hope for the best. ====

         INFO = 0

         // ==== Set up job flags for ILAENV. ====

         if ( WANTT ) {
            JBCMPZ( 1: 1 ) = 'S'
         } else {
            JBCMPZ( 1: 1 ) = 'E'
         }
         if ( WANTZ ) {
            JBCMPZ( 2: 2 ) = 'V'
         } else {
            JBCMPZ( 2: 2 ) = 'N'
         }

         // ==== NWR = recommended deflation window size.  At this
         // .    point,  N .GT. NTINY = 15, so there is enough
         // .    subdiagonal workspace for NWR.GE.2 as required.
         // .    (In fact, there is enough subdiagonal space for
         // .    NWR.GE.4.) ====

         NWR = ILAENV( 13, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )

         // ==== NSR = recommended number of simultaneous shifts.
         // .    At this point N .GT. NTINY = 15, so there is at
         // .    enough subdiagonal workspace for NSR to be even
         // .    and greater than or equal to two as required. ====

         NSR = ILAENV( 15, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N-3 ) / 6, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )

         // ==== Estimate optimal workspace ====

         // ==== Workspace query call to ZLAQR2 ====

         zlaqr2(WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ, IHIZ, Z, LDZ, LS, LD, W, H, LDH, N, H, LDH, N, H, LDH, WORK, -1 );

         // ==== Optimal workspace = MAX(ZLAQR5, ZLAQR2) ====

         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )

         // ==== Quick return in case of workspace query. ====

         if ( LWORK == -1 ) {
            WORK( 1 ) = DCMPLX( LWKOPT, 0 )
            RETURN
         }

         // ==== ZLAHQR/ZLAQR0 crossover point ====

         NMIN = ILAENV( 12, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )

         // ==== Nibble crossover point ====

         NIBBLE = ILAENV( 14, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )

         // ==== Accumulate reflections during ttswp?  Use block
         // .    2-by-2 structure during matrix-matrix multiply? ====

         KACC22 = ILAENV( 16, 'ZLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
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

         for (IT = 1; IT <= ITMAX; IT++) { // 70

            // ==== Done when KBOT falls below ILO ====

            if (KBOT.LT.ILO) GO TO 80;

            // ==== Locate active block ====

            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ) == ZERO ) GO TO 20
            } // 10
            K = ILO
            } // 20
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
            if ( NDFL.LT.KEXNW ) {
               NW = MIN( NWUPBD, NWR )
            } else {
               NW = MIN( NWUPBD, 2*NW )
            }
            if ( NW.LT.NWMAX ) {
               if ( NW.GE.NH-1 ) {
                  NW = NH
               } else {
                  KWTOP = KBOT - NW + 1
                  IF( CABS1( H( KWTOP, KWTOP-1 ) ).GT. CABS1( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               }
            }
            if ( NDFL.LT.KEXNW ) {
               NDEC = -1
            } else if ( NDEC.GE.0 .OR. NW.GE.NWUPBD ) {
               NDEC = NDEC + 1
               if (NW-NDEC.LT.2) NDEC = 0;
               NW = NW - NDEC
            }

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

            zlaqr2(WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, IHIZ, Z, LDZ, LS, LD, W, H( KV, 1 ), LDH, NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, WORK, LWORK );

            // ==== Adjust KBOT accounting for new deflations. ====

            KBOT = KBOT - LD

            // ==== KS points to the shifts. ====

            KS = KBOT - LS + 1

            // ==== Skip an expensive QR sweep if there is a (partly
            // .    heuristic) reason to expect that many eigenvalues
            // .    will deflate without it.  Here, the QR sweep is
            // .    skipped if many eigenvalues have just been deflated
            // .    or if the remaining active block is small.

            if ( ( LD == 0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) && ( KBOT- KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) {

               // ==== NS = nominal number of simultaneous shifts.
               // .    This may be lowered (slightly) if ZLAQR2
               // .    did not provide that many shifts. ====

               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )

               // ==== If there have been no deflations
               // .    in a multiple of KEXSH iterations,
               // .    then try exceptional shifts.
               // .    Otherwise use shifts provided by
               // .    ZLAQR2 above or from the eigenvalues
               // .    of a trailing principal submatrix. ====

               if ( MOD( NDFL, KEXSH ) == 0 ) {
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, KS + 1, -2
                     W( I ) = H( I, I ) + WILK1*CABS1( H( I, I-1 ) )
                     W( I-1 ) = W( I )
                  } // 30
               } else {

                  // ==== Got NS/2 or fewer shifts? Use ZLAHQR
                  // .    on a trailing principal submatrix to
                  // .    get more. (Since NS.LE.NSMAX.LE.(N-3)/6,
                  // .    there is enough space below the subdiagonal
                  // .    to fit an NS-by-NS scratch array.) ====

                  if ( KBOT-KS+1.LE.NS / 2 ) {
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     zlacpy('A', NS, NS, H( KS, KS ), LDH, H( KT, 1 ), LDH );
                     zlahqr( false , false , NS, 1, NS, H( KT, 1 ), LDH, W( KS ), 1, 1, ZDUM, 1, INF );
                     KS = KS + INF

                     // ==== In case of a rare QR failure use
                     // .    eigenvalues of the trailing 2-by-2
                     // .    principal submatrix.  Scale to avoid
                     // .    overflows, underflows and subnormals.
                     // .    (The scale factor S can not be zero,
                     // .    because H(KBOT,KBOT-1) is nonzero.) ====

                     if ( KS.GE.KBOT ) {
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
                     }
                  }

                  if ( KBOT-KS+1.GT.NS ) {

                     // ==== Sort the shifts (Helps a little) ====

                     SORTED = false;
                     DO 50 K = KBOT, KS + 1, -1
                        if (SORTED) GO TO 60;
                        SORTED = true;
                        for (I = KS; I <= K - 1; I++) { // 40
                           if ( CABS1( W( I ) ).LT.CABS1( W( I+1 ) ) ) {
                              SORTED = false;
                              SWAP = W( I )
                              W( I ) = W( I+1 )
                              W( I+1 ) = SWAP
                           }
                        } // 40
                     } // 50
                     } // 60
                  }
               }

               // ==== If there are only two shifts, then use
               // .    only one.  ====

               if ( KBOT-KS+1 == 2 ) {
                  if ( CABS1( W( KBOT )-H( KBOT, KBOT ) ).LT. CABS1( W( KBOT-1 )-H( KBOT, KBOT ) ) ) {
                     W( KBOT-1 ) = W( KBOT )
                  } else {
                     W( KBOT ) = W( KBOT-1 )
                  }
               }

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

               zlaqr5(WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS, W( KS ), H, LDH, ILOZ, IHIZ, Z, LDZ, WORK, 3, H( KU, 1 ), LDH, NVE, H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH );
            }

            // ==== Note progress (or the lack of it). ====

            if ( LD.GT.0 ) {
               NDFL = 1
            } else {
               NDFL = NDFL + 1
            }

            // ==== End of main loop ====
         } // 70

         // ==== Iteration limit exceeded.  Set INFO to show where
         // .    the problem occurred and exit. ====

         INFO = KBOT
         } // 80
      }

      // ==== Return the optimal value of LWORK. ====

      WORK( 1 ) = DCMPLX( LWKOPT, 0 )

      // ==== End of ZLAQR4 ====

      }
