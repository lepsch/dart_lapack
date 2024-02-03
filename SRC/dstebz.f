      SUBROUTINE DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             ORDER, RANGE;
      int                IL, INFO, IU, M, N, NSPLIT;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IBLOCK( * ), ISPLIT( * ), IWORK( * );
      double             D( * ), E( * ), W( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 1.0D0 / TWO ;
      double             FUDGE, RELFAC;
      const              FUDGE = 2.1D0, RELFAC = 2.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               NCNVRG, TOOFEW;
      int                IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO, IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX, ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL, NWU;
      double             ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN, TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH;
      // EXTERNAL LSAME, ILAENV, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAEBZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0

      // Decode RANGE

      if ( LSAME( RANGE, 'A' ) ) {
         IRANGE = 1
      } else if ( LSAME( RANGE, 'V' ) ) {
         IRANGE = 2
      } else if ( LSAME( RANGE, 'I' ) ) {
         IRANGE = 3
      } else {
         IRANGE = 0
      }

      // Decode ORDER

      if ( LSAME( ORDER, 'B' ) ) {
         IORDER = 2
      } else if ( LSAME( ORDER, 'E' ) ) {
         IORDER = 1
      } else {
         IORDER = 0
      }

      // Check for Errors

      if ( IRANGE.LE.0 ) {
         INFO = -1
      } else if ( IORDER.LE.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( IRANGE.EQ.2 ) {
         if ( VL.GE.VU ) INFO = -5       ELSE IF( IRANGE.EQ.3 .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) ) {
         INFO = -6
      } else if ( IRANGE.EQ.3 .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) {
         INFO = -7
      }

      if ( INFO.NE.0 ) {
         xerbla('DSTEBZ', -INFO );
         RETURN
      }

      // Initialize error flags

      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.

      // Quick return if possible

      M = 0
      IF( N.EQ.0 ) RETURN

      // Simplifications:

      IF( IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N ) IRANGE = 1

      // Get machine constants
      // NB is the minimum vector length for vector bisection, or 0
      // if only scalar is to be done.

      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'P' )
      RTOLI = ULP*RELFAC
      NB = ILAENV( 1, 'DSTEBZ', ' ', N, -1, -1, -1 )
      IF( NB.LE.1 ) NB = 0

      // Special Case when N=1

      if ( N.EQ.1 ) {
         NSPLIT = 1
         ISPLIT( 1 ) = 1
         if ( IRANGE.EQ.2 .AND. ( VL.GE.D( 1 ) .OR. VU.LT.D( 1 ) ) ) {
            M = 0
         } else {
            W( 1 ) = D( 1 )
            IBLOCK( 1 ) = 1
            M = 1
         }
         RETURN
      }

      // Compute Splitting Points

      NSPLIT = 1
      WORK( N ) = ZERO
      PIVMIN = ONE

      for (J = 2; J <= N; J++) { // 10
         TMP1 = E( J-1 )**2
         if ( ABS( D( J )*D( J-1 ) )*ULP**2+SAFEMN.GT.TMP1 ) {
            ISPLIT( NSPLIT ) = J - 1
            NSPLIT = NSPLIT + 1
            WORK( J-1 ) = ZERO
         } else {
            WORK( J-1 ) = TMP1
            PIVMIN = MAX( PIVMIN, TMP1 )
         }
      } // 10
      ISPLIT( NSPLIT ) = N
      PIVMIN = PIVMIN*SAFEMN

      // Compute Interval and ATOLI

      if ( IRANGE.EQ.3 ) {

         // RANGE='I': Compute the interval containing eigenvalues
                    // IL through IU.

         // Compute Gershgorin interval for entire (split) matrix
         // and use it as the initial interval

         GU = D( 1 )
         GL = D( 1 )
         TMP1 = ZERO

         DO 20 J = 1, N - 1
            TMP2 = SQRT( WORK( J ) )
            GU = MAX( GU, D( J )+TMP1+TMP2 )
            GL = MIN( GL, D( J )-TMP1-TMP2 )
            TMP1 = TMP2
         } // 20

         GU = MAX( GU, D( N )+TMP1 )
         GL = MIN( GL, D( N )-TMP1 )
         TNORM = MAX( ABS( GL ), ABS( GU ) )
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN

         // Compute Iteration parameters

         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2
         if ( ABSTOL.LE.ZERO ) {
            ATOLI = ULP*TNORM
         } else {
            ATOLI = ABSTOL
         }

         WORK( N+1 ) = GL
         WORK( N+2 ) = GL
         WORK( N+3 ) = GU
         WORK( N+4 ) = GU
         WORK( N+5 ) = GL
         WORK( N+6 ) = GU
         IWORK( 1 ) = -1
         IWORK( 2 ) = -1
         IWORK( 3 ) = N + 1
         IWORK( 4 ) = N + 1
         IWORK( 5 ) = IL - 1
         IWORK( 6 ) = IU

         dlaebz(3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E, WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT, IWORK, W, IBLOCK, IINFO );

         if ( IWORK( 6 ).EQ.IU ) {
            WL = WORK( N+1 )
            WLU = WORK( N+3 )
            NWL = IWORK( 1 )
            WU = WORK( N+4 )
            WUL = WORK( N+2 )
            NWU = IWORK( 4 )
         } else {
            WL = WORK( N+2 )
            WLU = WORK( N+4 )
            NWL = IWORK( 2 )
            WU = WORK( N+3 )
            WUL = WORK( N+1 )
            NWU = IWORK( 3 )
         }

         if ( NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) {
            INFO = 4
            RETURN
         }
      } else {

         // RANGE='A' or 'V' -- Set ATOLI

         TNORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ), ABS( D( N ) )+ABS( E( N-1 ) ) )

         DO 30 J = 2, N - 1
            TNORM = MAX( TNORM, ABS( D( J ) )+ABS( E( J-1 ) )+ ABS( E( J ) ) )
         } // 30

         if ( ABSTOL.LE.ZERO ) {
            ATOLI = ULP*TNORM
         } else {
            ATOLI = ABSTOL
         }

         if ( IRANGE.EQ.2 ) {
            WL = VL
            WU = VU
         } else {
            WL = ZERO
            WU = ZERO
         }
      }

      // Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
      // NWL accumulates the number of eigenvalues .le. WL,
      // NWU accumulates the number of eigenvalues .le. WU

      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0

      for (JB = 1; JB <= NSPLIT; JB++) { // 70
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JB )
         IN = IEND - IOFF

         if ( IN.EQ.1 ) {

            // Special Case -- IN=1

            if ( IRANGE.EQ.1 .OR. WL.GE.D( IBEGIN )-PIVMIN ) NWL = NWL + 1             IF( IRANGE.EQ.1 .OR. WU.GE.D( IBEGIN )-PIVMIN ) NWU = NWU + 1             IF( IRANGE.EQ.1 .OR. ( WL.LT.D( IBEGIN )-PIVMIN .AND. WU.GE. D( IBEGIN )-PIVMIN ) ) {
               M = M + 1
               W( M ) = D( IBEGIN )
               IBLOCK( M ) = JB
            }
         } else {

            // General Case -- IN > 1

            // Compute Gershgorin Interval
            // and use it as the initial interval

            GU = D( IBEGIN )
            GL = D( IBEGIN )
            TMP1 = ZERO

            DO 40 J = IBEGIN, IEND - 1
               TMP2 = ABS( E( J ) )
               GU = MAX( GU, D( J )+TMP1+TMP2 )
               GL = MIN( GL, D( J )-TMP1-TMP2 )
               TMP1 = TMP2
            } // 40

            GU = MAX( GU, D( IEND )+TMP1 )
            GL = MIN( GL, D( IEND )-TMP1 )
            BNORM = MAX( ABS( GL ), ABS( GU ) )
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN

            // Compute ATOLI for the current submatrix

            if ( ABSTOL.LE.ZERO ) {
               ATOLI = ULP*MAX( ABS( GL ), ABS( GU ) )
            } else {
               ATOLI = ABSTOL
            }

            if ( IRANGE.GT.1 ) {
               if ( GU.LT.WL ) {
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GO TO 70
               }
               GL = MAX( GL, WL )
               GU = MIN( GU, WU )
               IF( GL.GE.GU ) GO TO 70
            }

            // Set Up Initial Interval

            WORK( N+1 ) = GL
            WORK( N+IN+1 ) = GU
            dlaebz(1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM, IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO );

            NWL = NWL + IWORK( 1 )
            NWU = NWU + IWORK( IN+1 )
            IWOFF = M - IWORK( 1 )

            // Compute Eigenvalues

            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) / LOG( TWO ) ) + 2             CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ), IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT, IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )

            // Copy Eigenvalues Into W and IBLOCK
            // Use -JB for block number for unconverged eigenvalues.

            for (J = 1; J <= IOUT; J++) { // 60
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )

               // Flag non-convergence.

               if ( J.GT.IOUT-IINFO ) {
                  NCNVRG = .TRUE.
                  IB = -JB
               } else {
                  IB = JB
               }
               DO 50 JE = IWORK( J ) + 1 + IWOFF, IWORK( J+IN ) + IWOFF
                  W( JE ) = TMP1
                  IBLOCK( JE ) = IB
               } // 50
            } // 60

            M = M + IM
         }
      } // 70

      // If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
      // If NWL+1 < IL or NWU > IU, discard extra eigenvalues.

      if ( IRANGE.EQ.3 ) {
         IM = 0
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU

         if ( IDISCL.GT.0 .OR. IDISCU.GT.0 ) {
            for (JE = 1; JE <= M; JE++) { // 80
               if ( W( JE ).LE.WLU .AND. IDISCL.GT.0 ) {
                  IDISCL = IDISCL - 1
               } else if ( W( JE ).GE.WUL .AND. IDISCU.GT.0 ) {
                  IDISCU = IDISCU - 1
               } else {
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               }
            } // 80
            M = IM
         }
         if ( IDISCL.GT.0 .OR. IDISCU.GT.0 ) {

            // Code to deal with effects of bad arithmetic:
            // Some low eigenvalues to be discarded are not in (WL,WLU],
            // or high eigenvalues to be discarded are not in (WUL,WU]
            // so just kill off the smallest IDISCL/largest IDISCU
            // eigenvalues, by simply finding the smallest/largest
            // eigenvalue(s).

            // (If N(w) is monotone non-decreasing, this should never
                // happen.)

            if ( IDISCL.GT.0 ) {
               WKILL = WU
               for (JDISC = 1; JDISC <= IDISCL; JDISC++) { // 100
                  IW = 0
                  for (JE = 1; JE <= M; JE++) { // 90
                     if ( IBLOCK( JE ).NE.0 .AND. ( W( JE ).LT.WKILL .OR. IW.EQ.0 ) ) {
                        IW = JE
                        WKILL = W( JE )
                     }
                  } // 90
                  IBLOCK( IW ) = 0
               } // 100
            }
            if ( IDISCU.GT.0 ) {

               WKILL = WL
               for (JDISC = 1; JDISC <= IDISCU; JDISC++) { // 120
                  IW = 0
                  for (JE = 1; JE <= M; JE++) { // 110
                     if ( IBLOCK( JE ).NE.0 .AND. ( W( JE ).GT.WKILL .OR. IW.EQ.0 ) ) {
                        IW = JE
                        WKILL = W( JE )
                     }
                  } // 110
                  IBLOCK( IW ) = 0
               } // 120
            }
            IM = 0
            for (JE = 1; JE <= M; JE++) { // 130
               if ( IBLOCK( JE ).NE.0 ) {
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               }
            } // 130
            M = IM
         }
         if ( IDISCL.LT.0 .OR. IDISCU.LT.0 ) {
            TOOFEW = .TRUE.
         }
      }

      // If ORDER='B', do nothing -- the eigenvalues are already sorted
         // by block.
      // If ORDER='E', sort the eigenvalues from smallest to largest

      if ( IORDER.EQ.1 .AND. NSPLIT.GT.1 ) {
         DO 150 JE = 1, M - 1
            IE = 0
            TMP1 = W( JE )
            DO 140 J = JE + 1, M
               if ( W( J ).LT.TMP1 ) {
                  IE = J
                  TMP1 = W( J )
               }
            } // 140

            if ( IE.NE.0 ) {
               ITMP1 = IBLOCK( IE )
               W( IE ) = W( JE )
               IBLOCK( IE ) = IBLOCK( JE )
               W( JE ) = TMP1
               IBLOCK( JE ) = ITMP1
            }
         } // 150
      }

      INFO = 0
      IF( NCNVRG ) INFO = INFO + 1       IF( TOOFEW ) INFO = INFO + 2
      RETURN

      // End of DSTEBZ

      }
