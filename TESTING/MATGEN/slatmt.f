      SUBROUTINE SLATMT( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RANK, KL, KU, PACK, A, LDA, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               COND, DMAX
      int                INFO, KL, KU, LDA, M, MODE, N, RANK;
      String             DIST, PACK, SYM;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), WORK( * )
      int                ISEED( 4 );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      REAL               ONE
      const              ONE = 1.0E0 ;
      REAL               TWOPI
      const      TWOPI = 6.28318530717958647692528676655900576839E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               ALPHA, ANGLE, C, DUMMY, EXTRA, S, TEMP
      int                I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA, IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2, IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH, JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC, UUB;
      bool               GIVENS, ILEXTR, ILTEMP, TOPDWN;
      // ..
      // .. External Functions ..
      REAL               SLARND
      bool               LSAME;
      // EXTERNAL SLARND, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLATM7, SCOPY, SLAGGE, SLAGSY, SLAROT, SLARTG, SLASET, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, COS, MAX, MIN, MOD, REAL, SIN
      // ..
      // .. Executable Statements ..

      // 1)      Decode and Test the input parameters.
              // Initialize flags & seed.

      INFO = 0

      // Quick return if possible

      if (M.EQ.0 .OR. N.EQ.0) RETURN;

      // Decode DIST

      if ( LSAME( DIST, 'U' ) ) {
         IDIST = 1
      } else if ( LSAME( DIST, 'S' ) ) {
         IDIST = 2
      } else if ( LSAME( DIST, 'N' ) ) {
         IDIST = 3
      } else {
         IDIST = -1
      }

      // Decode SYM

      if ( LSAME( SYM, 'N' ) ) {
         ISYM = 1
         IRSIGN = 0
      } else if ( LSAME( SYM, 'P' ) ) {
         ISYM = 2
         IRSIGN = 0
      } else if ( LSAME( SYM, 'S' ) ) {
         ISYM = 2
         IRSIGN = 1
      } else if ( LSAME( SYM, 'H' ) ) {
         ISYM = 2
         IRSIGN = 1
      } else {
         ISYM = -1
      }

      // Decode PACK

      ISYMPK = 0
      if ( LSAME( PACK, 'N' ) ) {
         IPACK = 0
      } else if ( LSAME( PACK, 'U' ) ) {
         IPACK = 1
         ISYMPK = 1
      } else if ( LSAME( PACK, 'L' ) ) {
         IPACK = 2
         ISYMPK = 1
      } else if ( LSAME( PACK, 'C' ) ) {
         IPACK = 3
         ISYMPK = 2
      } else if ( LSAME( PACK, 'R' ) ) {
         IPACK = 4
         ISYMPK = 3
      } else if ( LSAME( PACK, 'B' ) ) {
         IPACK = 5
         ISYMPK = 3
      } else if ( LSAME( PACK, 'Q' ) ) {
         IPACK = 6
         ISYMPK = 2
      } else if ( LSAME( PACK, 'Z' ) ) {
         IPACK = 7
      } else {
         IPACK = -1
      }

      // Set certain internal parameters

      MNMIN = MIN( M, N )
      LLB = MIN( KL, M-1 )
      UUB = MIN( KU, N-1 )
      MR = MIN( M, N+LLB )
      NC = MIN( N, M+UUB )

      if ( IPACK.EQ.5 .OR. IPACK.EQ.6 ) {
         MINLDA = UUB + 1
      } else if ( IPACK.EQ.7 ) {
         MINLDA = LLB + UUB + 1
      } else {
         MINLDA = M
      }

      // Use Givens rotation method if bandwidth small enough,
      // or if LDA is too small to store the matrix unpacked.

      GIVENS = .FALSE.
      if ( ISYM.EQ.1 ) {
         IF( REAL( LLB+UUB ).LT.0.3*REAL( MAX( 1, MR+NC ) ) ) GIVENS = .TRUE.
      } else {
         if (2*LLB.LT.M) GIVENS = .TRUE.;
      }
      if (LDA.LT.M .AND. LDA.GE.MINLDA) GIVENS = .TRUE.;

      // Set INFO if an error

      if ( M.LT.0 ) {
         INFO = -1
      } else if ( M.NE.N .AND. ISYM.NE.1 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( IDIST.EQ.-1 ) {
         INFO = -3
      } else if ( ISYM.EQ.-1 ) {
         INFO = -5
      } else if ( ABS( MODE ).GT.6 ) {
         INFO = -7
      } else if ( ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) .AND. COND.LT.ONE ) {
         INFO = -8
      } else if ( KL.LT.0 ) {
         INFO = -10
      } else if ( KU.LT.0 .OR. ( ISYM.NE.1 .AND. KL.NE.KU ) ) {
         INFO = -11
      } else if ( IPACK.EQ.-1 .OR. ( ISYMPK.EQ.1 .AND. ISYM.EQ.1 ) .OR. ( ISYMPK.EQ.2 .AND. ISYM.EQ.1 .AND. KL.GT.0 ) .OR. ( ISYMPK.EQ.3 .AND. ISYM.EQ.1 .AND. KU.GT.0 ) .OR. ( ISYMPK.NE.0 .AND. M.NE.N ) ) {
         INFO = -12
      } else if ( LDA.LT.MAX( 1, MINLDA ) ) {
         INFO = -14
      }

      if ( INFO.NE.0 ) {
         xerbla('SLATMT', -INFO );
         RETURN
      }

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 100
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
      } // 100

      IF( MOD( ISEED( 4 ), 2 ).NE.1 ) ISEED( 4 ) = ISEED( 4 ) + 1

      // 2)      Set up D  if indicated.

              // Compute D according to COND and MODE

      slatm7(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, RANK, IINFO );
      if ( IINFO.NE.0 ) {
         INFO = 1
         RETURN
      }

      // Choose Top-Down if D is (apparently) increasing,
      // Bottom-Up if D is (apparently) decreasing.

      if ( ABS( D( 1 ) ).LE.ABS( D( RANK ) ) ) {
         TOPDWN = .TRUE.
      } else {
         TOPDWN = .FALSE.
      }

      if ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) {

         // Scale by DMAX

         TEMP = ABS( D( 1 ) )
         for (I = 2; I <= RANK; I++) { // 110
            TEMP = MAX( TEMP, ABS( D( I ) ) )
         } // 110

         if ( TEMP.GT.ZERO ) {
            ALPHA = DMAX / TEMP
         } else {
            INFO = 2
            RETURN
         }

         sscal(RANK, ALPHA, D, 1 );

      }

      // 3)      Generate Banded Matrix using Givens rotations.
              // Also the special case of UUB=LLB=0

                // Compute Addressing constants to cover all
                // storage formats.  Whether GE, SY, GB, or SB,
                // upper or lower triangle or both,
                // the (i,j)-th element is in
                // A( i - ISKEW*j + IOFFST, j )

      if ( IPACK.GT.4 ) {
         ILDA = LDA - 1
         ISKEW = 1
         if ( IPACK.GT.5 ) {
            IOFFST = UUB + 1
         } else {
            IOFFST = 1
         }
      } else {
         ILDA = LDA
         ISKEW = 0
         IOFFST = 0
      }

      // IPACKG is the format that the matrix is generated in. If this is
      // different from IPACK, then the matrix must be repacked at the
      // end.  It also signals how to compute the norm, for scaling.

      IPACKG = 0
      slaset('Full', LDA, N, ZERO, ZERO, A, LDA );

      // Diagonal Matrix -- We are done, unless it
      // is to be stored SP/PP/TP (PACK='R' or 'C')

      if ( LLB.EQ.0 .AND. UUB.EQ.0 ) {
         scopy(MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 );
         if (IPACK.LE.2 .OR. IPACK.GE.5) IPACKG = IPACK;

      } else if ( GIVENS ) {

         // Check whether to use Givens rotations,
         // Householder transformations, or nothing.

         if ( ISYM.EQ.1 ) {

            // Non-symmetric -- A = U D V

            if ( IPACK.GT.4 ) {
               IPACKG = IPACK
            } else {
               IPACKG = 0
            }

            scopy(MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 );

            if ( TOPDWN ) {
               JKL = 0
               for (JKU = 1; JKU <= UUB; JKU++) { // 140

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // Last row actually rotated is M
                  // Last column actually rotated is MIN( M+JKU, N )

                  DO 130 JR = 1, MIN( M+JKU, N ) + JKL - 1
                     EXTRA = ZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     ICOL = MAX( 1, JR-JKL )
                     if ( JR.LT.M ) {
                        IL = MIN( N, JR+JKU ) + 1 - ICOL
                        slarot(.TRUE., JR.GT.JKL, .FALSE., IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IR = JR
                     IC = ICOL
                     DO 120 JCH = JR - JKL, 1, -JKL - JKU
                        if ( IR.LT.M ) {
                           slartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, C, S, DUMMY );
                        }
                        IROW = MAX( 1, JCH-JKU )
                        IL = IR + 2 - IROW
                        TEMP = ZERO
                        ILTEMP = JCH.GT.JKU
                        slarot(.FALSE., ILTEMP, .TRUE., IL, C, -S, A( IROW-ISKEW*IC+IOFFST, IC ), ILDA, TEMP, EXTRA );
                        if ( ILTEMP ) {
                           slartg(A( IROW+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), TEMP, C, S, DUMMY );
                           ICOL = MAX( 1, JCH-JKU-JKL )
                           IL = IC + 2 - ICOL
                           EXTRA = ZERO
                           slarot(.TRUE., JCH.GT.JKU+JKL, .TRUE., IL, C, -S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, TEMP );
                           IC = ICOL
                           IR = IROW
                        }
                     } // 120
                  } // 130
               } // 140

               JKU = UUB
               for (JKL = 1; JKL <= LLB; JKL++) { // 170

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  DO 160 JC = 1, MIN( N+JKL, M ) + JKU - 1
                     EXTRA = ZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     IROW = MAX( 1, JC-JKU )
                     if ( JC.LT.N ) {
                        IL = MIN( M, JC+JKL ) + 1 - IROW
                        slarot(.FALSE., JC.GT.JKU, .FALSE., IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IC = JC
                     IR = IROW
                     DO 150 JCH = JC - JKU, 1, -JKL - JKU
                        if ( IC.LT.N ) {
                           slartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, C, S, DUMMY );
                        }
                        ICOL = MAX( 1, JCH-JKL )
                        IL = IC + 2 - ICOL
                        TEMP = ZERO
                        ILTEMP = JCH.GT.JKL
                        slarot(.TRUE., ILTEMP, .TRUE., IL, C, -S, A( IR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, TEMP, EXTRA );
                        if ( ILTEMP ) {
                           slartg(A( IR+1-ISKEW*( ICOL+1 )+IOFFST, ICOL+1 ), TEMP, C, S, DUMMY );
                           IROW = MAX( 1, JCH-JKL-JKU )
                           IL = IR + 2 - IROW
                           EXTRA = ZERO
                           slarot(.FALSE., JCH.GT.JKL+JKU, .TRUE., IL, C, -S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, TEMP );
                           IC = ICOL
                           IR = IROW
                        }
                     } // 150
                  } // 160
               } // 170

            } else {

               // Bottom-Up -- Start at the bottom right.

               JKL = 0
               for (JKU = 1; JKU <= UUB; JKU++) { // 200

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // First row actually rotated is M
                  // First column actually rotated is MIN( M+JKU, N )

                  IENDCH = MIN( M, N+JKL ) - 1
                  DO 190 JC = MIN( M+JKU, N ) - 1, 1 - JKL, -1
                     EXTRA = ZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     IROW = MAX( 1, JC-JKU+1 )
                     if ( JC.GT.0 ) {
                        IL = MIN( M, JC+JKL+1 ) + 1 - IROW
                        slarot(.FALSE., .FALSE., JC+JKL.LT.M, IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IC = JC
                     DO 180 JCH = JC + JKL, IENDCH, JKL + JKU
                        ILEXTR = IC.GT.0
                        if ( ILEXTR ) {
                           slartg(A( JCH-ISKEW*IC+IOFFST, IC ), EXTRA, C, S, DUMMY );
                        }
                        IC = MAX( 1, IC )
                        ICOL = MIN( N-1, JCH+JKU )
                        ILTEMP = JCH + JKU.LT.N
                        TEMP = ZERO
                        slarot(.TRUE., ILEXTR, ILTEMP, ICOL+2-IC, C, S, A( JCH-ISKEW*IC+IOFFST, IC ), ILDA, EXTRA, TEMP );
                        if ( ILTEMP ) {
                           slartg(A( JCH-ISKEW*ICOL+IOFFST, ICOL ), TEMP, C, S, DUMMY );
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = ZERO
                           slarot(.FALSE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( JCH-ISKEW*ICOL+IOFFST, ICOL ), ILDA, TEMP, EXTRA );
                           IC = ICOL
                        }
                     } // 180
                  } // 190
               } // 200

               JKU = UUB
               for (JKL = 1; JKL <= LLB; JKL++) { // 230

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  // First row actually rotated is MIN( N+JKL, M )
                  // First column actually rotated is N

                  IENDCH = MIN( N, M+JKU ) - 1
                  DO 220 JR = MIN( N+JKL, M ) - 1, 1 - JKU, -1
                     EXTRA = ZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     ICOL = MAX( 1, JR-JKL+1 )
                     if ( JR.GT.0 ) {
                        IL = MIN( N, JR+JKU+1 ) + 1 - ICOL
                        slarot(.TRUE., .FALSE., JR+JKU.LT.N, IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IR = JR
                     DO 210 JCH = JR + JKU, IENDCH, JKL + JKU
                        ILEXTR = IR.GT.0
                        if ( ILEXTR ) {
                           slartg(A( IR-ISKEW*JCH+IOFFST, JCH ), EXTRA, C, S, DUMMY );
                        }
                        IR = MAX( 1, IR )
                        IROW = MIN( M-1, JCH+JKL )
                        ILTEMP = JCH + JKL.LT.M
                        TEMP = ZERO
                        slarot(.FALSE., ILEXTR, ILTEMP, IROW+2-IR, C, S, A( IR-ISKEW*JCH+IOFFST, JCH ), ILDA, EXTRA, TEMP );
                        if ( ILTEMP ) {
                           slartg(A( IROW-ISKEW*JCH+IOFFST, JCH ), TEMP, C, S, DUMMY );
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = ZERO
                           slarot(.TRUE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( IROW-ISKEW*JCH+IOFFST, JCH ), ILDA, TEMP, EXTRA );
                           IR = IROW
                        }
                     } // 210
                  } // 220
               } // 230
            }

         } else {

            // Symmetric -- A = U D U'

            IPACKG = IPACK
            IOFFG = IOFFST

            if ( TOPDWN ) {

               // Top-Down -- Generate Upper triangle only

               if ( IPACK.GE.5 ) {
                  IPACKG = 6
                  IOFFG = UUB + 1
               } else {
                  IPACKG = 1
               }
               scopy(MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 );

               for (K = 1; K <= UUB; K++) { // 260
                  for (JC = 1; JC <= N - 1; JC++) { // 250
                     IROW = MAX( 1, JC-K )
                     IL = MIN( JC+1, K+2 )
                     EXTRA = ZERO
                     TEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 )
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     slarot(.FALSE., JC.GT.K, .TRUE., IL, C, S, A( IROW-ISKEW*JC+IOFFG, JC ), ILDA, EXTRA, TEMP );
                     slarot(.TRUE., .TRUE., .FALSE., MIN( K, N-JC )+1, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, TEMP, DUMMY );

                     // Chase EXTRA back up the matrix

                     ICOL = JC
                     DO 240 JCH = JC - K, 1, -K
                        slartg(A( JCH+1-ISKEW*( ICOL+1 )+IOFFG, ICOL+1 ), EXTRA, C, S, DUMMY );
                        TEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 )
                        slarot(.TRUE., .TRUE., .TRUE., K+2, C, -S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, TEMP, EXTRA );
                        IROW = MAX( 1, JCH-K )
                        IL = MIN( JCH+1, K+2 )
                        EXTRA = ZERO
                        slarot(.FALSE., JCH.GT.K, .TRUE., IL, C, -S, A( IROW-ISKEW*JCH+IOFFG, JCH ), ILDA, EXTRA, TEMP );
                        ICOL = JCH
                     } // 240
                  } // 250
               } // 260

               // If we need lower triangle, copy from upper. Note that
               // the order of copying is chosen to work for 'q' -> 'b'

               if ( IPACK.NE.IPACKG .AND. IPACK.NE.3 ) {
                  for (JC = 1; JC <= N; JC++) { // 280
                     IROW = IOFFST - ISKEW*JC
                     DO 270 JR = JC, MIN( N, JC+UUB )
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
                     } // 270
                  } // 280
                  if ( IPACK.EQ.5 ) {
                     for (JC = N - UUB + 1; JC <= N; JC++) { // 300
                        for (JR = N + 2 - JC; JR <= UUB + 1; JR++) { // 290
                           A( JR, JC ) = ZERO
                        } // 290
                     } // 300
                  }
                  if ( IPACKG.EQ.6 ) {
                     IPACKG = IPACK
                  } else {
                     IPACKG = 0
                  }
               }
            } else {

               // Bottom-Up -- Generate Lower triangle only

               if ( IPACK.GE.5 ) {
                  IPACKG = 5
                  if (IPACK.EQ.6) IOFFG = 1;
               } else {
                  IPACKG = 2
               }
               scopy(MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 );

               for (K = 1; K <= UUB; K++) { // 330
                  DO 320 JC = N - 1, 1, -1
                     IL = MIN( N+1-JC, K+2 )
                     EXTRA = ZERO
                     TEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC )
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = -SIN( ANGLE )
                     slarot(.FALSE., .TRUE., N-JC.GT.K, IL, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, TEMP, EXTRA );
                     ICOL = MAX( 1, JC-K+1 )
                     slarot(.TRUE., .FALSE., .TRUE., JC+2-ICOL, C, S, A( JC-ISKEW*ICOL+IOFFG, ICOL ), ILDA, DUMMY, TEMP );

                     // Chase EXTRA back down the matrix

                     ICOL = JC
                     DO 310 JCH = JC + K, N - 1, K
                        slartg(A( JCH-ISKEW*ICOL+IOFFG, ICOL ), EXTRA, C, S, DUMMY );
                        TEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH )
                        slarot(.TRUE., .TRUE., .TRUE., K+2, C, S, A( JCH-ISKEW*ICOL+IOFFG, ICOL ), ILDA, EXTRA, TEMP );
                        IL = MIN( N+1-JCH, K+2 )
                        EXTRA = ZERO
                        slarot(.FALSE., .TRUE., N-JCH.GT.K, IL, C, S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, TEMP, EXTRA );
                        ICOL = JCH
                     } // 310
                  } // 320
               } // 330

               // If we need upper triangle, copy from lower. Note that
               // the order of copying is chosen to work for 'b' -> 'q'

               if ( IPACK.NE.IPACKG .AND. IPACK.NE.4 ) {
                  DO 350 JC = N, 1, -1
                     IROW = IOFFST - ISKEW*JC
                     DO 340 JR = JC, MAX( 1, JC-UUB ), -1
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
                     } // 340
                  } // 350
                  if ( IPACK.EQ.6 ) {
                     for (JC = 1; JC <= UUB; JC++) { // 370
                        for (JR = 1; JR <= UUB + 1 - JC; JR++) { // 360
                           A( JR, JC ) = ZERO
                        } // 360
                     } // 370
                  }
                  if ( IPACKG.EQ.5 ) {
                     IPACKG = IPACK
                  } else {
                     IPACKG = 0
                  }
               }
            }
         }

      } else {

         // 4)      Generate Banded Matrix by first
                 // Rotating by random Unitary matrices,
                 // then reducing the bandwidth using Householder
                 // transformations.

                 // Note: we should get here only if LDA .ge. N

         if ( ISYM.EQ.1 ) {

            // Non-symmetric -- A = U D V

            slagge(MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO );
         } else {

            // Symmetric -- A = U D U'

            slagsy(M, LLB, D, A, LDA, ISEED, WORK, IINFO );

         }
         if ( IINFO.NE.0 ) {
            INFO = 3
            RETURN
         }
      }

      // 5)      Pack the matrix

      if ( IPACK.NE.IPACKG ) {
         if ( IPACK.EQ.1 ) {

            // 'U' -- Upper triangular, not packed

            for (J = 1; J <= M; J++) { // 390
               for (I = J + 1; I <= M; I++) { // 380
                  A( I, J ) = ZERO
               } // 380
            } // 390

         } else if ( IPACK.EQ.2 ) {

            // 'L' -- Lower triangular, not packed

            for (J = 2; J <= M; J++) { // 410
               for (I = 1; I <= J - 1; I++) { // 400
                  A( I, J ) = ZERO
               } // 400
            } // 410

         } else if ( IPACK.EQ.3 ) {

            // 'C' -- Upper triangle packed Columnwise.

            ICOL = 1
            IROW = 0
            for (J = 1; J <= M; J++) { // 430
               for (I = 1; I <= J; I++) { // 420
                  IROW = IROW + 1
                  if ( IROW.GT.LDA ) {
                     IROW = 1
                     ICOL = ICOL + 1
                  }
                  A( IROW, ICOL ) = A( I, J )
               } // 420
            } // 430

         } else if ( IPACK.EQ.4 ) {

            // 'R' -- Lower triangle packed Columnwise.

            ICOL = 1
            IROW = 0
            for (J = 1; J <= M; J++) { // 450
               for (I = J; I <= M; I++) { // 440
                  IROW = IROW + 1
                  if ( IROW.GT.LDA ) {
                     IROW = 1
                     ICOL = ICOL + 1
                  }
                  A( IROW, ICOL ) = A( I, J )
               } // 440
            } // 450

         } else if ( IPACK.GE.5 ) {

            // 'B' -- The lower triangle is packed as a band matrix.
            // 'Q' -- The upper triangle is packed as a band matrix.
            // 'Z' -- The whole matrix is packed as a band matrix.

            if (IPACK.EQ.5) UUB = 0             IF( IPACK.EQ.6 ) LLB = 0;

            for (J = 1; J <= UUB; J++) { // 470
               DO 460 I = MIN( J+LLB, M ), 1, -1
                  A( I-J+UUB+1, J ) = A( I, J )
               } // 460
            } // 470

            for (J = UUB + 2; J <= N; J++) { // 490
               DO 480 I = J - UUB, MIN( J+LLB, M )
                  A( I-J+UUB+1, J ) = A( I, J )
               } // 480
            } // 490
         }

         // If packed, zero out extraneous elements.

         // Symmetric/Triangular Packed --
         // zero out everything after A(IROW,ICOL)

         if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {
            for (JC = ICOL; JC <= M; JC++) { // 510
               for (JR = IROW + 1; JR <= LDA; JR++) { // 500
                  A( JR, JC ) = ZERO
               } // 500
               IROW = 0
            } // 510

         } else if ( IPACK.GE.5 ) {

            // Packed Band --
               // 1st row is now in A( UUB+2-j, j), zero above it
               // m-th row is now in A( M+UUB-j,j), zero below it
               // last non-zero diagonal is now in A( UUB+LLB+1,j ),
                  // zero below it, too.

            IR1 = UUB + LLB + 2
            IR2 = UUB + M + 2
            for (JC = 1; JC <= N; JC++) { // 540
               for (JR = 1; JR <= UUB + 1 - JC; JR++) { // 520
                  A( JR, JC ) = ZERO
               } // 520
               DO 530 JR = MAX( 1, MIN( IR1, IR2-JC ) ), LDA
                  A( JR, JC ) = ZERO
               } // 530
            } // 540
         }
      }

      RETURN

      // End of SLATMT

      }
