      SUBROUTINE DLATMS( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, KL, KU, PACK, A, LDA, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, PACK, SYM;
      int                INFO, KL, KU, LDA, M, MODE, N;
      double             COND, DMAX;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * ), D( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      double             ONE;
      const              ONE = 1.0D0 ;
      double             TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839D+0 ;
      // ..
      // .. Local Scalars ..
      bool               GIVENS, ILEXTR, ILTEMP, TOPDWN;
      int                I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA, IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2, IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH, JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC, UUB;
      double             ALPHA, ANGLE, C, DUMMY, EXTRA, S, TEMP;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLARND;
      // EXTERNAL LSAME, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAGGE, DLAGSY, DLAROT, DLARTG, DLASET, DLATM1, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, COS, DBLE, MAX, MIN, MOD, SIN
      // ..
      // .. Executable Statements ..

      // 1)      Decode and Test the input parameters.
              // Initialize flags & seed.

      INFO = 0

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

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
         IF( DBLE( LLB+UUB ).LT.0.3D0*DBLE( MAX( 1, MR+NC ) ) ) GIVENS = .TRUE.
      } else {
         IF( 2*LLB.LT.M ) GIVENS = .TRUE.
      }
      IF( LDA.LT.M .AND. LDA.GE.MINLDA ) GIVENS = .TRUE.

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
         xerbla('DLATMS', -INFO );
         RETURN
      }

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
      } // 10

      IF( MOD( ISEED( 4 ), 2 ).NE.1 ) ISEED( 4 ) = ISEED( 4 ) + 1

      // 2)      Set up D  if indicated.

              // Compute D according to COND and MODE

      dlatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, IINFO );
      if ( IINFO.NE.0 ) {
         INFO = 1
         RETURN
      }

      // Choose Top-Down if D is (apparently) increasing,
      // Bottom-Up if D is (apparently) decreasing.

      if ( ABS( D( 1 ) ).LE.ABS( D( MNMIN ) ) ) {
         TOPDWN = .TRUE.
      } else {
         TOPDWN = .FALSE.
      }

      if ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) {

         // Scale by DMAX

         TEMP = ABS( D( 1 ) )
         for (I = 2; I <= MNMIN; I++) { // 20
            TEMP = MAX( TEMP, ABS( D( I ) ) )
         } // 20

         if ( TEMP.GT.ZERO ) {
            ALPHA = DMAX / TEMP
         } else {
            INFO = 2
            RETURN
         }

         dscal(MNMIN, ALPHA, D, 1 );

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
      dlaset('Full', LDA, N, ZERO, ZERO, A, LDA );

      // Diagonal Matrix -- We are done, unless it
      // is to be stored SP/PP/TP (PACK='R' or 'C')

      if ( LLB.EQ.0 .AND. UUB.EQ.0 ) {
         dcopy(MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 );
         IF( IPACK.LE.2 .OR. IPACK.GE.5 ) IPACKG = IPACK

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

            dcopy(MNMIN, D, 1, A( 1-ISKEW+IOFFST, 1 ), ILDA+1 );

            if ( TOPDWN ) {
               JKL = 0
               for (JKU = 1; JKU <= UUB; JKU++) { // 50

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // Last row actually rotated is M
                  // Last column actually rotated is MIN( M+JKU, N )

                  DO 40 JR = 1, MIN( M+JKU, N ) + JKL - 1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     ICOL = MAX( 1, JR-JKL )
                     if ( JR.LT.M ) {
                        IL = MIN( N, JR+JKU ) + 1 - ICOL
                        dlarot(.TRUE., JR.GT.JKL, .FALSE., IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IR = JR
                     IC = ICOL
                     DO 30 JCH = JR - JKL, 1, -JKL - JKU
                        if ( IR.LT.M ) {
                           dlartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, C, S, DUMMY );
                        }
                        IROW = MAX( 1, JCH-JKU )
                        IL = IR + 2 - IROW
                        TEMP = ZERO
                        ILTEMP = JCH.GT.JKU
                        dlarot(.FALSE., ILTEMP, .TRUE., IL, C, -S, A( IROW-ISKEW*IC+IOFFST, IC ), ILDA, TEMP, EXTRA );
                        if ( ILTEMP ) {
                           dlartg(A( IROW+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), TEMP, C, S, DUMMY );
                           ICOL = MAX( 1, JCH-JKU-JKL )
                           IL = IC + 2 - ICOL
                           EXTRA = ZERO
                           dlarot(.TRUE., JCH.GT.JKU+JKL, .TRUE., IL, C, -S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, TEMP );
                           IC = ICOL
                           IR = IROW
                        }
                     } // 30
                  } // 40
               } // 50

               JKU = UUB
               for (JKL = 1; JKL <= LLB; JKL++) { // 80

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  DO 70 JC = 1, MIN( N+JKL, M ) + JKU - 1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     IROW = MAX( 1, JC-JKU )
                     if ( JC.LT.N ) {
                        IL = MIN( M, JC+JKL ) + 1 - IROW
                        dlarot(.FALSE., JC.GT.JKU, .FALSE., IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IC = JC
                     IR = IROW
                     DO 60 JCH = JC - JKU, 1, -JKL - JKU
                        if ( IC.LT.N ) {
                           dlartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, C, S, DUMMY );
                        }
                        ICOL = MAX( 1, JCH-JKL )
                        IL = IC + 2 - ICOL
                        TEMP = ZERO
                        ILTEMP = JCH.GT.JKL
                        dlarot(.TRUE., ILTEMP, .TRUE., IL, C, -S, A( IR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, TEMP, EXTRA );
                        if ( ILTEMP ) {
                           dlartg(A( IR+1-ISKEW*( ICOL+1 )+IOFFST, ICOL+1 ), TEMP, C, S, DUMMY );
                           IROW = MAX( 1, JCH-JKL-JKU )
                           IL = IR + 2 - IROW
                           EXTRA = ZERO
                           dlarot(.FALSE., JCH.GT.JKL+JKU, .TRUE., IL, C, -S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, TEMP );
                           IC = ICOL
                           IR = IROW
                        }
                     } // 60
                  } // 70
               } // 80

            } else {

               // Bottom-Up -- Start at the bottom right.

               JKL = 0
               for (JKU = 1; JKU <= UUB; JKU++) { // 110

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // First row actually rotated is M
                  // First column actually rotated is MIN( M+JKU, N )

                  IENDCH = MIN( M, N+JKL ) - 1
                  DO 100 JC = MIN( M+JKU, N ) - 1, 1 - JKL, -1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     IROW = MAX( 1, JC-JKU+1 )
                     if ( JC.GT.0 ) {
                        IL = MIN( M, JC+JKL+1 ) + 1 - IROW
                        dlarot(.FALSE., .FALSE., JC+JKL.LT.M, IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IC = JC
                     DO 90 JCH = JC + JKL, IENDCH, JKL + JKU
                        ILEXTR = IC.GT.0
                        if ( ILEXTR ) {
                           dlartg(A( JCH-ISKEW*IC+IOFFST, IC ), EXTRA, C, S, DUMMY );
                        }
                        IC = MAX( 1, IC )
                        ICOL = MIN( N-1, JCH+JKU )
                        ILTEMP = JCH + JKU.LT.N
                        TEMP = ZERO
                        dlarot(.TRUE., ILEXTR, ILTEMP, ICOL+2-IC, C, S, A( JCH-ISKEW*IC+IOFFST, IC ), ILDA, EXTRA, TEMP );
                        if ( ILTEMP ) {
                           dlartg(A( JCH-ISKEW*ICOL+IOFFST, ICOL ), TEMP, C, S, DUMMY );
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = ZERO
                           dlarot(.FALSE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( JCH-ISKEW*ICOL+IOFFST, ICOL ), ILDA, TEMP, EXTRA );
                           IC = ICOL
                        }
                     } // 90
                  } // 100
               } // 110

               JKU = UUB
               for (JKL = 1; JKL <= LLB; JKL++) { // 140

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  // First row actually rotated is MIN( N+JKL, M )
                  // First column actually rotated is N

                  IENDCH = MIN( N, M+JKU ) - 1
                  DO 130 JR = MIN( N+JKL, M ) - 1, 1 - JKU, -1
                     EXTRA = ZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     ICOL = MAX( 1, JR-JKL+1 )
                     if ( JR.GT.0 ) {
                        IL = MIN( N, JR+JKU+1 ) + 1 - ICOL
                        dlarot(.TRUE., .FALSE., JR+JKU.LT.N, IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IR = JR
                     DO 120 JCH = JR + JKU, IENDCH, JKL + JKU
                        ILEXTR = IR.GT.0
                        if ( ILEXTR ) {
                           dlartg(A( IR-ISKEW*JCH+IOFFST, JCH ), EXTRA, C, S, DUMMY );
                        }
                        IR = MAX( 1, IR )
                        IROW = MIN( M-1, JCH+JKL )
                        ILTEMP = JCH + JKL.LT.M
                        TEMP = ZERO
                        dlarot(.FALSE., ILEXTR, ILTEMP, IROW+2-IR, C, S, A( IR-ISKEW*JCH+IOFFST, JCH ), ILDA, EXTRA, TEMP );
                        if ( ILTEMP ) {
                           dlartg(A( IROW-ISKEW*JCH+IOFFST, JCH ), TEMP, C, S, DUMMY );
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = ZERO
                           dlarot(.TRUE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( IROW-ISKEW*JCH+IOFFST, JCH ), ILDA, TEMP, EXTRA );
                           IR = IROW
                        }
                     } // 120
                  } // 130
               } // 140
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
               dcopy(MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 );

               for (K = 1; K <= UUB; K++) { // 170
                  DO 160 JC = 1, N - 1
                     IROW = MAX( 1, JC-K )
                     IL = MIN( JC+1, K+2 )
                     EXTRA = ZERO
                     TEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 )
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = SIN( ANGLE )
                     dlarot(.FALSE., JC.GT.K, .TRUE., IL, C, S, A( IROW-ISKEW*JC+IOFFG, JC ), ILDA, EXTRA, TEMP )                      CALL DLAROT( .TRUE., .TRUE., .FALSE., MIN( K, N-JC )+1, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, TEMP, DUMMY );

                     // Chase EXTRA back up the matrix

                     ICOL = JC
                     DO 150 JCH = JC - K, 1, -K
                        dlartg(A( JCH+1-ISKEW*( ICOL+1 )+IOFFG, ICOL+1 ), EXTRA, C, S, DUMMY );
                        TEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 )
                        dlarot(.TRUE., .TRUE., .TRUE., K+2, C, -S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, TEMP, EXTRA );
                        IROW = MAX( 1, JCH-K )
                        IL = MIN( JCH+1, K+2 )
                        EXTRA = ZERO
                        dlarot(.FALSE., JCH.GT.K, .TRUE., IL, C, -S, A( IROW-ISKEW*JCH+IOFFG, JCH ), ILDA, EXTRA, TEMP );
                        ICOL = JCH
                     } // 150
                  } // 160
               } // 170

               // If we need lower triangle, copy from upper. Note that
               // the order of copying is chosen to work for 'q' -> 'b'

               if ( IPACK.NE.IPACKG .AND. IPACK.NE.3 ) {
                  for (JC = 1; JC <= N; JC++) { // 190
                     IROW = IOFFST - ISKEW*JC
                     DO 180 JR = JC, MIN( N, JC+UUB )
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
                     } // 180
                  } // 190
                  if ( IPACK.EQ.5 ) {
                     DO 210 JC = N - UUB + 1, N
                        DO 200 JR = N + 2 - JC, UUB + 1
                           A( JR, JC ) = ZERO
                        } // 200
                     } // 210
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
                  IF( IPACK.EQ.6 ) IOFFG = 1
               } else {
                  IPACKG = 2
               }
               dcopy(MNMIN, D, 1, A( 1-ISKEW+IOFFG, 1 ), ILDA+1 );

               for (K = 1; K <= UUB; K++) { // 240
                  DO 230 JC = N - 1, 1, -1
                     IL = MIN( N+1-JC, K+2 )
                     EXTRA = ZERO
                     TEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC )
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )
                     S = -SIN( ANGLE )
                     dlarot(.FALSE., .TRUE., N-JC.GT.K, IL, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, TEMP, EXTRA );
                     ICOL = MAX( 1, JC-K+1 )
                     dlarot(.TRUE., .FALSE., .TRUE., JC+2-ICOL, C, S, A( JC-ISKEW*ICOL+IOFFG, ICOL ), ILDA, DUMMY, TEMP );

                     // Chase EXTRA back down the matrix

                     ICOL = JC
                     DO 220 JCH = JC + K, N - 1, K
                        dlartg(A( JCH-ISKEW*ICOL+IOFFG, ICOL ), EXTRA, C, S, DUMMY );
                        TEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH )
                        dlarot(.TRUE., .TRUE., .TRUE., K+2, C, S, A( JCH-ISKEW*ICOL+IOFFG, ICOL ), ILDA, EXTRA, TEMP );
                        IL = MIN( N+1-JCH, K+2 )
                        EXTRA = ZERO
                        dlarot(.FALSE., .TRUE., N-JCH.GT.K, IL, C, S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, TEMP, EXTRA );
                        ICOL = JCH
                     } // 220
                  } // 230
               } // 240

               // If we need upper triangle, copy from lower. Note that
               // the order of copying is chosen to work for 'b' -> 'q'

               if ( IPACK.NE.IPACKG .AND. IPACK.NE.4 ) {
                  DO 260 JC = N, 1, -1
                     IROW = IOFFST - ISKEW*JC
                     DO 250 JR = JC, MAX( 1, JC-UUB ), -1
                        A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
                     } // 250
                  } // 260
                  if ( IPACK.EQ.6 ) {
                     for (JC = 1; JC <= UUB; JC++) { // 280
                        DO 270 JR = 1, UUB + 1 - JC
                           A( JR, JC ) = ZERO
                        } // 270
                     } // 280
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

            dlagge(MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO );
         } else {

            // Symmetric -- A = U D U'

            dlagsy(M, LLB, D, A, LDA, ISEED, WORK, IINFO );

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

            for (J = 1; J <= M; J++) { // 300
               DO 290 I = J + 1, M
                  A( I, J ) = ZERO
               } // 290
            } // 300

         } else if ( IPACK.EQ.2 ) {

            // 'L' -- Lower triangular, not packed

            for (J = 2; J <= M; J++) { // 320
               DO 310 I = 1, J - 1
                  A( I, J ) = ZERO
               } // 310
            } // 320

         } else if ( IPACK.EQ.3 ) {

            // 'C' -- Upper triangle packed Columnwise.

            ICOL = 1
            IROW = 0
            for (J = 1; J <= M; J++) { // 340
               for (I = 1; I <= J; I++) { // 330
                  IROW = IROW + 1
                  if ( IROW.GT.LDA ) {
                     IROW = 1
                     ICOL = ICOL + 1
                  }
                  A( IROW, ICOL ) = A( I, J )
               } // 330
            } // 340

         } else if ( IPACK.EQ.4 ) {

            // 'R' -- Lower triangle packed Columnwise.

            ICOL = 1
            IROW = 0
            for (J = 1; J <= M; J++) { // 360
               for (I = J; I <= M; I++) { // 350
                  IROW = IROW + 1
                  if ( IROW.GT.LDA ) {
                     IROW = 1
                     ICOL = ICOL + 1
                  }
                  A( IROW, ICOL ) = A( I, J )
               } // 350
            } // 360

         } else if ( IPACK.GE.5 ) {

            // 'B' -- The lower triangle is packed as a band matrix.
            // 'Q' -- The upper triangle is packed as a band matrix.
            // 'Z' -- The whole matrix is packed as a band matrix.

            IF( IPACK.EQ.5 ) UUB = 0             IF( IPACK.EQ.6 ) LLB = 0

            for (J = 1; J <= UUB; J++) { // 380
               DO 370 I = MIN( J+LLB, M ), 1, -1
                  A( I-J+UUB+1, J ) = A( I, J )
               } // 370
            } // 380

            DO 400 J = UUB + 2, N
               DO 390 I = J - UUB, MIN( J+LLB, M )
                  A( I-J+UUB+1, J ) = A( I, J )
               } // 390
            } // 400
         }

         // If packed, zero out extraneous elements.

         // Symmetric/Triangular Packed --
         // zero out everything after A(IROW,ICOL)

         if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {
            for (JC = ICOL; JC <= M; JC++) { // 420
               DO 410 JR = IROW + 1, LDA
                  A( JR, JC ) = ZERO
               } // 410
               IROW = 0
            } // 420

         } else if ( IPACK.GE.5 ) {

            // Packed Band --
               // 1st row is now in A( UUB+2-j, j), zero above it
               // m-th row is now in A( M+UUB-j,j), zero below it
               // last non-zero diagonal is now in A( UUB+LLB+1,j ),
                  // zero below it, too.

            IR1 = UUB + LLB + 2
            IR2 = UUB + M + 2
            for (JC = 1; JC <= N; JC++) { // 450
               DO 430 JR = 1, UUB + 1 - JC
                  A( JR, JC ) = ZERO
               } // 430
               DO 440 JR = MAX( 1, MIN( IR1, IR2-JC ) ), LDA
                  A( JR, JC ) = ZERO
               } // 440
            } // 450
         }
      }

      RETURN

      // End of DLATMS

      }
