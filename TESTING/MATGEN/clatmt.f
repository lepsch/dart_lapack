      SUBROUTINE CLATMT( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RANK, KL, KU, PACK, A, LDA, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               COND, DMAX
      int                INFO, KL, KU, LDA, M, MODE, N, RANK;
      String             DIST, PACK, SYM;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), WORK( * )
      REAL               D( * )
      int                ISEED( 4 );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      REAL               ONE
      const              ONE = 1.0E+0 ;
      COMPLEX            CZERO
      const              CZERO = ( 0.0E+0, 0.0E+0 ) ;
      REAL               TWOPI
      const      TWOPI = 6.28318530717958647692528676655900576839E+0 ;
      // ..
      // .. Local Scalars ..
      COMPLEX            C, CT, CTEMP, DUMMY, EXTRA, S, ST
      REAL               ALPHA, ANGLE, REALC, TEMP
      int                I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA, IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2, IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH, JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC, UUB;
      bool               CSYM, GIVENS, ILEXTR, ILTEMP, TOPDWN;
      // ..
      // .. External Functions ..
      COMPLEX            CLARND
      REAL               SLARND
      bool               LSAME;
      // EXTERNAL CLARND, SLARND, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLAGGE, CLAGHE, CLAGSY, CLAROT, CLARTG, CLASET, SLATM7, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, CONJG, COS, MAX, MIN, MOD, REAL, SIN
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
         CSYM = .FALSE.
      } else if ( LSAME( SYM, 'P' ) ) {
         ISYM = 2
         IRSIGN = 0
         CSYM = .FALSE.
      } else if ( LSAME( SYM, 'S' ) ) {
         ISYM = 2
         IRSIGN = 0
         CSYM = .TRUE.
      } else if ( LSAME( SYM, 'H' ) ) {
         ISYM = 2
         IRSIGN = 1
         CSYM = .FALSE.
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
         xerbla('CLATMT', -INFO );
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

      claset('Full', LDA, N, CZERO, CZERO, A, LDA );

      // 3)      Generate Banded Matrix using Givens rotations.
              // Also the special case of UUB=LLB=0

                // Compute Addressing constants to cover all
                // storage formats.  Whether GE, HE, SY, GB, HB, or SB,
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

      // Diagonal Matrix -- We are done, unless it
      // is to be stored HP/SP/PP/TP (PACK='R' or 'C')

      if ( LLB.EQ.0 .AND. UUB.EQ.0 ) {
         for (J = 1; J <= MNMIN; J++) { // 120
            A( ( 1-ISKEW )*J+IOFFST, J ) = CMPLX( D( J ) )
         } // 120

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

            for (J = 1; J <= MNMIN; J++) { // 130
               A( ( 1-ISKEW )*J+IOFFST, J ) = CMPLX( D( J ) )
            } // 130

            if ( TOPDWN ) {
               JKL = 0
               for (JKU = 1; JKU <= UUB; JKU++) { // 160

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // Last row actually rotated is M
                  // Last column actually rotated is MIN( M+JKU, N )

                  DO 150 JR = 1, MIN( M+JKU, N ) + JKL - 1
                     EXTRA = CZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     ICOL = MAX( 1, JR-JKL )
                     if ( JR.LT.M ) {
                        IL = MIN( N, JR+JKU ) + 1 - ICOL
                        clarot(.TRUE., JR.GT.JKL, .FALSE., IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IR = JR
                     IC = ICOL
                     DO 140 JCH = JR - JKL, 1, -JKL - JKU
                        if ( IR.LT.M ) {
                           clartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED )
                           C = CONJG( REALC*DUMMY )
                           S = CONJG( -S*DUMMY )
                        }
                        IROW = MAX( 1, JCH-JKU )
                        IL = IR + 2 - IROW
                        CTEMP = CZERO
                        ILTEMP = JCH.GT.JKU
                        clarot(.FALSE., ILTEMP, .TRUE., IL, C, S, A( IROW-ISKEW*IC+IOFFST, IC ), ILDA, CTEMP, EXTRA );
                        if ( ILTEMP ) {
                           clartg(A( IROW+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), CTEMP, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED )
                           C = CONJG( REALC*DUMMY )
                           S = CONJG( -S*DUMMY )

                           ICOL = MAX( 1, JCH-JKU-JKL )
                           IL = IC + 2 - ICOL
                           EXTRA = CZERO
                           clarot(.TRUE., JCH.GT.JKU+JKL, .TRUE., IL, C, S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, CTEMP );
                           IC = ICOL
                           IR = IROW
                        }
                     } // 140
                  } // 150
               } // 160

               JKU = UUB
               for (JKL = 1; JKL <= LLB; JKL++) { // 190

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  DO 180 JC = 1, MIN( N+JKL, M ) + JKU - 1
                     EXTRA = CZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     IROW = MAX( 1, JC-JKU )
                     if ( JC.LT.N ) {
                        IL = MIN( M, JC+JKL ) + 1 - IROW
                        clarot(.FALSE., JC.GT.JKU, .FALSE., IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, EXTRA, DUMMY );
                     }

                     // Chase "EXTRA" back up

                     IC = JC
                     IR = IROW
                     DO 170 JCH = JC - JKU, 1, -JKL - JKU
                        if ( IC.LT.N ) {
                           clartg(A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED )
                           C = CONJG( REALC*DUMMY )
                           S = CONJG( -S*DUMMY )
                        }
                        ICOL = MAX( 1, JCH-JKL )
                        IL = IC + 2 - ICOL
                        CTEMP = CZERO
                        ILTEMP = JCH.GT.JKL
                        clarot(.TRUE., ILTEMP, .TRUE., IL, C, S, A( IR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, CTEMP, EXTRA );
                        if ( ILTEMP ) {
                           clartg(A( IR+1-ISKEW*( ICOL+1 )+IOFFST, ICOL+1 ), CTEMP, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED )
                           C = CONJG( REALC*DUMMY )
                           S = CONJG( -S*DUMMY )
                           IROW = MAX( 1, JCH-JKL-JKU )
                           IL = IR + 2 - IROW
                           EXTRA = CZERO
                           clarot(.FALSE., JCH.GT.JKL+JKU, .TRUE., IL, C, S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, CTEMP );
                           IC = ICOL
                           IR = IROW
                        }
                     } // 170
                  } // 180
               } // 190

            } else {

               // Bottom-Up -- Start at the bottom right.

               JKL = 0
               for (JKU = 1; JKU <= UUB; JKU++) { // 220

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // First row actually rotated is M
                  // First column actually rotated is MIN( M+JKU, N )

                  IENDCH = MIN( M, N+JKL ) - 1
                  DO 210 JC = MIN( M+JKU, N ) - 1, 1 - JKL, -1
                     EXTRA = CZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     IROW = MAX( 1, JC-JKU+1 )
                     if ( JC.GT.0 ) {
                        IL = MIN( M, JC+JKL+1 ) + 1 - IROW
                        clarot(.FALSE., .FALSE., JC+JKL.LT.M, IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IC = JC
                     DO 200 JCH = JC + JKL, IENDCH, JKL + JKU
                        ILEXTR = IC.GT.0
                        if ( ILEXTR ) {
                           clartg(A( JCH-ISKEW*IC+IOFFST, IC ), EXTRA, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                        }
                        IC = MAX( 1, IC )
                        ICOL = MIN( N-1, JCH+JKU )
                        ILTEMP = JCH + JKU.LT.N
                        CTEMP = CZERO
                        clarot(.TRUE., ILEXTR, ILTEMP, ICOL+2-IC, C, S, A( JCH-ISKEW*IC+IOFFST, IC ), ILDA, EXTRA, CTEMP );
                        if ( ILTEMP ) {
                           clartg(A( JCH-ISKEW*ICOL+IOFFST, ICOL ), CTEMP, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = CZERO
                           clarot(.FALSE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( JCH-ISKEW*ICOL+IOFFST, ICOL ), ILDA, CTEMP, EXTRA );
                           IC = ICOL
                        }
                     } // 200
                  } // 210
               } // 220

               JKU = UUB
               for (JKL = 1; JKL <= LLB; JKL++) { // 250

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  // First row actually rotated is MIN( N+JKL, M )
                  // First column actually rotated is N

                  IENDCH = MIN( N, M+JKU ) - 1
                  DO 240 JR = MIN( N+JKL, M ) - 1, 1 - JKU, -1
                     EXTRA = CZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     ICOL = MAX( 1, JR-JKL+1 )
                     if ( JR.GT.0 ) {
                        IL = MIN( N, JR+JKU+1 ) + 1 - ICOL
                        clarot(.TRUE., .FALSE., JR+JKU.LT.N, IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, DUMMY, EXTRA );
                     }

                     // Chase "EXTRA" back down

                     IR = JR
                     DO 230 JCH = JR + JKU, IENDCH, JKL + JKU
                        ILEXTR = IR.GT.0
                        if ( ILEXTR ) {
                           clartg(A( IR-ISKEW*JCH+IOFFST, JCH ), EXTRA, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                        }
                        IR = MAX( 1, IR )
                        IROW = MIN( M-1, JCH+JKL )
                        ILTEMP = JCH + JKL.LT.M
                        CTEMP = CZERO
                        clarot(.FALSE., ILEXTR, ILTEMP, IROW+2-IR, C, S, A( IR-ISKEW*JCH+IOFFST, JCH ), ILDA, EXTRA, CTEMP );
                        if ( ILTEMP ) {
                           clartg(A( IROW-ISKEW*JCH+IOFFST, JCH ), CTEMP, REALC, S, DUMMY );
                           DUMMY = CLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = CZERO
                           clarot(.TRUE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( IROW-ISKEW*JCH+IOFFST, JCH ), ILDA, CTEMP, EXTRA );
                           IR = IROW
                        }
                     } // 230
                  } // 240
               } // 250

            }

         } else {

            // Symmetric -- A = U D U'
            // Hermitian -- A = U D U*

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

               for (J = 1; J <= MNMIN; J++) { // 260
                  A( ( 1-ISKEW )*J+IOFFG, J ) = CMPLX( D( J ) )
               } // 260

               for (K = 1; K <= UUB; K++) { // 290
                  for (JC = 1; JC <= N - 1; JC++) { // 280
                     IROW = MAX( 1, JC-K )
                     IL = MIN( JC+1, K+2 )
                     EXTRA = CZERO
                     CTEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 )
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     if ( CSYM ) {
                        CT = C
                        ST = S
                     } else {
                        CTEMP = CONJG( CTEMP )
                        CT = CONJG( C )
                        ST = CONJG( S )
                     }
                     clarot(.FALSE., JC.GT.K, .TRUE., IL, C, S, A( IROW-ISKEW*JC+IOFFG, JC ), ILDA, EXTRA, CTEMP );
                     clarot(.TRUE., .TRUE., .FALSE., MIN( K, N-JC )+1, CT, ST, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, CTEMP, DUMMY );

                     // Chase EXTRA back up the matrix

                     ICOL = JC
                     DO 270 JCH = JC - K, 1, -K
                        clartg(A( JCH+1-ISKEW*( ICOL+1 )+IOFFG, ICOL+1 ), EXTRA, REALC, S, DUMMY );
                        DUMMY = CLARND( 5, ISEED )
                        C = CONJG( REALC*DUMMY )
                        S = CONJG( -S*DUMMY )
                        CTEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 )
                        if ( CSYM ) {
                           CT = C
                           ST = S
                        } else {
                           CTEMP = CONJG( CTEMP )
                           CT = CONJG( C )
                           ST = CONJG( S )
                        }
                        clarot(.TRUE., .TRUE., .TRUE., K+2, C, S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, CTEMP, EXTRA );
                        IROW = MAX( 1, JCH-K )
                        IL = MIN( JCH+1, K+2 )
                        EXTRA = CZERO
                        clarot(.FALSE., JCH.GT.K, .TRUE., IL, CT, ST, A( IROW-ISKEW*JCH+IOFFG, JCH ), ILDA, EXTRA, CTEMP );
                        ICOL = JCH
                     } // 270
                  } // 280
               } // 290

               // If we need lower triangle, copy from upper. Note that
               // the order of copying is chosen to work for 'q' -> 'b'

               if ( IPACK.NE.IPACKG .AND. IPACK.NE.3 ) {
                  for (JC = 1; JC <= N; JC++) { // 320
                     IROW = IOFFST - ISKEW*JC
                     if ( CSYM ) {
                        DO 300 JR = JC, MIN( N, JC+UUB )
                           A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
                        } // 300
                     } else {
                        DO 310 JR = JC, MIN( N, JC+UUB )
                           A( JR+IROW, JC ) = CONJG( A( JC-ISKEW*JR+ IOFFG, JR ) )
                        } // 310
                     }
                  } // 320
                  if ( IPACK.EQ.5 ) {
                     for (JC = N - UUB + 1; JC <= N; JC++) { // 340
                        for (JR = N + 2 - JC; JR <= UUB + 1; JR++) { // 330
                           A( JR, JC ) = CZERO
                        } // 330
                     } // 340
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

               for (J = 1; J <= MNMIN; J++) { // 350
                  A( ( 1-ISKEW )*J+IOFFG, J ) = CMPLX( D( J ) )
               } // 350

               for (K = 1; K <= UUB; K++) { // 380
                  DO 370 JC = N - 1, 1, -1
                     IL = MIN( N+1-JC, K+2 )
                     EXTRA = CZERO
                     CTEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC )
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     if ( CSYM ) {
                        CT = C
                        ST = S
                     } else {
                        CTEMP = CONJG( CTEMP )
                        CT = CONJG( C )
                        ST = CONJG( S )
                     }
                     clarot(.FALSE., .TRUE., N-JC.GT.K, IL, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, CTEMP, EXTRA );
                     ICOL = MAX( 1, JC-K+1 )
                     clarot(.TRUE., .FALSE., .TRUE., JC+2-ICOL, CT, ST, A( JC-ISKEW*ICOL+IOFFG, ICOL ), ILDA, DUMMY, CTEMP );

                     // Chase EXTRA back down the matrix

                     ICOL = JC
                     DO 360 JCH = JC + K, N - 1, K
                        clartg(A( JCH-ISKEW*ICOL+IOFFG, ICOL ), EXTRA, REALC, S, DUMMY );
                        DUMMY = CLARND( 5, ISEED )
                        C = REALC*DUMMY
                        S = S*DUMMY
                        CTEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH )
                        if ( CSYM ) {
                           CT = C
                           ST = S
                        } else {
                           CTEMP = CONJG( CTEMP )
                           CT = CONJG( C )
                           ST = CONJG( S )
                        }
                        clarot(.TRUE., .TRUE., .TRUE., K+2, C, S, A( JCH-ISKEW*ICOL+IOFFG, ICOL ), ILDA, EXTRA, CTEMP );
                        IL = MIN( N+1-JCH, K+2 )
                        EXTRA = CZERO
                        clarot(.FALSE., .TRUE., N-JCH.GT.K, IL, CT, ST, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, CTEMP, EXTRA );
                        ICOL = JCH
                     } // 360
                  } // 370
               } // 380

               // If we need upper triangle, copy from lower. Note that
               // the order of copying is chosen to work for 'b' -> 'q'

               if ( IPACK.NE.IPACKG .AND. IPACK.NE.4 ) {
                  DO 410 JC = N, 1, -1
                     IROW = IOFFST - ISKEW*JC
                     if ( CSYM ) {
                        DO 390 JR = JC, MAX( 1, JC-UUB ), -1
                           A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
                        } // 390
                     } else {
                        DO 400 JR = JC, MAX( 1, JC-UUB ), -1
                           A( JR+IROW, JC ) = CONJG( A( JC-ISKEW*JR+ IOFFG, JR ) )
                        } // 400
                     }
                  } // 410
                  if ( IPACK.EQ.6 ) {
                     for (JC = 1; JC <= UUB; JC++) { // 430
                        for (JR = 1; JR <= UUB + 1 - JC; JR++) { // 420
                           A( JR, JC ) = CZERO
                        } // 420
                     } // 430
                  }
                  if ( IPACKG.EQ.5 ) {
                     IPACKG = IPACK
                  } else {
                     IPACKG = 0
                  }
               }
            }

            // Ensure that the diagonal is real if Hermitian

            if ( .NOT.CSYM ) {
               for (JC = 1; JC <= N; JC++) { // 440
                  IROW = IOFFST + ( 1-ISKEW )*JC
                  A( IROW, JC ) = CMPLX( REAL( A( IROW, JC ) ) )
               } // 440
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

            clagge(MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO );
         } else {

            // Symmetric -- A = U D U' or
            // Hermitian -- A = U D U*

            if ( CSYM ) {
               clagsy(M, LLB, D, A, LDA, ISEED, WORK, IINFO );
            } else {
               claghe(M, LLB, D, A, LDA, ISEED, WORK, IINFO );
            }
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

            for (J = 1; J <= M; J++) { // 460
               for (I = J + 1; I <= M; I++) { // 450
                  A( I, J ) = CZERO
               } // 450
            } // 460

         } else if ( IPACK.EQ.2 ) {

            // 'L' -- Lower triangular, not packed

            for (J = 2; J <= M; J++) { // 480
               for (I = 1; I <= J - 1; I++) { // 470
                  A( I, J ) = CZERO
               } // 470
            } // 480

         } else if ( IPACK.EQ.3 ) {

            // 'C' -- Upper triangle packed Columnwise.

            ICOL = 1
            IROW = 0
            for (J = 1; J <= M; J++) { // 500
               for (I = 1; I <= J; I++) { // 490
                  IROW = IROW + 1
                  if ( IROW.GT.LDA ) {
                     IROW = 1
                     ICOL = ICOL + 1
                  }
                  A( IROW, ICOL ) = A( I, J )
               } // 490
            } // 500

         } else if ( IPACK.EQ.4 ) {

            // 'R' -- Lower triangle packed Columnwise.

            ICOL = 1
            IROW = 0
            for (J = 1; J <= M; J++) { // 520
               for (I = J; I <= M; I++) { // 510
                  IROW = IROW + 1
                  if ( IROW.GT.LDA ) {
                     IROW = 1
                     ICOL = ICOL + 1
                  }
                  A( IROW, ICOL ) = A( I, J )
               } // 510
            } // 520

         } else if ( IPACK.GE.5 ) {

            // 'B' -- The lower triangle is packed as a band matrix.
            // 'Q' -- The upper triangle is packed as a band matrix.
            // 'Z' -- The whole matrix is packed as a band matrix.

            if (IPACK.EQ.5) UUB = 0             IF( IPACK.EQ.6 ) LLB = 0;

            for (J = 1; J <= UUB; J++) { // 540
               DO 530 I = MIN( J+LLB, M ), 1, -1
                  A( I-J+UUB+1, J ) = A( I, J )
               } // 530
            } // 540

            for (J = UUB + 2; J <= N; J++) { // 560
               DO 550 I = J - UUB, MIN( J+LLB, M )
                  A( I-J+UUB+1, J ) = A( I, J )
               } // 550
            } // 560
         }

         // If packed, zero out extraneous elements.

         // Symmetric/Triangular Packed --
         // zero out everything after A(IROW,ICOL)

         if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {
            for (JC = ICOL; JC <= M; JC++) { // 580
               for (JR = IROW + 1; JR <= LDA; JR++) { // 570
                  A( JR, JC ) = CZERO
               } // 570
               IROW = 0
            } // 580

         } else if ( IPACK.GE.5 ) {

            // Packed Band --
               // 1st row is now in A( UUB+2-j, j), zero above it
               // m-th row is now in A( M+UUB-j,j), zero below it
               // last non-zero diagonal is now in A( UUB+LLB+1,j ),
                  // zero below it, too.

            IR1 = UUB + LLB + 2
            IR2 = UUB + M + 2
            for (JC = 1; JC <= N; JC++) { // 610
               for (JR = 1; JR <= UUB + 1 - JC; JR++) { // 590
                  A( JR, JC ) = CZERO
               } // 590
               DO 600 JR = MAX( 1, MIN( IR1, IR2-JC ) ), LDA
                  A( JR, JC ) = CZERO
               } // 600
            } // 610
         }
      }

      RETURN

      // End of CLATMT

      }
