      SUBROUTINE ZLATMS( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, KL, KU, PACK, A, LDA, WORK, INFO )

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
      double             D( * );
      COMPLEX*16         A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      double             ONE;
      const              ONE = 1.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      double             TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839D+0 ;
      // ..
      // .. Local Scalars ..
      bool               GIVENS, ILEXTR, ILTEMP, TOPDWN, ZSYM;
      int                I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA, IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2, IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH, JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC, UUB;
      double             ALPHA, ANGLE, REALC, TEMP;
      COMPLEX*16         C, CT, CTEMP, DUMMY, EXTRA, S, ST
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLARND;
      COMPLEX*16         ZLARND
      // EXTERNAL LSAME, DLARND, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLATM1, DSCAL, XERBLA, ZLAGGE, ZLAGHE, ZLAGSY, ZLAROT, ZLARTG, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, COS, DBLE, DCMPLX, DCONJG, MAX, MIN, MOD, SIN
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
         ZSYM = .FALSE.
      } else if ( LSAME( SYM, 'P' ) ) {
         ISYM = 2
         IRSIGN = 0
         ZSYM = .FALSE.
      } else if ( LSAME( SYM, 'S' ) ) {
         ISYM = 2
         IRSIGN = 0
         ZSYM = .TRUE.
      } else if ( LSAME( SYM, 'H' ) ) {
         ISYM = 2
         IRSIGN = 1
         ZSYM = .FALSE.
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
         CALL XERBLA( 'ZLATMS', -INFO )
         RETURN
      }

      // Initialize random number generator

      DO 10 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   10 CONTINUE

      IF( MOD( ISEED( 4 ), 2 ).NE.1 ) ISEED( 4 ) = ISEED( 4 ) + 1

      // 2)      Set up D  if indicated.

              // Compute D according to COND and MODE

      CALL DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, IINFO )
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
         DO 20 I = 2, MNMIN
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   20    CONTINUE

         if ( TEMP.GT.ZERO ) {
            ALPHA = DMAX / TEMP
         } else {
            INFO = 2
            RETURN
         }

         CALL DSCAL( MNMIN, ALPHA, D, 1 )

      }

      CALL ZLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )

      // 3)      Generate Banded Matrix using Givens rotations.
              // Also the special case of UUB=LLB=0

                // Compute Addressing constants to cover all
                // storage formats.  Whether GE, HE, SY, GB, HB, or SB,
                // upper or lower triangle or both,
               t // he (i,j)-th element is in
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
         DO 30 J = 1, MNMIN
            A( ( 1-ISKEW )*J+IOFFST, J ) = DCMPLX( D( J ) )
   30    CONTINUE

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

            DO 40 J = 1, MNMIN
               A( ( 1-ISKEW )*J+IOFFST, J ) = DCMPLX( D( J ) )
   40       CONTINUE

            if ( TOPDWN ) {
               JKL = 0
               DO 70 JKU = 1, UUB

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // Last row actually rotated is M
                  // Last column actually rotated is MIN( M+JKU, N )

                  DO 60 JR = 1, MIN( M+JKU, N ) + JKL - 1
                     EXTRA = CZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )*ZLARND( 5, ISEED )
                     S = SIN( ANGLE )*ZLARND( 5, ISEED )
                     ICOL = MAX( 1, JR-JKL )
                     if ( JR.LT.M ) {
                        IL = MIN( N, JR+JKU ) + 1 - ICOL
                        CALL ZLAROT( .TRUE., JR.GT.JKL, .FALSE., IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, EXTRA, DUMMY )
                     }

                     // Chase "EXTRA" back up

                     IR = JR
                     IC = ICOL
                     DO 50 JCH = JR - JKL, 1, -JKL - JKU
                        if ( IR.LT.M ) {
                           CALL ZLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, REALC, S, DUMMY )
                           DUMMY = ZLARND( 5, ISEED )
                           C = DCONJG( REALC*DUMMY )
                           S = DCONJG( -S*DUMMY )
                        }
                        IROW = MAX( 1, JCH-JKU )
                        IL = IR + 2 - IROW
                        CTEMP = CZERO
                        ILTEMP = JCH.GT.JKU
                        CALL ZLAROT( .FALSE., ILTEMP, .TRUE., IL, C, S, A( IROW-ISKEW*IC+IOFFST, IC ), ILDA, CTEMP, EXTRA )
                        if ( ILTEMP ) {
                           CALL ZLARTG( A( IROW+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), CTEMP, REALC, S, DUMMY )
                           DUMMY = ZLARND( 5, ISEED )
                           C = DCONJG( REALC*DUMMY )
                           S = DCONJG( -S*DUMMY )

                           ICOL = MAX( 1, JCH-JKU-JKL )
                           IL = IC + 2 - ICOL
                           EXTRA = CZERO
                           CALL ZLAROT( .TRUE., JCH.GT.JKU+JKL, .TRUE., IL, C, S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, CTEMP )
                           IC = ICOL
                           IR = IROW
                        }
   50                CONTINUE
   60             CONTINUE
   70          CONTINUE

               JKU = UUB
               DO 100 JKL = 1, LLB

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  DO 90 JC = 1, MIN( N+JKL, M ) + JKU - 1
                     EXTRA = CZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )*ZLARND( 5, ISEED )
                     S = SIN( ANGLE )*ZLARND( 5, ISEED )
                     IROW = MAX( 1, JC-JKU )
                     if ( JC.LT.N ) {
                        IL = MIN( M, JC+JKL ) + 1 - IROW
                        CALL ZLAROT( .FALSE., JC.GT.JKU, .FALSE., IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, EXTRA, DUMMY )
                     }

                     // Chase "EXTRA" back up

                     IC = JC
                     IR = IROW
                     DO 80 JCH = JC - JKU, 1, -JKL - JKU
                        if ( IC.LT.N ) {
                           CALL ZLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, REALC, S, DUMMY )
                           DUMMY = ZLARND( 5, ISEED )
                           C = DCONJG( REALC*DUMMY )
                           S = DCONJG( -S*DUMMY )
                        }
                        ICOL = MAX( 1, JCH-JKL )
                        IL = IC + 2 - ICOL
                        CTEMP = CZERO
                        ILTEMP = JCH.GT.JKL
                        CALL ZLAROT( .TRUE., ILTEMP, .TRUE., IL, C, S, A( IR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, CTEMP, EXTRA )
                        if ( ILTEMP ) {
                           CALL ZLARTG( A( IR+1-ISKEW*( ICOL+1 )+IOFFST, ICOL+1 ), CTEMP, REALC, S, DUMMY )
                           DUMMY = ZLARND( 5, ISEED )
                           C = DCONJG( REALC*DUMMY )
                           S = DCONJG( -S*DUMMY )
                           IROW = MAX( 1, JCH-JKL-JKU )
                           IL = IR + 2 - IROW
                           EXTRA = CZERO
                           CALL ZLAROT( .FALSE., JCH.GT.JKL+JKU, .TRUE., IL, C, S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, CTEMP )
                           IC = ICOL
                           IR = IROW
                        }
   80                CONTINUE
   90             CONTINUE
  100          CONTINUE

            } else {

               // Bottom-Up -- Start at the bottom right.

               JKL = 0
               DO 130 JKU = 1, UUB

                  // Transform from bandwidth JKL, JKU-1 to JKL, JKU

                  // First row actually rotated is M
                  // First column actually rotated is MIN( M+JKU, N )

                  IENDCH = MIN( M, N+JKL ) - 1
                  DO 120 JC = MIN( M+JKU, N ) - 1, 1 - JKL, -1
                     EXTRA = CZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )*ZLARND( 5, ISEED )
                     S = SIN( ANGLE )*ZLARND( 5, ISEED )
                     IROW = MAX( 1, JC-JKU+1 )
                     if ( JC.GT.0 ) {
                        IL = MIN( M, JC+JKL+1 ) + 1 - IROW
                        CALL ZLAROT( .FALSE., .FALSE., JC+JKL.LT.M, IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, DUMMY, EXTRA )
                     }

                     // Chase "EXTRA" back down

                     IC = JC
                     DO 110 JCH = JC + JKL, IENDCH, JKL + JKU
                        ILEXTR = IC.GT.0
                        if ( ILEXTR ) {
                           CALL ZLARTG( A( JCH-ISKEW*IC+IOFFST, IC ), EXTRA, REALC, S, DUMMY )
                           DUMMY = ZLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                        }
                        IC = MAX( 1, IC )
                        ICOL = MIN( N-1, JCH+JKU )
                        ILTEMP = JCH + JKU.LT.N
                        CTEMP = CZERO
                        CALL ZLAROT( .TRUE., ILEXTR, ILTEMP, ICOL+2-IC, C, S, A( JCH-ISKEW*IC+IOFFST, IC ), ILDA, EXTRA, CTEMP )
                        if ( ILTEMP ) {
                           CALL ZLARTG( A( JCH-ISKEW*ICOL+IOFFST, ICOL ), CTEMP, REALC, S, DUMMY )
                           DUMMY = ZLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = CZERO
                           CALL ZLAROT( .FALSE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( JCH-ISKEW*ICOL+IOFFST, ICOL ), ILDA, CTEMP, EXTRA )
                           IC = ICOL
                        }
  110                CONTINUE
  120             CONTINUE
  130          CONTINUE

               JKU = UUB
               DO 160 JKL = 1, LLB

                  // Transform from bandwidth JKL-1, JKU to JKL, JKU

                  // First row actually rotated is MIN( N+JKL, M )
                  // First column actually rotated is N

                  IENDCH = MIN( N, M+JKU ) - 1
                  DO 150 JR = MIN( N+JKL, M ) - 1, 1 - JKU, -1
                     EXTRA = CZERO
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )*ZLARND( 5, ISEED )
                     S = SIN( ANGLE )*ZLARND( 5, ISEED )
                     ICOL = MAX( 1, JR-JKL+1 )
                     if ( JR.GT.0 ) {
                        IL = MIN( N, JR+JKU+1 ) + 1 - ICOL
                        CALL ZLAROT( .TRUE., .FALSE., JR+JKU.LT.N, IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, DUMMY, EXTRA )
                     }

                     // Chase "EXTRA" back down

                     IR = JR
                     DO 140 JCH = JR + JKU, IENDCH, JKL + JKU
                        ILEXTR = IR.GT.0
                        if ( ILEXTR ) {
                           CALL ZLARTG( A( IR-ISKEW*JCH+IOFFST, JCH ), EXTRA, REALC, S, DUMMY )
                           DUMMY = ZLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                        }
                        IR = MAX( 1, IR )
                        IROW = MIN( M-1, JCH+JKL )
                        ILTEMP = JCH + JKL.LT.M
                        CTEMP = CZERO
                        CALL ZLAROT( .FALSE., ILEXTR, ILTEMP, IROW+2-IR, C, S, A( IR-ISKEW*JCH+IOFFST, JCH ), ILDA, EXTRA, CTEMP )
                        if ( ILTEMP ) {
                           CALL ZLARTG( A( IROW-ISKEW*JCH+IOFFST, JCH ), CTEMP, REALC, S, DUMMY )
                           DUMMY = ZLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = CZERO
                           CALL ZLAROT( .TRUE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( IROW-ISKEW*JCH+IOFFST, JCH ), ILDA, CTEMP, EXTRA )
                           IR = IROW
                        }
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE

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

               DO 170 J = 1, MNMIN
                  A( ( 1-ISKEW )*J+IOFFG, J ) = DCMPLX( D( J ) )
  170          CONTINUE

               DO 200 K = 1, UUB
                  DO 190 JC = 1, N - 1
                     IROW = MAX( 1, JC-K )
                     IL = MIN( JC+1, K+2 )
                     EXTRA = CZERO
                     CTEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 )
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )*ZLARND( 5, ISEED )
                     S = SIN( ANGLE )*ZLARND( 5, ISEED )
                     if ( ZSYM ) {
                        CT = C
                        ST = S
                     } else {
                        CTEMP = DCONJG( CTEMP )
                        CT = DCONJG( C )
                        ST = DCONJG( S )
                     }
                     CALL ZLAROT( .FALSE., JC.GT.K, .TRUE., IL, C, S, A( IROW-ISKEW*JC+IOFFG, JC ), ILDA, EXTRA, CTEMP )                      CALL ZLAROT( .TRUE., .TRUE., .FALSE., MIN( K, N-JC )+1, CT, ST, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, CTEMP, DUMMY )

                     // Chase EXTRA back up the matrix

                     ICOL = JC
                     DO 180 JCH = JC - K, 1, -K
                        CALL ZLARTG( A( JCH+1-ISKEW*( ICOL+1 )+IOFFG, ICOL+1 ), EXTRA, REALC, S, DUMMY )
                        DUMMY = ZLARND( 5, ISEED )
                        C = DCONJG( REALC*DUMMY )
                        S = DCONJG( -S*DUMMY )
                        CTEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 )
                        if ( ZSYM ) {
                           CT = C
                           ST = S
                        } else {
                           CTEMP = DCONJG( CTEMP )
                           CT = DCONJG( C )
                           ST = DCONJG( S )
                        }
                        CALL ZLAROT( .TRUE., .TRUE., .TRUE., K+2, C, S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, CTEMP, EXTRA )
                        IROW = MAX( 1, JCH-K )
                        IL = MIN( JCH+1, K+2 )
                        EXTRA = CZERO
                        CALL ZLAROT( .FALSE., JCH.GT.K, .TRUE., IL, CT, ST, A( IROW-ISKEW*JCH+IOFFG, JCH ), ILDA, EXTRA, CTEMP )
                        ICOL = JCH
  180                CONTINUE
  190             CONTINUE
  200          CONTINUE

               // If we need lower triangle, copy from upper. Note that
              t // he order of copying is chosen to work for 'q' -> 'b'

               if ( IPACK.NE.IPACKG .AND. IPACK.NE.3 ) {
                  DO 230 JC = 1, N
                     IROW = IOFFST - ISKEW*JC
                     if ( ZSYM ) {
                        DO 210 JR = JC, MIN( N, JC+UUB )
                           A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
  210                   CONTINUE
                     } else {
                        DO 220 JR = JC, MIN( N, JC+UUB )
                           A( JR+IROW, JC ) = DCONJG( A( JC-ISKEW*JR+ IOFFG, JR ) )
  220                   CONTINUE
                     }
  230             CONTINUE
                  if ( IPACK.EQ.5 ) {
                     DO 250 JC = N - UUB + 1, N
                        DO 240 JR = N + 2 - JC, UUB + 1
                           A( JR, JC ) = CZERO
  240                   CONTINUE
  250                CONTINUE
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

               DO 260 J = 1, MNMIN
                  A( ( 1-ISKEW )*J+IOFFG, J ) = DCMPLX( D( J ) )
  260          CONTINUE

               DO 290 K = 1, UUB
                  DO 280 JC = N - 1, 1, -1
                     IL = MIN( N+1-JC, K+2 )
                     EXTRA = CZERO
                     CTEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC )
                     ANGLE = TWOPI*DLARND( 1, ISEED )
                     C = COS( ANGLE )*ZLARND( 5, ISEED )
                     S = SIN( ANGLE )*ZLARND( 5, ISEED )
                     if ( ZSYM ) {
                        CT = C
                        ST = S
                     } else {
                        CTEMP = DCONJG( CTEMP )
                        CT = DCONJG( C )
                        ST = DCONJG( S )
                     }
                     CALL ZLAROT( .FALSE., .TRUE., N-JC.GT.K, IL, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, CTEMP, EXTRA )
                     ICOL = MAX( 1, JC-K+1 )
                     CALL ZLAROT( .TRUE., .FALSE., .TRUE., JC+2-ICOL, CT, ST, A( JC-ISKEW*ICOL+IOFFG, ICOL ), ILDA, DUMMY, CTEMP )

                     // Chase EXTRA back down the matrix

                     ICOL = JC
                     DO 270 JCH = JC + K, N - 1, K
                        CALL ZLARTG( A( JCH-ISKEW*ICOL+IOFFG, ICOL ), EXTRA, REALC, S, DUMMY )
                        DUMMY = ZLARND( 5, ISEED )
                        C = REALC*DUMMY
                        S = S*DUMMY
                        CTEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH )
                        if ( ZSYM ) {
                           CT = C
                           ST = S
                        } else {
                           CTEMP = DCONJG( CTEMP )
                           CT = DCONJG( C )
                           ST = DCONJG( S )
                        }
                        CALL ZLAROT( .TRUE., .TRUE., .TRUE., K+2, C, S, A( JCH-ISKEW*ICOL+IOFFG, ICOL ), ILDA, EXTRA, CTEMP )
                        IL = MIN( N+1-JCH, K+2 )
                        EXTRA = CZERO
                        CALL ZLAROT( .FALSE., .TRUE., N-JCH.GT.K, IL, CT, ST, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, CTEMP, EXTRA )
                        ICOL = JCH
  270                CONTINUE
  280             CONTINUE
  290          CONTINUE

               // If we need upper triangle, copy from lower. Note that
              t // he order of copying is chosen to work for 'b' -> 'q'

               if ( IPACK.NE.IPACKG .AND. IPACK.NE.4 ) {
                  DO 320 JC = N, 1, -1
                     IROW = IOFFST - ISKEW*JC
                     if ( ZSYM ) {
                        DO 300 JR = JC, MAX( 1, JC-UUB ), -1
                           A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
  300                   CONTINUE
                     } else {
                        DO 310 JR = JC, MAX( 1, JC-UUB ), -1
                           A( JR+IROW, JC ) = DCONJG( A( JC-ISKEW*JR+ IOFFG, JR ) )
  310                   CONTINUE
                     }
  320             CONTINUE
                  if ( IPACK.EQ.6 ) {
                     DO 340 JC = 1, UUB
                        DO 330 JR = 1, UUB + 1 - JC
                           A( JR, JC ) = CZERO
  330                   CONTINUE
  340                CONTINUE
                  }
                  if ( IPACKG.EQ.5 ) {
                     IPACKG = IPACK
                  } else {
                     IPACKG = 0
                  }
               }
            }

            // Ensure that the diagonal is real if Hermitian

            if ( .NOT.ZSYM ) {
               DO 350 JC = 1, N
                  IROW = IOFFST + ( 1-ISKEW )*JC
                  A( IROW, JC ) = DCMPLX( DBLE( A( IROW, JC ) ) )
  350          CONTINUE
            }

         }

      } else {

         // 4)      Generate Banded Matrix by first
                 // Rotating by random Unitary matrices,
                t // hen reducing the bandwidth using Householder
                t // ransformations.

                 // Note: we should get here only if LDA .ge. N

         if ( ISYM.EQ.1 ) {

            // Non-symmetric -- A = U D V

            CALL ZLAGGE( MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO )
         } else {

            // Symmetric -- A = U D U' or
            // Hermitian -- A = U D U*

            if ( ZSYM ) {
               CALL ZLAGSY( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
            } else {
               CALL ZLAGHE( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
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

            DO 370 J = 1, M
               DO 360 I = J + 1, M
                  A( I, J ) = CZERO
  360          CONTINUE
  370       CONTINUE

         } else if ( IPACK.EQ.2 ) {

            // 'L' -- Lower triangular, not packed

            DO 390 J = 2, M
               DO 380 I = 1, J - 1
                  A( I, J ) = CZERO
  380          CONTINUE
  390       CONTINUE

         } else if ( IPACK.EQ.3 ) {

            // 'C' -- Upper triangle packed Columnwise.

            ICOL = 1
            IROW = 0
            DO 410 J = 1, M
               DO 400 I = 1, J
                  IROW = IROW + 1
                  if ( IROW.GT.LDA ) {
                     IROW = 1
                     ICOL = ICOL + 1
                  }
                  A( IROW, ICOL ) = A( I, J )
  400          CONTINUE
  410       CONTINUE

         } else if ( IPACK.EQ.4 ) {

            // 'R' -- Lower triangle packed Columnwise.

            ICOL = 1
            IROW = 0
            DO 430 J = 1, M
               DO 420 I = J, M
                  IROW = IROW + 1
                  if ( IROW.GT.LDA ) {
                     IROW = 1
                     ICOL = ICOL + 1
                  }
                  A( IROW, ICOL ) = A( I, J )
  420          CONTINUE
  430       CONTINUE

         } else if ( IPACK.GE.5 ) {

            // 'B' -- The lower triangle is packed as a band matrix.
            // 'Q' -- The upper triangle is packed as a band matrix.
            // 'Z' -- The whole matrix is packed as a band matrix.

            IF( IPACK.EQ.5 ) UUB = 0             IF( IPACK.EQ.6 ) LLB = 0

            DO 450 J = 1, UUB
               DO 440 I = MIN( J+LLB, M ), 1, -1
                  A( I-J+UUB+1, J ) = A( I, J )
  440          CONTINUE
  450       CONTINUE

            DO 470 J = UUB + 2, N
               DO 460 I = J - UUB, MIN( J+LLB, M )
                  A( I-J+UUB+1, J ) = A( I, J )
  460          CONTINUE
  470       CONTINUE
         }

         // If packed, zero out extraneous elements.

         // Symmetric/Triangular Packed --
         // zero out everything after A(IROW,ICOL)

         if ( IPACK.EQ.3 .OR. IPACK.EQ.4 ) {
            DO 490 JC = ICOL, M
               DO 480 JR = IROW + 1, LDA
                  A( JR, JC ) = CZERO
  480          CONTINUE
               IROW = 0
  490       CONTINUE

         } else if ( IPACK.GE.5 ) {

            // Packed Band --
               // 1st row is now in A( UUB+2-j, j), zero above it
               // m-th row is now in A( M+UUB-j,j), zero below it
               // last non-zero diagonal is now in A( UUB+LLB+1,j ),
                  // zero below it, too.

            IR1 = UUB + LLB + 2
            IR2 = UUB + M + 2
            DO 520 JC = 1, N
               DO 500 JR = 1, UUB + 1 - JC
                  A( JR, JC ) = CZERO
  500          CONTINUE
               DO 510 JR = MAX( 1, MIN( IR1, IR2-JC ) ), LDA
                  A( JR, JC ) = CZERO
  510          CONTINUE
  520       CONTINUE
         }
      }

      RETURN

      // End of ZLATMS

      }
