      SUBROUTINE STGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, LDVL, VR, LDVR, MM, M, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDP, LDS, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      REAL               P( LDP, * ), S( LDS, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..


*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, SAFETY
      const              ZERO = 0.0E+0, ONE = 1.0E+0, SAFETY = 1.0E+2 ;
      // ..
      // .. Local Scalars ..
      bool               COMPL, COMPR, IL2BY2, ILABAD, ILALL, ILBACK, ILBBAD, ILCOMP, ILCPLX, LSA, LSB;
      int                I, IBEG, IEIG, IEND, IHWMNY, IINFO, IM, ISIDE, J, JA, JC, JE, JR, JW, NA, NW;
      REAL               ACOEF, ACOEFA, ANORM, ASCALE, BCOEFA, BCOEFI, BCOEFR, BIG, BIGNUM, BNORM, BSCALE, CIM2A, CIM2B, CIMAGA, CIMAGB, CRE2A, CRE2B, CREALA, CREALB, DMIN, SAFMIN, SALFAR, SBETA, SCALE, SMALL, TEMP, TEMP2, TEMP2I, TEMP2R, ULP, XMAX, XSCALE;
      // ..
      // .. Local Arrays ..
      REAL               BDIAG( 2 ), SUM( 2, 2 ), SUMS( 2, 2 ), SUMP( 2, 2 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SLACPY, SLAG2, SLALN2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and Test the input parameters

      if ( LSAME( HOWMNY, 'A' ) ) {
         IHWMNY = 1
         ILALL = .TRUE.
         ILBACK = .FALSE.
      } else if ( LSAME( HOWMNY, 'S' ) ) {
         IHWMNY = 2
         ILALL = .FALSE.
         ILBACK = .FALSE.
      } else if ( LSAME( HOWMNY, 'B' ) ) {
         IHWMNY = 3
         ILALL = .TRUE.
         ILBACK = .TRUE.
      } else {
         IHWMNY = -1
         ILALL = .TRUE.
      }

      if ( LSAME( SIDE, 'R' ) ) {
         ISIDE = 1
         COMPL = .FALSE.
         COMPR = .TRUE.
      } else if ( LSAME( SIDE, 'L' ) ) {
         ISIDE = 2
         COMPL = .TRUE.
         COMPR = .FALSE.
      } else if ( LSAME( SIDE, 'B' ) ) {
         ISIDE = 3
         COMPL = .TRUE.
         COMPR = .TRUE.
      } else {
         ISIDE = -1
      }

      INFO = 0
      if ( ISIDE.LT.0 ) {
         INFO = -1
      } else if ( IHWMNY.LT.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDS.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDP.LT.MAX( 1, N ) ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('STGEVC', -INFO );
         RETURN
      }

      // Count the number of eigenvectors to be computed

      if ( .NOT.ILALL ) {
         IM = 0
         ILCPLX = .FALSE.
         for (J = 1; J <= N; J++) { // 10
            if ( ILCPLX ) {
               ILCPLX = .FALSE.
               GO TO 10
            }
            if ( J.LT.N ) {
               IF( S( J+1, J ).NE.ZERO ) ILCPLX = .TRUE.
            }
            if ( ILCPLX ) {
               IF( SELECT( J ) .OR. SELECT( J+1 ) ) IM = IM + 2
            } else {
               IF( SELECT( J ) ) IM = IM + 1
            }
   10    CONTINUE
      } else {
         IM = N
      }

      // Check 2-by-2 diagonal blocks of A, B

      ILABAD = .FALSE.
      ILBBAD = .FALSE.
      DO 20 J = 1, N - 1
         if ( S( J+1, J ).NE.ZERO ) {
            IF( P( J, J ).EQ.ZERO .OR. P( J+1, J+1 ).EQ.ZERO .OR. P( J, J+1 ).NE.ZERO )ILBBAD = .TRUE.
            if ( J.LT.N-1 ) {
               IF( S( J+2, J+1 ).NE.ZERO ) ILABAD = .TRUE.
            }
         }
   20 CONTINUE

      if ( ILABAD ) {
         INFO = -5
      } else if ( ILBBAD ) {
         INFO = -7
      } else if ( COMPL .AND. LDVL.LT.N .OR. LDVL.LT.1 ) {
         INFO = -10
      } else if ( COMPR .AND. LDVR.LT.N .OR. LDVR.LT.1 ) {
         INFO = -12
      } else if ( MM.LT.IM ) {
         INFO = -13
      }
      if ( INFO.NE.0 ) {
         xerbla('STGEVC', -INFO );
         RETURN
      }

      // Quick return if possible

      M = IM
      IF( N.EQ.0 ) RETURN

      // Machine Constants

      SAFMIN = SLAMCH( 'Safe minimum' )
      BIG = ONE / SAFMIN
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      SMALL = SAFMIN*N / ULP
      BIG = ONE / SMALL
      BIGNUM = ONE / ( SAFMIN*N )

      // Compute the 1-norm of each column of the strictly upper triangular
      // part (i.e., excluding all elements belonging to the diagonal
      // blocks) of A and B to check for possible overflow in the
      // triangular solver.

      ANORM = ABS( S( 1, 1 ) )
      IF( N.GT.1 ) ANORM = ANORM + ABS( S( 2, 1 ) )
      BNORM = ABS( P( 1, 1 ) )
      WORK( 1 ) = ZERO
      WORK( N+1 ) = ZERO

      for (J = 2; J <= N; J++) { // 50
         TEMP = ZERO
         TEMP2 = ZERO
         if ( S( J, J-1 ).EQ.ZERO ) {
            IEND = J - 1
         } else {
            IEND = J - 2
         }
         for (I = 1; I <= IEND; I++) { // 30
            TEMP = TEMP + ABS( S( I, J ) )
            TEMP2 = TEMP2 + ABS( P( I, J ) )
   30    CONTINUE
         WORK( J ) = TEMP
         WORK( N+J ) = TEMP2
         DO 40 I = IEND + 1, MIN( J+1, N )
            TEMP = TEMP + ABS( S( I, J ) )
            TEMP2 = TEMP2 + ABS( P( I, J ) )
   40    CONTINUE
         ANORM = MAX( ANORM, TEMP )
         BNORM = MAX( BNORM, TEMP2 )
   50 CONTINUE

      ASCALE = ONE / MAX( ANORM, SAFMIN )
      BSCALE = ONE / MAX( BNORM, SAFMIN )

      // Left eigenvectors

      if ( COMPL ) {
         IEIG = 0

         // Main loop over eigenvalues

         ILCPLX = .FALSE.
         for (JE = 1; JE <= N; JE++) { // 220

            // Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
            // (b) this would be the second of a complex pair.
            // Check for complex eigenvalue, so as to be sure of which
            // entry(-ies) of SELECT to look at.

            if ( ILCPLX ) {
               ILCPLX = .FALSE.
               GO TO 220
            }
            NW = 1
            if ( JE.LT.N ) {
               if ( S( JE+1, JE ).NE.ZERO ) {
                  ILCPLX = .TRUE.
                  NW = 2
               }
            }
            if ( ILALL ) {
               ILCOMP = .TRUE.
            } else if ( ILCPLX ) {
               ILCOMP = SELECT( JE ) .OR. SELECT( JE+1 )
            } else {
               ILCOMP = SELECT( JE )
            }
            IF( .NOT.ILCOMP ) GO TO 220

            // Decide if (a) singular pencil, (b) real eigenvalue, or
            // (c) complex eigenvalue.

            if ( .NOT.ILCPLX ) {
               if ( ABS( S( JE, JE ) ).LE.SAFMIN .AND. ABS( P( JE, JE ) ).LE.SAFMIN ) {

                  // Singular matrix pencil -- return unit eigenvector

                  IEIG = IEIG + 1
                  for (JR = 1; JR <= N; JR++) { // 60
                     VL( JR, IEIG ) = ZERO
   60             CONTINUE
                  VL( IEIG, IEIG ) = ONE
                  GO TO 220
               }
            }

            // Clear vector

            DO 70 JR = 1, NW*N
               WORK( 2*N+JR ) = ZERO
   70       CONTINUE
                                                  // T
            // Compute coefficients in  ( a A - b B )  y = 0
               // a  is  ACOEF
               // b  is  BCOEFR + i*BCOEFI

            if ( .NOT.ILCPLX ) {

               // Real eigenvalue

               TEMP = ONE / MAX( ABS( S( JE, JE ) )*ASCALE, ABS( P( JE, JE ) )*BSCALE, SAFMIN )
               SALFAR = ( TEMP*S( JE, JE ) )*ASCALE
               SBETA = ( TEMP*P( JE, JE ) )*BSCALE
               ACOEF = SBETA*ASCALE
               BCOEFR = SALFAR*BSCALE
               BCOEFI = ZERO

               // Scale to avoid underflow

               SCALE = ONE
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEF ).LT.SMALL
               LSB = ABS( SALFAR ).GE.SAFMIN .AND. ABS( BCOEFR ).LT. SMALL                IF( LSA ) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )                IF( LSB ) SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )* MIN( BNORM, BIG ) )
               if ( LSA .OR. LSB ) {
                  SCALE = MIN( SCALE, ONE / ( SAFMIN*MAX( ONE, ABS( ACOEF ), ABS( BCOEFR ) ) ) )
                  if ( LSA ) {
                     ACOEF = ASCALE*( SCALE*SBETA )
                  } else {
                     ACOEF = SCALE*ACOEF
                  }
                  if ( LSB ) {
                     BCOEFR = BSCALE*( SCALE*SALFAR )
                  } else {
                     BCOEFR = SCALE*BCOEFR
                  }
               }
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR )

               // First component is 1

               WORK( 2*N+JE ) = ONE
               XMAX = ONE
            } else {

               // Complex eigenvalue

               slag2(S( JE, JE ), LDS, P( JE, JE ), LDP, SAFMIN*SAFETY, ACOEF, TEMP, BCOEFR, TEMP2, BCOEFI );
               BCOEFI = -BCOEFI
               if ( BCOEFI.EQ.ZERO ) {
                  INFO = JE
                  RETURN
               }

               // Scale to avoid over/underflow

               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               SCALE = ONE
               IF( ACOEFA*ULP.LT.SAFMIN .AND. ACOEFA.GE.SAFMIN ) SCALE = ( SAFMIN / ULP ) / ACOEFA                IF( BCOEFA*ULP.LT.SAFMIN .AND. BCOEFA.GE.SAFMIN ) SCALE = MAX( SCALE, ( SAFMIN / ULP ) / BCOEFA )                IF( SAFMIN*ACOEFA.GT.ASCALE ) SCALE = ASCALE / ( SAFMIN*ACOEFA )                IF( SAFMIN*BCOEFA.GT.BSCALE ) SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) )
               if ( SCALE.NE.ONE ) {
                  ACOEF = SCALE*ACOEF
                  ACOEFA = ABS( ACOEF )
                  BCOEFR = SCALE*BCOEFR
                  BCOEFI = SCALE*BCOEFI
                  BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               }

               // Compute first two components of eigenvector

               TEMP = ACOEF*S( JE+1, JE )
               TEMP2R = ACOEF*S( JE, JE ) - BCOEFR*P( JE, JE )
               TEMP2I = -BCOEFI*P( JE, JE )
               if ( ABS( TEMP ).GT.ABS( TEMP2R )+ABS( TEMP2I ) ) {
                  WORK( 2*N+JE ) = ONE
                  WORK( 3*N+JE ) = ZERO
                  WORK( 2*N+JE+1 ) = -TEMP2R / TEMP
                  WORK( 3*N+JE+1 ) = -TEMP2I / TEMP
               } else {
                  WORK( 2*N+JE+1 ) = ONE
                  WORK( 3*N+JE+1 ) = ZERO
                  TEMP = ACOEF*S( JE, JE+1 )
                  WORK( 2*N+JE ) = ( BCOEFR*P( JE+1, JE+1 )-ACOEF* S( JE+1, JE+1 ) ) / TEMP
                  WORK( 3*N+JE ) = BCOEFI*P( JE+1, JE+1 ) / TEMP
               }
               XMAX = MAX( ABS( WORK( 2*N+JE ) )+ABS( WORK( 3*N+JE ) ), ABS( WORK( 2*N+JE+1 ) )+ABS( WORK( 3*N+JE+1 ) ) )
            }

            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )

                                            // T
            // Triangular solve of  (a A - b B)  y = 0

                                    // T
            // (rowwise in  (a A - b B) , or columnwise in (a A - b B) )

            IL2BY2 = .FALSE.

            DO 160 J = JE + NW, N
               if ( IL2BY2 ) {
                  IL2BY2 = .FALSE.
                  GO TO 160
               }

               NA = 1
               BDIAG( 1 ) = P( J, J )
               if ( J.LT.N ) {
                  if ( S( J+1, J ).NE.ZERO ) {
                     IL2BY2 = .TRUE.
                     BDIAG( 2 ) = P( J+1, J+1 )
                     NA = 2
                  }
               }

               // Check whether scaling is necessary for dot products

               XSCALE = ONE / MAX( ONE, XMAX )
               TEMP = MAX( WORK( J ), WORK( N+J ), ACOEFA*WORK( J )+BCOEFA*WORK( N+J ) )                IF( IL2BY2 ) TEMP = MAX( TEMP, WORK( J+1 ), WORK( N+J+1 ), ACOEFA*WORK( J+1 )+BCOEFA*WORK( N+J+1 ) )
               if ( TEMP.GT.BIGNUM*XSCALE ) {
                  DO 90 JW = 0, NW - 1
                     DO 80 JR = JE, J - 1
                        WORK( ( JW+2 )*N+JR ) = XSCALE* WORK( ( JW+2 )*N+JR )
   80                CONTINUE
   90             CONTINUE
                  XMAX = XMAX*XSCALE
               }

               // Compute dot products

                     // j-1
               // SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
                     // k=je

               // To reduce the op count, this is done as

               // _        j-1                  _        j-1
               // a*conjg( sum  S(k,j)*x(k) ) - b*conjg( sum  P(k,j)*x(k) )
                        // k=je                          k=je

               // which may cause underflow problems if A or B are close
               // to underflow.  (E.g., less than SMALL.)


               for (JW = 1; JW <= NW; JW++) { // 120
                  for (JA = 1; JA <= NA; JA++) { // 110
                     SUMS( JA, JW ) = ZERO
                     SUMP( JA, JW ) = ZERO

                     DO 100 JR = JE, J - 1
                        SUMS( JA, JW ) = SUMS( JA, JW ) + S( JR, J+JA-1 )* WORK( ( JW+1 )*N+JR )                         SUMP( JA, JW ) = SUMP( JA, JW ) + P( JR, J+JA-1 )* WORK( ( JW+1 )*N+JR )
  100                CONTINUE
  110             CONTINUE
  120          CONTINUE

               for (JA = 1; JA <= NA; JA++) { // 130
                  if ( ILCPLX ) {
                     SUM( JA, 1 ) = -ACOEF*SUMS( JA, 1 ) + BCOEFR*SUMP( JA, 1 ) - BCOEFI*SUMP( JA, 2 )                      SUM( JA, 2 ) = -ACOEF*SUMS( JA, 2 ) + BCOEFR*SUMP( JA, 2 ) + BCOEFI*SUMP( JA, 1 )
                  } else {
                     SUM( JA, 1 ) = -ACOEF*SUMS( JA, 1 ) + BCOEFR*SUMP( JA, 1 )
                  }
  130          CONTINUE

                                   // T
               // Solve  ( a A - b B )  y = SUM(,)
               // with scaling and perturbation of the denominator

               slaln2(.TRUE., NA, NW, DMIN, ACOEF, S( J, J ), LDS, BDIAG( 1 ), BDIAG( 2 ), SUM, 2, BCOEFR, BCOEFI, WORK( 2*N+J ), N, SCALE, TEMP, IINFO );
               if ( SCALE.LT.ONE ) {
                  DO 150 JW = 0, NW - 1
                     DO 140 JR = JE, J - 1
                        WORK( ( JW+2 )*N+JR ) = SCALE* WORK( ( JW+2 )*N+JR )
  140                CONTINUE
  150             CONTINUE
                  XMAX = SCALE*XMAX
               }
               XMAX = MAX( XMAX, TEMP )
  160       CONTINUE

            // Copy eigenvector to VL, back transforming if
            // HOWMNY='B'.

            IEIG = IEIG + 1
            if ( ILBACK ) {
               DO 170 JW = 0, NW - 1
                  sgemv('N', N, N+1-JE, ONE, VL( 1, JE ), LDVL, WORK( ( JW+2 )*N+JE ), 1, ZERO, WORK( ( JW+4 )*N+1 ), 1 );
  170          CONTINUE
               slacpy(' ', N, NW, WORK( 4*N+1 ), N, VL( 1, JE ), LDVL );
               IBEG = 1
            } else {
               slacpy(' ', N, NW, WORK( 2*N+1 ), N, VL( 1, IEIG ), LDVL );
               IBEG = JE
            }

            // Scale eigenvector

            XMAX = ZERO
            if ( ILCPLX ) {
               for (J = IBEG; J <= N; J++) { // 180
                  XMAX = MAX( XMAX, ABS( VL( J, IEIG ) )+ ABS( VL( J, IEIG+1 ) ) )
  180          CONTINUE
            } else {
               for (J = IBEG; J <= N; J++) { // 190
                  XMAX = MAX( XMAX, ABS( VL( J, IEIG ) ) )
  190          CONTINUE
            }

            if ( XMAX.GT.SAFMIN ) {
               XSCALE = ONE / XMAX

               DO 210 JW = 0, NW - 1
                  for (JR = IBEG; JR <= N; JR++) { // 200
                     VL( JR, IEIG+JW ) = XSCALE*VL( JR, IEIG+JW )
  200             CONTINUE
  210          CONTINUE
            }
            IEIG = IEIG + NW - 1

  220    CONTINUE
      }

      // Right eigenvectors

      if ( COMPR ) {
         IEIG = IM + 1

         // Main loop over eigenvalues

         ILCPLX = .FALSE.
         DO 500 JE = N, 1, -1

            // Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
            // (b) this would be the second of a complex pair.
            // Check for complex eigenvalue, so as to be sure of which
            // entry(-ies) of SELECT to look at -- if complex, SELECT(JE)
            // or SELECT(JE-1).
            // If this is a complex pair, the 2-by-2 diagonal block
            // corresponding to the eigenvalue is in rows/columns JE-1:JE

            if ( ILCPLX ) {
               ILCPLX = .FALSE.
               GO TO 500
            }
            NW = 1
            if ( JE.GT.1 ) {
               if ( S( JE, JE-1 ).NE.ZERO ) {
                  ILCPLX = .TRUE.
                  NW = 2
               }
            }
            if ( ILALL ) {
               ILCOMP = .TRUE.
            } else if ( ILCPLX ) {
               ILCOMP = SELECT( JE ) .OR. SELECT( JE-1 )
            } else {
               ILCOMP = SELECT( JE )
            }
            IF( .NOT.ILCOMP ) GO TO 500

            // Decide if (a) singular pencil, (b) real eigenvalue, or
            // (c) complex eigenvalue.

            if ( .NOT.ILCPLX ) {
               if ( ABS( S( JE, JE ) ).LE.SAFMIN .AND. ABS( P( JE, JE ) ).LE.SAFMIN ) {

                  // Singular matrix pencil -- unit eigenvector

                  IEIG = IEIG - 1
                  for (JR = 1; JR <= N; JR++) { // 230
                     VR( JR, IEIG ) = ZERO
  230             CONTINUE
                  VR( IEIG, IEIG ) = ONE
                  GO TO 500
               }
            }

            // Clear vector

            DO 250 JW = 0, NW - 1
               for (JR = 1; JR <= N; JR++) { // 240
                  WORK( ( JW+2 )*N+JR ) = ZERO
  240          CONTINUE
  250       CONTINUE

            // Compute coefficients in  ( a A - b B ) x = 0
               // a  is  ACOEF
               // b  is  BCOEFR + i*BCOEFI

            if ( .NOT.ILCPLX ) {

               // Real eigenvalue

               TEMP = ONE / MAX( ABS( S( JE, JE ) )*ASCALE, ABS( P( JE, JE ) )*BSCALE, SAFMIN )
               SALFAR = ( TEMP*S( JE, JE ) )*ASCALE
               SBETA = ( TEMP*P( JE, JE ) )*BSCALE
               ACOEF = SBETA*ASCALE
               BCOEFR = SALFAR*BSCALE
               BCOEFI = ZERO

               // Scale to avoid underflow

               SCALE = ONE
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEF ).LT.SMALL
               LSB = ABS( SALFAR ).GE.SAFMIN .AND. ABS( BCOEFR ).LT. SMALL                IF( LSA ) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )                IF( LSB ) SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )* MIN( BNORM, BIG ) )
               if ( LSA .OR. LSB ) {
                  SCALE = MIN( SCALE, ONE / ( SAFMIN*MAX( ONE, ABS( ACOEF ), ABS( BCOEFR ) ) ) )
                  if ( LSA ) {
                     ACOEF = ASCALE*( SCALE*SBETA )
                  } else {
                     ACOEF = SCALE*ACOEF
                  }
                  if ( LSB ) {
                     BCOEFR = BSCALE*( SCALE*SALFAR )
                  } else {
                     BCOEFR = SCALE*BCOEFR
                  }
               }
               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR )

               // First component is 1

               WORK( 2*N+JE ) = ONE
               XMAX = ONE

               // Compute contribution from column JE of A and B to sum
               // (See "Further Details", above.)

               DO 260 JR = 1, JE - 1
                  WORK( 2*N+JR ) = BCOEFR*P( JR, JE ) - ACOEF*S( JR, JE )
  260          CONTINUE
            } else {

               // Complex eigenvalue

               slag2(S( JE-1, JE-1 ), LDS, P( JE-1, JE-1 ), LDP, SAFMIN*SAFETY, ACOEF, TEMP, BCOEFR, TEMP2, BCOEFI );
               if ( BCOEFI.EQ.ZERO ) {
                  INFO = JE - 1
                  RETURN
               }

               // Scale to avoid over/underflow

               ACOEFA = ABS( ACOEF )
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               SCALE = ONE
               IF( ACOEFA*ULP.LT.SAFMIN .AND. ACOEFA.GE.SAFMIN ) SCALE = ( SAFMIN / ULP ) / ACOEFA                IF( BCOEFA*ULP.LT.SAFMIN .AND. BCOEFA.GE.SAFMIN ) SCALE = MAX( SCALE, ( SAFMIN / ULP ) / BCOEFA )                IF( SAFMIN*ACOEFA.GT.ASCALE ) SCALE = ASCALE / ( SAFMIN*ACOEFA )                IF( SAFMIN*BCOEFA.GT.BSCALE ) SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) )
               if ( SCALE.NE.ONE ) {
                  ACOEF = SCALE*ACOEF
                  ACOEFA = ABS( ACOEF )
                  BCOEFR = SCALE*BCOEFR
                  BCOEFI = SCALE*BCOEFI
                  BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI )
               }

               // Compute first two components of eigenvector
               // and contribution to sums

               TEMP = ACOEF*S( JE, JE-1 )
               TEMP2R = ACOEF*S( JE, JE ) - BCOEFR*P( JE, JE )
               TEMP2I = -BCOEFI*P( JE, JE )
               if ( ABS( TEMP ).GE.ABS( TEMP2R )+ABS( TEMP2I ) ) {
                  WORK( 2*N+JE ) = ONE
                  WORK( 3*N+JE ) = ZERO
                  WORK( 2*N+JE-1 ) = -TEMP2R / TEMP
                  WORK( 3*N+JE-1 ) = -TEMP2I / TEMP
               } else {
                  WORK( 2*N+JE-1 ) = ONE
                  WORK( 3*N+JE-1 ) = ZERO
                  TEMP = ACOEF*S( JE-1, JE )
                  WORK( 2*N+JE ) = ( BCOEFR*P( JE-1, JE-1 )-ACOEF* S( JE-1, JE-1 ) ) / TEMP
                  WORK( 3*N+JE ) = BCOEFI*P( JE-1, JE-1 ) / TEMP
               }

               XMAX = MAX( ABS( WORK( 2*N+JE ) )+ABS( WORK( 3*N+JE ) ), ABS( WORK( 2*N+JE-1 ) )+ABS( WORK( 3*N+JE-1 ) ) )

               // Compute contribution from columns JE and JE-1
               // of A and B to the sums.

               CREALA = ACOEF*WORK( 2*N+JE-1 )
               CIMAGA = ACOEF*WORK( 3*N+JE-1 )
               CREALB = BCOEFR*WORK( 2*N+JE-1 ) - BCOEFI*WORK( 3*N+JE-1 )                CIMAGB = BCOEFI*WORK( 2*N+JE-1 ) + BCOEFR*WORK( 3*N+JE-1 )
               CRE2A = ACOEF*WORK( 2*N+JE )
               CIM2A = ACOEF*WORK( 3*N+JE )
               CRE2B = BCOEFR*WORK( 2*N+JE ) - BCOEFI*WORK( 3*N+JE )
               CIM2B = BCOEFI*WORK( 2*N+JE ) + BCOEFR*WORK( 3*N+JE )
               DO 270 JR = 1, JE - 2
                  WORK( 2*N+JR ) = -CREALA*S( JR, JE-1 ) + CREALB*P( JR, JE-1 ) - CRE2A*S( JR, JE ) + CRE2B*P( JR, JE )                   WORK( 3*N+JR ) = -CIMAGA*S( JR, JE-1 ) + CIMAGB*P( JR, JE-1 ) - CIM2A*S( JR, JE ) + CIM2B*P( JR, JE )
  270          CONTINUE
            }

            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )

            // Columnwise triangular solve of  (a A - b B)  x = 0

            IL2BY2 = .FALSE.
            DO 370 J = JE - NW, 1, -1

               // If a 2-by-2 block, is in position j-1:j, wait until
               // next iteration to process it (when it will be j:j+1)

               if ( .NOT.IL2BY2 .AND. J.GT.1 ) {
                  if ( S( J, J-1 ).NE.ZERO ) {
                     IL2BY2 = .TRUE.
                     GO TO 370
                  }
               }
               BDIAG( 1 ) = P( J, J )
               if ( IL2BY2 ) {
                  NA = 2
                  BDIAG( 2 ) = P( J+1, J+1 )
               } else {
                  NA = 1
               }

               // Compute x(j) (and x(j+1), if 2-by-2 block)

               slaln2(.FALSE., NA, NW, DMIN, ACOEF, S( J, J ), LDS, BDIAG( 1 ), BDIAG( 2 ), WORK( 2*N+J ), N, BCOEFR, BCOEFI, SUM, 2, SCALE, TEMP, IINFO );
               if ( SCALE.LT.ONE ) {

                  DO 290 JW = 0, NW - 1
                     for (JR = 1; JR <= JE; JR++) { // 280
                        WORK( ( JW+2 )*N+JR ) = SCALE* WORK( ( JW+2 )*N+JR )
  280                CONTINUE
  290             CONTINUE
               }
               XMAX = MAX( SCALE*XMAX, TEMP )

               for (JW = 1; JW <= NW; JW++) { // 310
                  for (JA = 1; JA <= NA; JA++) { // 300
                     WORK( ( JW+1 )*N+J+JA-1 ) = SUM( JA, JW )
  300             CONTINUE
  310          CONTINUE

               // w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling

               if ( J.GT.1 ) {

                  // Check whether scaling is necessary for sum.

                  XSCALE = ONE / MAX( ONE, XMAX )
                  TEMP = ACOEFA*WORK( J ) + BCOEFA*WORK( N+J )
                  IF( IL2BY2 ) TEMP = MAX( TEMP, ACOEFA*WORK( J+1 )+BCOEFA* WORK( N+J+1 ) )
                  TEMP = MAX( TEMP, ACOEFA, BCOEFA )
                  if ( TEMP.GT.BIGNUM*XSCALE ) {

                     DO 330 JW = 0, NW - 1
                        for (JR = 1; JR <= JE; JR++) { // 320
                           WORK( ( JW+2 )*N+JR ) = XSCALE* WORK( ( JW+2 )*N+JR )
  320                   CONTINUE
  330                CONTINUE
                     XMAX = XMAX*XSCALE
                  }

                  // Compute the contributions of the off-diagonals of
                  // column j (and j+1, if 2-by-2 block) of A and B to the
                  // sums.


                  for (JA = 1; JA <= NA; JA++) { // 360
                     if ( ILCPLX ) {
                        CREALA = ACOEF*WORK( 2*N+J+JA-1 )
                        CIMAGA = ACOEF*WORK( 3*N+J+JA-1 )
                        CREALB = BCOEFR*WORK( 2*N+J+JA-1 ) - BCOEFI*WORK( 3*N+J+JA-1 )                         CIMAGB = BCOEFI*WORK( 2*N+J+JA-1 ) + BCOEFR*WORK( 3*N+J+JA-1 )
                        DO 340 JR = 1, J - 1
                           WORK( 2*N+JR ) = WORK( 2*N+JR ) - CREALA*S( JR, J+JA-1 ) + CREALB*P( JR, J+JA-1 )                            WORK( 3*N+JR ) = WORK( 3*N+JR ) - CIMAGA*S( JR, J+JA-1 ) + CIMAGB*P( JR, J+JA-1 )
  340                   CONTINUE
                     } else {
                        CREALA = ACOEF*WORK( 2*N+J+JA-1 )
                        CREALB = BCOEFR*WORK( 2*N+J+JA-1 )
                        DO 350 JR = 1, J - 1
                           WORK( 2*N+JR ) = WORK( 2*N+JR ) - CREALA*S( JR, J+JA-1 ) + CREALB*P( JR, J+JA-1 )
  350                   CONTINUE
                     }
  360             CONTINUE
               }

               IL2BY2 = .FALSE.
  370       CONTINUE

            // Copy eigenvector to VR, back transforming if
            // HOWMNY='B'.

            IEIG = IEIG - NW
            if ( ILBACK ) {

               DO 410 JW = 0, NW - 1
                  for (JR = 1; JR <= N; JR++) { // 380
                     WORK( ( JW+4 )*N+JR ) = WORK( ( JW+2 )*N+1 )* VR( JR, 1 )
  380             CONTINUE

                  // A series of compiler directives to defeat
                  // vectorization for the next loop


                  for (JC = 2; JC <= JE; JC++) { // 400
                     for (JR = 1; JR <= N; JR++) { // 390
                        WORK( ( JW+4 )*N+JR ) = WORK( ( JW+4 )*N+JR ) + WORK( ( JW+2 )*N+JC )*VR( JR, JC )
  390                CONTINUE
  400             CONTINUE
  410          CONTINUE

               DO 430 JW = 0, NW - 1
                  for (JR = 1; JR <= N; JR++) { // 420
                     VR( JR, IEIG+JW ) = WORK( ( JW+4 )*N+JR )
  420             CONTINUE
  430          CONTINUE

               IEND = N
            } else {
               DO 450 JW = 0, NW - 1
                  for (JR = 1; JR <= N; JR++) { // 440
                     VR( JR, IEIG+JW ) = WORK( ( JW+2 )*N+JR )
  440             CONTINUE
  450          CONTINUE

               IEND = JE
            }

            // Scale eigenvector

            XMAX = ZERO
            if ( ILCPLX ) {
               for (J = 1; J <= IEND; J++) { // 460
                  XMAX = MAX( XMAX, ABS( VR( J, IEIG ) )+ ABS( VR( J, IEIG+1 ) ) )
  460          CONTINUE
            } else {
               for (J = 1; J <= IEND; J++) { // 470
                  XMAX = MAX( XMAX, ABS( VR( J, IEIG ) ) )
  470          CONTINUE
            }

            if ( XMAX.GT.SAFMIN ) {
               XSCALE = ONE / XMAX
               DO 490 JW = 0, NW - 1
                  for (JR = 1; JR <= IEND; JR++) { // 480
                     VR( JR, IEIG+JW ) = XSCALE*VR( JR, IEIG+JW )
  480             CONTINUE
  490          CONTINUE
            }
  500    CONTINUE
      }

      RETURN

      // End of STGEVC

      }
