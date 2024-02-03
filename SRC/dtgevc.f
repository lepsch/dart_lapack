      SUBROUTINE DTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, LDVL, VR, LDVR, MM, M, WORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDP, LDS, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      double             P( LDP, * ), S( LDS, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..


// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, SAFETY;
      const              ZERO = 0.0, ONE = 1.0, SAFETY = 1.0e+2 ;
      // ..
      // .. Local Scalars ..
      bool               COMPL, COMPR, IL2BY2, ILABAD, ILALL, ILBACK, ILBBAD, ILCOMP, ILCPLX, LSA, LSB;
      int                I, IBEG, IEIG, IEND, IHWMNY, IINFO, IM, ISIDE, J, JA, JC, JE, JR, JW, NA, NW;
      double             ACOEF, ACOEFA, ANORM, ASCALE, BCOEFA, BCOEFI, BCOEFR, BIG, BIGNUM, BNORM, BSCALE, CIM2A, CIM2B, CIMAGA, CIMAGB, CRE2A, CRE2B, CREALA, CREALB, DMIN, SAFMIN, SALFAR, SBETA, SCALE, SMALL, TEMP, TEMP2, TEMP2I, TEMP2R, ULP, XMAX, XSCALE;
      // ..
      // .. Local Arrays ..
      double             BDIAG( 2 ), SUM( 2, 2 ), SUMS( 2, 2 ), SUMP( 2, 2 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DLACPY, DLAG2, DLALN2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and Test the input parameters

      if ( LSAME( HOWMNY, 'A' ) ) {
         IHWMNY = 1;
         ILALL = true;
         ILBACK = false;
      } else if ( LSAME( HOWMNY, 'S' ) ) {
         IHWMNY = 2;
         ILALL = false;
         ILBACK = false;
      } else if ( LSAME( HOWMNY, 'B' ) ) {
         IHWMNY = 3;
         ILALL = true;
         ILBACK = true;
      } else {
         IHWMNY = -1;
         ILALL = true;
      }

      if ( LSAME( SIDE, 'R' ) ) {
         ISIDE = 1;
         COMPL = false;
         COMPR = true;
      } else if ( LSAME( SIDE, 'L' ) ) {
         ISIDE = 2;
         COMPL = true;
         COMPR = false;
      } else if ( LSAME( SIDE, 'B' ) ) {
         ISIDE = 3;
         COMPL = true;
         COMPR = true;
      } else {
         ISIDE = -1;
      }

      INFO = 0;
      if ( ISIDE < 0 ) {
         INFO = -1;
      } else if ( IHWMNY < 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDS < MAX( 1, N ) ) {
         INFO = -6;
      } else if ( LDP < MAX( 1, N ) ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('DTGEVC', -INFO );
         return;
      }

      // Count the number of eigenvectors to be computed

      if ( !ILALL ) {
         IM = 0;
         ILCPLX = false;
         for (J = 1; J <= N; J++) { // 10
            if ( ILCPLX ) {
               ILCPLX = false;
               GO TO 10;
            }
            if ( J < N ) {
               if( S( J+1, J ) != ZERO ) ILCPLX = true;
            }
            if ( ILCPLX ) {
               if( SELECT( J ) || SELECT( J+1 ) ) IM = IM + 2;
            } else {
               if( SELECT( J ) ) IM = IM + 1;
            }
         } // 10
      } else {
         IM = N;
      }

      // Check 2-by-2 diagonal blocks of A, B

      ILABAD = false;
      ILBBAD = false;
      for (J = 1; J <= N - 1; J++) { // 20
         if ( S( J+1, J ) != ZERO ) {
            if( P( J, J ) == ZERO || P( J+1, J+1 ) == ZERO || P( J, J+1 ) != ZERO )ILBBAD = true;
            if ( J < N-1 ) {
               if( S( J+2, J+1 ) != ZERO ) ILABAD = true;
            }
         }
      } // 20

      if ( ILABAD ) {
         INFO = -5;
      } else if ( ILBBAD ) {
         INFO = -7;
      } else if ( COMPL && LDVL < N || LDVL < 1 ) {
         INFO = -10;
      } else if ( COMPR && LDVR < N || LDVR < 1 ) {
         INFO = -12;
      } else if ( MM < IM ) {
         INFO = -13;
      }
      if ( INFO != 0 ) {
         xerbla('DTGEVC', -INFO );
         return;
      }

      // Quick return if possible

      M = IM;
      if (N == 0) RETURN;

      // Machine Constants

      SAFMIN = DLAMCH( 'Safe minimum' );
      BIG = ONE / SAFMIN;
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );
      SMALL = SAFMIN*N / ULP;
      BIG = ONE / SMALL;
      BIGNUM = ONE / ( SAFMIN*N );

      // Compute the 1-norm of each column of the strictly upper triangular
      // part (i.e., excluding all elements belonging to the diagonal
      // blocks) of A and B to check for possible overflow in the
      // triangular solver.

      ANORM = ABS( S( 1, 1 ) );
      if (N > 1) ANORM = ANORM + ABS( S( 2, 1 ) );
      BNORM = ABS( P( 1, 1 ) );
      WORK( 1 ) = ZERO;
      WORK( N+1 ) = ZERO;

      for (J = 2; J <= N; J++) { // 50
         TEMP = ZERO;
         TEMP2 = ZERO;
         if ( S( J, J-1 ) == ZERO ) {
            IEND = J - 1;
         } else {
            IEND = J - 2;
         }
         for (I = 1; I <= IEND; I++) { // 30
            TEMP = TEMP + ABS( S( I, J ) );
            TEMP2 = TEMP2 + ABS( P( I, J ) );
         } // 30
         WORK( J ) = TEMP;
         WORK( N+J ) = TEMP2;
         DO 40 I = IEND + 1, MIN( J+1, N );
            TEMP = TEMP + ABS( S( I, J ) );
            TEMP2 = TEMP2 + ABS( P( I, J ) );
         } // 40
         ANORM = MAX( ANORM, TEMP );
         BNORM = MAX( BNORM, TEMP2 );
      } // 50

      ASCALE = ONE / MAX( ANORM, SAFMIN );
      BSCALE = ONE / MAX( BNORM, SAFMIN );

      // Left eigenvectors

      if ( COMPL ) {
         IEIG = 0;

         // Main loop over eigenvalues

         ILCPLX = false;
         for (JE = 1; JE <= N; JE++) { // 220

            // Skip this iteration if (a) HOWMNY='S' and SELECT= false , or
            // (b) this would be the second of a complex pair.
            // Check for complex eigenvalue, so as to be sure of which
            // entry(-ies) of SELECT to look at.

            if ( ILCPLX ) {
               ILCPLX = false;
               GO TO 220;
            }
            NW = 1;
            if ( JE < N ) {
               if ( S( JE+1, JE ) != ZERO ) {
                  ILCPLX = true;
                  NW = 2;
               }
            }
            if ( ILALL ) {
               ILCOMP = true;
            } else if ( ILCPLX ) {
               ILCOMP = SELECT( JE ) || SELECT( JE+1 );
            } else {
               ILCOMP = SELECT( JE );
            }
            if ( !ILCOMP) GO TO 220;

            // Decide if (a) singular pencil, (b) real eigenvalue, or
            // (c) complex eigenvalue.

            if ( !ILCPLX ) {
               if ( ABS( S( JE, JE ) ) <= SAFMIN && ABS( P( JE, JE ) ) <= SAFMIN ) {

                  // Singular matrix pencil -- return unit eigenvector

                  IEIG = IEIG + 1;
                  for (JR = 1; JR <= N; JR++) { // 60
                     VL( JR, IEIG ) = ZERO;
                  } // 60
                  VL( IEIG, IEIG ) = ONE;
                  GO TO 220;
               }
            }

            // Clear vector

            for (JR = 1; JR <= NW*N; JR++) { // 70
               WORK( 2*N+JR ) = ZERO;
            } // 70
                                                  // T
            // Compute coefficients in  ( a A - b B )  y = 0
               // a  is  ACOEF
               // b  is  BCOEFR + i*BCOEFI

            if ( !ILCPLX ) {

               // Real eigenvalue

               TEMP = ONE / MAX( ABS( S( JE, JE ) )*ASCALE, ABS( P( JE, JE ) )*BSCALE, SAFMIN );
               SALFAR = ( TEMP*S( JE, JE ) )*ASCALE;
               SBETA = ( TEMP*P( JE, JE ) )*BSCALE;
               ACOEF = SBETA*ASCALE;
               BCOEFR = SALFAR*BSCALE;
               BCOEFI = ZERO;

               // Scale to avoid underflow

               SCALE = ONE;
               LSA = ABS( SBETA ) >= SAFMIN && ABS( ACOEF ) < SMALL;
               LSB = ABS( SALFAR ) >= SAFMIN && ABS( BCOEFR ) < SMALL                IF( LSA ) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )                IF( LSB ) SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )* MIN( BNORM, BIG ) );
               if ( LSA || LSB ) {
                  SCALE = MIN( SCALE, ONE / ( SAFMIN*MAX( ONE, ABS( ACOEF ), ABS( BCOEFR ) ) ) );
                  if ( LSA ) {
                     ACOEF = ASCALE*( SCALE*SBETA );
                  } else {
                     ACOEF = SCALE*ACOEF;
                  }
                  if ( LSB ) {
                     BCOEFR = BSCALE*( SCALE*SALFAR );
                  } else {
                     BCOEFR = SCALE*BCOEFR;
                  }
               }
               ACOEFA = ABS( ACOEF );
               BCOEFA = ABS( BCOEFR );

               // First component is 1

               WORK( 2*N+JE ) = ONE;
               XMAX = ONE;
            } else {

               // Complex eigenvalue

               dlag2(S( JE, JE ), LDS, P( JE, JE ), LDP, SAFMIN*SAFETY, ACOEF, TEMP, BCOEFR, TEMP2, BCOEFI );
               BCOEFI = -BCOEFI;
               if ( BCOEFI == ZERO ) {
                  INFO = JE;
                  return;
               }

               // Scale to avoid over/underflow

               ACOEFA = ABS( ACOEF );
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI );
               SCALE = ONE;
               if (ACOEFA*ULP < SAFMIN && ACOEFA >= SAFMIN) SCALE = ( SAFMIN / ULP ) / ACOEFA;
               if( BCOEFA*ULP < SAFMIN && BCOEFA >= SAFMIN ) SCALE = MAX( SCALE, ( SAFMIN / ULP ) / BCOEFA );
               if( SAFMIN*ACOEFA > ASCALE ) SCALE = ASCALE / ( SAFMIN*ACOEFA );
               IF( SAFMIN*BCOEFA > BSCALE ) SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) );
               if ( SCALE != ONE ) {
                  ACOEF = SCALE*ACOEF;
                  ACOEFA = ABS( ACOEF );
                  BCOEFR = SCALE*BCOEFR;
                  BCOEFI = SCALE*BCOEFI;
                  BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI );
               }

               // Compute first two components of eigenvector

               TEMP = ACOEF*S( JE+1, JE );
               TEMP2R = ACOEF*S( JE, JE ) - BCOEFR*P( JE, JE );
               TEMP2I = -BCOEFI*P( JE, JE );
               if ( ABS( TEMP ) > ABS( TEMP2R )+ABS( TEMP2I ) ) {
                  WORK( 2*N+JE ) = ONE;
                  WORK( 3*N+JE ) = ZERO;
                  WORK( 2*N+JE+1 ) = -TEMP2R / TEMP;
                  WORK( 3*N+JE+1 ) = -TEMP2I / TEMP;
               } else {
                  WORK( 2*N+JE+1 ) = ONE;
                  WORK( 3*N+JE+1 ) = ZERO;
                  TEMP = ACOEF*S( JE, JE+1 );
                  WORK( 2*N+JE ) = ( BCOEFR*P( JE+1, JE+1 )-ACOEF* S( JE+1, JE+1 ) ) / TEMP;
                  WORK( 3*N+JE ) = BCOEFI*P( JE+1, JE+1 ) / TEMP;
               }
               XMAX = MAX( ABS( WORK( 2*N+JE ) )+ABS( WORK( 3*N+JE ) ), ABS( WORK( 2*N+JE+1 ) )+ABS( WORK( 3*N+JE+1 ) ) );
            }

            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN );

                                            // T
            // Triangular solve of  (a A - b B)  y = 0

                                    // T
            // (rowwise in  (a A - b B) , or columnwise in (a A - b B) )

            IL2BY2 = false;

            for (J = JE + NW; J <= N; J++) { // 160
               if ( IL2BY2 ) {
                  IL2BY2 = false;
                  GO TO 160;
               }

               NA = 1;
               BDIAG( 1 ) = P( J, J );
               if ( J < N ) {
                  if ( S( J+1, J ) != ZERO ) {
                     IL2BY2 = true;
                     BDIAG( 2 ) = P( J+1, J+1 );
                     NA = 2;
                  }
               }

               // Check whether scaling is necessary for dot products

               XSCALE = ONE / MAX( ONE, XMAX );
               TEMP = MAX( WORK( J ), WORK( N+J ), ACOEFA*WORK( J )+BCOEFA*WORK( N+J ) )                IF( IL2BY2 ) TEMP = MAX( TEMP, WORK( J+1 ), WORK( N+J+1 ), ACOEFA*WORK( J+1 )+BCOEFA*WORK( N+J+1 ) );
               if ( TEMP > BIGNUM*XSCALE ) {
                  for (JW = 0; JW <= NW - 1; JW++) { // 90
                     for (JR = JE; JR <= J - 1; JR++) { // 80
                        WORK( ( JW+2 )*N+JR ) = XSCALE* WORK( ( JW+2 )*N+JR );
                     } // 80
                  } // 90
                  XMAX = XMAX*XSCALE;
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
                     SUMS( JA, JW ) = ZERO;
                     SUMP( JA, JW ) = ZERO;

                     for (JR = JE; JR <= J - 1; JR++) { // 100
                        SUMS( JA, JW ) = SUMS( JA, JW ) + S( JR, J+JA-1 )* WORK( ( JW+1 )*N+JR )                         SUMP( JA, JW ) = SUMP( JA, JW ) + P( JR, J+JA-1 )* WORK( ( JW+1 )*N+JR );
                     } // 100
                  } // 110
               } // 120

               for (JA = 1; JA <= NA; JA++) { // 130
                  if ( ILCPLX ) {
                     SUM( JA, 1 ) = -ACOEF*SUMS( JA, 1 ) + BCOEFR*SUMP( JA, 1 ) - BCOEFI*SUMP( JA, 2 )                      SUM( JA, 2 ) = -ACOEF*SUMS( JA, 2 ) + BCOEFR*SUMP( JA, 2 ) + BCOEFI*SUMP( JA, 1 );
                  } else {
                     SUM( JA, 1 ) = -ACOEF*SUMS( JA, 1 ) + BCOEFR*SUMP( JA, 1 );
                  }
               } // 130

                                   // T
               // Solve  ( a A - b B )  y = SUM(,)
               // with scaling and perturbation of the denominator

               dlaln2( true , NA, NW, DMIN, ACOEF, S( J, J ), LDS, BDIAG( 1 ), BDIAG( 2 ), SUM, 2, BCOEFR, BCOEFI, WORK( 2*N+J ), N, SCALE, TEMP, IINFO );
               if ( SCALE < ONE ) {
                  for (JW = 0; JW <= NW - 1; JW++) { // 150
                     for (JR = JE; JR <= J - 1; JR++) { // 140
                        WORK( ( JW+2 )*N+JR ) = SCALE* WORK( ( JW+2 )*N+JR );
                     } // 140
                  } // 150
                  XMAX = SCALE*XMAX;
               }
               XMAX = MAX( XMAX, TEMP );
            } // 160

            // Copy eigenvector to VL, back transforming if
            // HOWMNY='B'.

            IEIG = IEIG + 1;
            if ( ILBACK ) {
               for (JW = 0; JW <= NW - 1; JW++) { // 170
                  dgemv('N', N, N+1-JE, ONE, VL( 1, JE ), LDVL, WORK( ( JW+2 )*N+JE ), 1, ZERO, WORK( ( JW+4 )*N+1 ), 1 );
               } // 170
               dlacpy(' ', N, NW, WORK( 4*N+1 ), N, VL( 1, JE ), LDVL );
               IBEG = 1;
            } else {
               dlacpy(' ', N, NW, WORK( 2*N+1 ), N, VL( 1, IEIG ), LDVL );
               IBEG = JE;
            }

            // Scale eigenvector

            XMAX = ZERO;
            if ( ILCPLX ) {
               for (J = IBEG; J <= N; J++) { // 180
                  XMAX = MAX( XMAX, ABS( VL( J, IEIG ) )+ ABS( VL( J, IEIG+1 ) ) );
               } // 180
            } else {
               for (J = IBEG; J <= N; J++) { // 190
                  XMAX = MAX( XMAX, ABS( VL( J, IEIG ) ) );
               } // 190
            }

            if ( XMAX > SAFMIN ) {
               XSCALE = ONE / XMAX;

               for (JW = 0; JW <= NW - 1; JW++) { // 210
                  for (JR = IBEG; JR <= N; JR++) { // 200
                     VL( JR, IEIG+JW ) = XSCALE*VL( JR, IEIG+JW );
                  } // 200
               } // 210
            }
            IEIG = IEIG + NW - 1;

         } // 220
      }

      // Right eigenvectors

      if ( COMPR ) {
         IEIG = IM + 1;

         // Main loop over eigenvalues

         ILCPLX = false;
         DO 500 JE = N, 1, -1;

            // Skip this iteration if (a) HOWMNY='S' and SELECT= false , or
            // (b) this would be the second of a complex pair.
            // Check for complex eigenvalue, so as to be sure of which
            // entry(-ies) of SELECT to look at -- if complex, SELECT(JE)
            // or SELECT(JE-1).
            // If this is a complex pair, the 2-by-2 diagonal block
            // corresponding to the eigenvalue is in rows/columns JE-1:JE

            if ( ILCPLX ) {
               ILCPLX = false;
               GO TO 500;
            }
            NW = 1;
            if ( JE > 1 ) {
               if ( S( JE, JE-1 ) != ZERO ) {
                  ILCPLX = true;
                  NW = 2;
               }
            }
            if ( ILALL ) {
               ILCOMP = true;
            } else if ( ILCPLX ) {
               ILCOMP = SELECT( JE ) || SELECT( JE-1 );
            } else {
               ILCOMP = SELECT( JE );
            }
            if ( !ILCOMP) GO TO 500;

            // Decide if (a) singular pencil, (b) real eigenvalue, or
            // (c) complex eigenvalue.

            if ( !ILCPLX ) {
               if ( ABS( S( JE, JE ) ) <= SAFMIN && ABS( P( JE, JE ) ) <= SAFMIN ) {

                  // Singular matrix pencil -- unit eigenvector

                  IEIG = IEIG - 1;
                  for (JR = 1; JR <= N; JR++) { // 230
                     VR( JR, IEIG ) = ZERO;
                  } // 230
                  VR( IEIG, IEIG ) = ONE;
                  GO TO 500;
               }
            }

            // Clear vector

            for (JW = 0; JW <= NW - 1; JW++) { // 250
               for (JR = 1; JR <= N; JR++) { // 240
                  WORK( ( JW+2 )*N+JR ) = ZERO;
               } // 240
            } // 250

            // Compute coefficients in  ( a A - b B ) x = 0
               // a  is  ACOEF
               // b  is  BCOEFR + i*BCOEFI

            if ( !ILCPLX ) {

               // Real eigenvalue

               TEMP = ONE / MAX( ABS( S( JE, JE ) )*ASCALE, ABS( P( JE, JE ) )*BSCALE, SAFMIN );
               SALFAR = ( TEMP*S( JE, JE ) )*ASCALE;
               SBETA = ( TEMP*P( JE, JE ) )*BSCALE;
               ACOEF = SBETA*ASCALE;
               BCOEFR = SALFAR*BSCALE;
               BCOEFI = ZERO;

               // Scale to avoid underflow

               SCALE = ONE;
               LSA = ABS( SBETA ) >= SAFMIN && ABS( ACOEF ) < SMALL;
               LSB = ABS( SALFAR ) >= SAFMIN && ABS( BCOEFR ) < SMALL                IF( LSA ) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )                IF( LSB ) SCALE = MAX( SCALE, ( SMALL / ABS( SALFAR ) )* MIN( BNORM, BIG ) );
               if ( LSA || LSB ) {
                  SCALE = MIN( SCALE, ONE / ( SAFMIN*MAX( ONE, ABS( ACOEF ), ABS( BCOEFR ) ) ) );
                  if ( LSA ) {
                     ACOEF = ASCALE*( SCALE*SBETA );
                  } else {
                     ACOEF = SCALE*ACOEF;
                  }
                  if ( LSB ) {
                     BCOEFR = BSCALE*( SCALE*SALFAR );
                  } else {
                     BCOEFR = SCALE*BCOEFR;
                  }
               }
               ACOEFA = ABS( ACOEF );
               BCOEFA = ABS( BCOEFR );

               // First component is 1

               WORK( 2*N+JE ) = ONE;
               XMAX = ONE;

               // Compute contribution from column JE of A and B to sum
               // (See "Further Details", above.)

               for (JR = 1; JR <= JE - 1; JR++) { // 260
                  WORK( 2*N+JR ) = BCOEFR*P( JR, JE ) - ACOEF*S( JR, JE );
               } // 260
            } else {

               // Complex eigenvalue

               dlag2(S( JE-1, JE-1 ), LDS, P( JE-1, JE-1 ), LDP, SAFMIN*SAFETY, ACOEF, TEMP, BCOEFR, TEMP2, BCOEFI );
               if ( BCOEFI == ZERO ) {
                  INFO = JE - 1;
                  return;
               }

               // Scale to avoid over/underflow

               ACOEFA = ABS( ACOEF );
               BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI );
               SCALE = ONE;
               if (ACOEFA*ULP < SAFMIN && ACOEFA >= SAFMIN) SCALE = ( SAFMIN / ULP ) / ACOEFA;
               if( BCOEFA*ULP < SAFMIN && BCOEFA >= SAFMIN ) SCALE = MAX( SCALE, ( SAFMIN / ULP ) / BCOEFA );
               if( SAFMIN*ACOEFA > ASCALE ) SCALE = ASCALE / ( SAFMIN*ACOEFA );
               IF( SAFMIN*BCOEFA > BSCALE ) SCALE = MIN( SCALE, BSCALE / ( SAFMIN*BCOEFA ) );
               if ( SCALE != ONE ) {
                  ACOEF = SCALE*ACOEF;
                  ACOEFA = ABS( ACOEF );
                  BCOEFR = SCALE*BCOEFR;
                  BCOEFI = SCALE*BCOEFI;
                  BCOEFA = ABS( BCOEFR ) + ABS( BCOEFI );
               }

               // Compute first two components of eigenvector
               // and contribution to sums

               TEMP = ACOEF*S( JE, JE-1 );
               TEMP2R = ACOEF*S( JE, JE ) - BCOEFR*P( JE, JE );
               TEMP2I = -BCOEFI*P( JE, JE );
               if ( ABS( TEMP ) >= ABS( TEMP2R )+ABS( TEMP2I ) ) {
                  WORK( 2*N+JE ) = ONE;
                  WORK( 3*N+JE ) = ZERO;
                  WORK( 2*N+JE-1 ) = -TEMP2R / TEMP;
                  WORK( 3*N+JE-1 ) = -TEMP2I / TEMP;
               } else {
                  WORK( 2*N+JE-1 ) = ONE;
                  WORK( 3*N+JE-1 ) = ZERO;
                  TEMP = ACOEF*S( JE-1, JE );
                  WORK( 2*N+JE ) = ( BCOEFR*P( JE-1, JE-1 )-ACOEF* S( JE-1, JE-1 ) ) / TEMP;
                  WORK( 3*N+JE ) = BCOEFI*P( JE-1, JE-1 ) / TEMP;
               }

               XMAX = MAX( ABS( WORK( 2*N+JE ) )+ABS( WORK( 3*N+JE ) ), ABS( WORK( 2*N+JE-1 ) )+ABS( WORK( 3*N+JE-1 ) ) );

               // Compute contribution from columns JE and JE-1
               // of A and B to the sums.

               CREALA = ACOEF*WORK( 2*N+JE-1 );
               CIMAGA = ACOEF*WORK( 3*N+JE-1 );
               CREALB = BCOEFR*WORK( 2*N+JE-1 ) - BCOEFI*WORK( 3*N+JE-1 )                CIMAGB = BCOEFI*WORK( 2*N+JE-1 ) + BCOEFR*WORK( 3*N+JE-1 );
               CRE2A = ACOEF*WORK( 2*N+JE );
               CIM2A = ACOEF*WORK( 3*N+JE );
               CRE2B = BCOEFR*WORK( 2*N+JE ) - BCOEFI*WORK( 3*N+JE );
               CIM2B = BCOEFI*WORK( 2*N+JE ) + BCOEFR*WORK( 3*N+JE );
               for (JR = 1; JR <= JE - 2; JR++) { // 270
                  WORK( 2*N+JR ) = -CREALA*S( JR, JE-1 ) + CREALB*P( JR, JE-1 ) - CRE2A*S( JR, JE ) + CRE2B*P( JR, JE )                   WORK( 3*N+JR ) = -CIMAGA*S( JR, JE-1 ) + CIMAGB*P( JR, JE-1 ) - CIM2A*S( JR, JE ) + CIM2B*P( JR, JE );
               } // 270
            }

            DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN );

            // Columnwise triangular solve of  (a A - b B)  x = 0

            IL2BY2 = false;
            DO 370 J = JE - NW, 1, -1;

               // If a 2-by-2 block, is in position j-1:j, wait until
               // next iteration to process it (when it will be j:j+1)

               if ( !IL2BY2 && J > 1 ) {
                  if ( S( J, J-1 ) != ZERO ) {
                     IL2BY2 = true;
                     GO TO 370;
                  }
               }
               BDIAG( 1 ) = P( J, J );
               if ( IL2BY2 ) {
                  NA = 2;
                  BDIAG( 2 ) = P( J+1, J+1 );
               } else {
                  NA = 1;
               }

               // Compute x(j) (and x(j+1), if 2-by-2 block)

               dlaln2( false , NA, NW, DMIN, ACOEF, S( J, J ), LDS, BDIAG( 1 ), BDIAG( 2 ), WORK( 2*N+J ), N, BCOEFR, BCOEFI, SUM, 2, SCALE, TEMP, IINFO );
               if ( SCALE < ONE ) {

                  for (JW = 0; JW <= NW - 1; JW++) { // 290
                     for (JR = 1; JR <= JE; JR++) { // 280
                        WORK( ( JW+2 )*N+JR ) = SCALE* WORK( ( JW+2 )*N+JR );
                     } // 280
                  } // 290
               }
               XMAX = MAX( SCALE*XMAX, TEMP );

               for (JW = 1; JW <= NW; JW++) { // 310
                  for (JA = 1; JA <= NA; JA++) { // 300
                     WORK( ( JW+1 )*N+J+JA-1 ) = SUM( JA, JW );
                  } // 300
               } // 310

               // w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling

               if ( J > 1 ) {

                  // Check whether scaling is necessary for sum.

                  XSCALE = ONE / MAX( ONE, XMAX );
                  TEMP = ACOEFA*WORK( J ) + BCOEFA*WORK( N+J );
                  if (IL2BY2) TEMP = MAX( TEMP, ACOEFA*WORK( J+1 )+BCOEFA* WORK( N+J+1 ) );
                  TEMP = MAX( TEMP, ACOEFA, BCOEFA );
                  if ( TEMP > BIGNUM*XSCALE ) {

                     for (JW = 0; JW <= NW - 1; JW++) { // 330
                        for (JR = 1; JR <= JE; JR++) { // 320
                           WORK( ( JW+2 )*N+JR ) = XSCALE* WORK( ( JW+2 )*N+JR );
                        } // 320
                     } // 330
                     XMAX = XMAX*XSCALE;
                  }

                  // Compute the contributions of the off-diagonals of
                  // column j (and j+1, if 2-by-2 block) of A and B to the
                  // sums.


                  for (JA = 1; JA <= NA; JA++) { // 360
                     if ( ILCPLX ) {
                        CREALA = ACOEF*WORK( 2*N+J+JA-1 );
                        CIMAGA = ACOEF*WORK( 3*N+J+JA-1 );
                        CREALB = BCOEFR*WORK( 2*N+J+JA-1 ) - BCOEFI*WORK( 3*N+J+JA-1 )                         CIMAGB = BCOEFI*WORK( 2*N+J+JA-1 ) + BCOEFR*WORK( 3*N+J+JA-1 );
                        for (JR = 1; JR <= J - 1; JR++) { // 340
                           WORK( 2*N+JR ) = WORK( 2*N+JR ) - CREALA*S( JR, J+JA-1 ) + CREALB*P( JR, J+JA-1 )                            WORK( 3*N+JR ) = WORK( 3*N+JR ) - CIMAGA*S( JR, J+JA-1 ) + CIMAGB*P( JR, J+JA-1 );
                        } // 340
                     } else {
                        CREALA = ACOEF*WORK( 2*N+J+JA-1 );
                        CREALB = BCOEFR*WORK( 2*N+J+JA-1 );
                        for (JR = 1; JR <= J - 1; JR++) { // 350
                           WORK( 2*N+JR ) = WORK( 2*N+JR ) - CREALA*S( JR, J+JA-1 ) + CREALB*P( JR, J+JA-1 );
                        } // 350
                     }
                  } // 360
               }

               IL2BY2 = false;
            } // 370

            // Copy eigenvector to VR, back transforming if
            // HOWMNY='B'.

            IEIG = IEIG - NW;
            if ( ILBACK ) {

               for (JW = 0; JW <= NW - 1; JW++) { // 410
                  for (JR = 1; JR <= N; JR++) { // 380
                     WORK( ( JW+4 )*N+JR ) = WORK( ( JW+2 )*N+1 )* VR( JR, 1 );
                  } // 380

                  // A series of compiler directives to defeat
                  // vectorization for the next loop


                  for (JC = 2; JC <= JE; JC++) { // 400
                     for (JR = 1; JR <= N; JR++) { // 390
                        WORK( ( JW+4 )*N+JR ) = WORK( ( JW+4 )*N+JR ) + WORK( ( JW+2 )*N+JC )*VR( JR, JC );
                     } // 390
                  } // 400
               } // 410

               for (JW = 0; JW <= NW - 1; JW++) { // 430
                  for (JR = 1; JR <= N; JR++) { // 420
                     VR( JR, IEIG+JW ) = WORK( ( JW+4 )*N+JR );
                  } // 420
               } // 430

               IEND = N;
            } else {
               for (JW = 0; JW <= NW - 1; JW++) { // 450
                  for (JR = 1; JR <= N; JR++) { // 440
                     VR( JR, IEIG+JW ) = WORK( ( JW+2 )*N+JR );
                  } // 440
               } // 450

               IEND = JE;
            }

            // Scale eigenvector

            XMAX = ZERO;
            if ( ILCPLX ) {
               for (J = 1; J <= IEND; J++) { // 460
                  XMAX = MAX( XMAX, ABS( VR( J, IEIG ) )+ ABS( VR( J, IEIG+1 ) ) );
               } // 460
            } else {
               for (J = 1; J <= IEND; J++) { // 470
                  XMAX = MAX( XMAX, ABS( VR( J, IEIG ) ) );
               } // 470
            }

            if ( XMAX > SAFMIN ) {
               XSCALE = ONE / XMAX;
               for (JW = 0; JW <= NW - 1; JW++) { // 490
                  for (JR = 1; JR <= IEND; JR++) { // 480
                     VR( JR, IEIG+JW ) = XSCALE*VR( JR, IEIG+JW );
                  } // 480
               } // 490
            }
         } // 500
      }

      return;

      // End of DTGEVC

      }
