      SUBROUTINE ZTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDP, LDS, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      double             RWORK( * );
      COMPLEX*16         P( LDP, * ), S( LDS, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..


*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               COMPL, COMPR, ILALL, ILBACK, ILBBAD, ILCOMP, LSA, LSB;
      int                I, IBEG, IEIG, IEND, IHWMNY, IM, ISIDE, ISRC, J, JE, JR;
      double             ACOEFA, ACOEFF, ANORM, ASCALE, BCOEFA, BIG, BIGNUM, BNORM, BSCALE, DMIN, SAFMIN, SBETA, SCALE, SMALL, TEMP, ULP, XMAX;
      COMPLEX*16         BCOEFF, CA, CB, D, SALPHA, SUM, SUMA, SUMB, X;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      COMPLEX*16         ZLADIV;
      // EXTERNAL LSAME, DLAMCH, ZLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) );
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
         xerbla('ZTGEVC', -INFO );
         RETURN;
      }

      // Count the number of eigenvectors

      if ( !ILALL ) {
         IM = 0;
         for (J = 1; J <= N; J++) { // 10
            IF( SELECT( J ) ) IM = IM + 1;
         } // 10
      } else {
         IM = N;
      }

      // Check diagonal of B

      ILBBAD = false;
      for (J = 1; J <= N; J++) { // 20
         IF( DIMAG( P( J, J ) ) != ZERO ) ILBBAD = true;
      } // 20

      if ( ILBBAD ) {
         INFO = -7;
      } else if ( COMPL && LDVL < N || LDVL < 1 ) {
         INFO = -10;
      } else if ( COMPR && LDVR < N || LDVR < 1 ) {
         INFO = -12;
      } else if ( MM < IM ) {
         INFO = -13;
      }
      if ( INFO != 0 ) {
         xerbla('ZTGEVC', -INFO );
         RETURN;
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
      // part of A and B to check for possible overflow in the triangular
      // solver.

      ANORM = ABS1( S( 1, 1 ) );
      BNORM = ABS1( P( 1, 1 ) );
      RWORK( 1 ) = ZERO;
      RWORK( N+1 ) = ZERO;
      for (J = 2; J <= N; J++) { // 40
         RWORK( J ) = ZERO;
         RWORK( N+J ) = ZERO;
         for (I = 1; I <= J - 1; I++) { // 30
            RWORK( J ) = RWORK( J ) + ABS1( S( I, J ) );
            RWORK( N+J ) = RWORK( N+J ) + ABS1( P( I, J ) );
         } // 30
         ANORM = MAX( ANORM, RWORK( J )+ABS1( S( J, J ) ) );
         BNORM = MAX( BNORM, RWORK( N+J )+ABS1( P( J, J ) ) );
      } // 40

      ASCALE = ONE / MAX( ANORM, SAFMIN );
      BSCALE = ONE / MAX( BNORM, SAFMIN );

      // Left eigenvectors

      if ( COMPL ) {
         IEIG = 0;

         // Main loop over eigenvalues

         for (JE = 1; JE <= N; JE++) { // 140
            if ( ILALL ) {
               ILCOMP = true;
            } else {
               ILCOMP = SELECT( JE );
            }
            if ( ILCOMP ) {
               IEIG = IEIG + 1;

               if ( ABS1( S( JE, JE ) ) <= SAFMIN && ABS( DBLE( P( JE, JE ) ) ) <= SAFMIN ) {

                  // Singular matrix pencil -- return unit eigenvector

                  for (JR = 1; JR <= N; JR++) { // 50
                     VL( JR, IEIG ) = CZERO;
                  } // 50
                  VL( IEIG, IEIG ) = CONE;
                  GO TO 140;
               }

               // Non-singular eigenvalue:
               // Compute coefficients  a  and  b  in
                    // H
                  // y  ( a A - b B ) = 0

               TEMP = ONE / MAX( ABS1( S( JE, JE ) )*ASCALE, ABS( DBLE( P( JE, JE ) ) )*BSCALE, SAFMIN );
               SALPHA = ( TEMP*S( JE, JE ) )*ASCALE;
               SBETA = ( TEMP*DBLE( P( JE, JE ) ) )*BSCALE;
               ACOEFF = SBETA*ASCALE;
               BCOEFF = SALPHA*BSCALE;

               // Scale to avoid underflow

               LSA = ABS( SBETA ) >= SAFMIN && ABS( ACOEFF ) < SMALL;
               LSB = ABS1( SALPHA ) >= SAFMIN && ABS1( BCOEFF ) < SMALL;

               SCALE = ONE;
               if (LSA) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )                IF( LSB ) SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )* MIN( BNORM, BIG ) );
               if ( LSA || LSB ) {
                  SCALE = MIN( SCALE, ONE / ( SAFMIN*MAX( ONE, ABS( ACOEFF ), ABS1( BCOEFF ) ) ) );
                  if ( LSA ) {
                     ACOEFF = ASCALE*( SCALE*SBETA );
                  } else {
                     ACOEFF = SCALE*ACOEFF;
                  }
                  if ( LSB ) {
                     BCOEFF = BSCALE*( SCALE*SALPHA );
                  } else {
                     BCOEFF = SCALE*BCOEFF;
                  }
               }

               ACOEFA = ABS( ACOEFF );
               BCOEFA = ABS1( BCOEFF );
               XMAX = ONE;
               for (JR = 1; JR <= N; JR++) { // 60
                  WORK( JR ) = CZERO;
               } // 60
               WORK( JE ) = CONE;
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN );

                                               // H
               // Triangular solve of  (a A - b B)  y = 0

                                       // H
               // (rowwise in  (a A - b B) , or columnwise in a A - b B)

               for (J = JE + 1; J <= N; J++) { // 100

                  // Compute
                        // j-1
                  // SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
                        // k=je
                  // (Scale if necessary)

                  TEMP = ONE / XMAX;
                  if ( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ) > BIGNUM* TEMP ) {
                     for (JR = JE; JR <= J - 1; JR++) { // 70
                        WORK( JR ) = TEMP*WORK( JR );
                     } // 70
                     XMAX = ONE;
                  }
                  SUMA = CZERO;
                  SUMB = CZERO;

                  for (JR = JE; JR <= J - 1; JR++) { // 80
                     SUMA = SUMA + DCONJG( S( JR, J ) )*WORK( JR );
                     SUMB = SUMB + DCONJG( P( JR, J ) )*WORK( JR );
                  } // 80
                  SUM = ACOEFF*SUMA - DCONJG( BCOEFF )*SUMB;

                  // Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) )

                  // with scaling and perturbation of the denominator

                  D = DCONJG( ACOEFF*S( J, J )-BCOEFF*P( J, J ) );
                  IF( ABS1( D ) <= DMIN ) D = DCMPLX( DMIN );

                  if ( ABS1( D ) < ONE ) {
                     if ( ABS1( SUM ) >= BIGNUM*ABS1( D ) ) {
                        TEMP = ONE / ABS1( SUM );
                        for (JR = JE; JR <= J - 1; JR++) { // 90
                           WORK( JR ) = TEMP*WORK( JR );
                        } // 90
                        XMAX = TEMP*XMAX;
                        SUM = TEMP*SUM;
                     }
                  }
                  WORK( J ) = ZLADIV( -SUM, D );
                  XMAX = MAX( XMAX, ABS1( WORK( J ) ) );
               } // 100

               // Back transform eigenvector if HOWMNY='B'.

               if ( ILBACK ) {
                  zgemv('N', N, N+1-JE, CONE, VL( 1, JE ), LDVL, WORK( JE ), 1, CZERO, WORK( N+1 ), 1 );
                  ISRC = 2;
                  IBEG = 1;
               } else {
                  ISRC = 1;
                  IBEG = JE;
               }

               // Copy and scale eigenvector into column of VL

               XMAX = ZERO;
               for (JR = IBEG; JR <= N; JR++) { // 110
                  XMAX = MAX( XMAX, ABS1( WORK( ( ISRC-1 )*N+JR ) ) );
               } // 110

               if ( XMAX > SAFMIN ) {
                  TEMP = ONE / XMAX;
                  for (JR = IBEG; JR <= N; JR++) { // 120
                     VL( JR, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+JR );
                  } // 120
               } else {
                  IBEG = N + 1;
               }

               for (JR = 1; JR <= IBEG - 1; JR++) { // 130
                  VL( JR, IEIG ) = CZERO;
               } // 130

            }
         } // 140
      }

      // Right eigenvectors

      if ( COMPR ) {
         IEIG = IM + 1;

         // Main loop over eigenvalues

         DO 250 JE = N, 1, -1;
            if ( ILALL ) {
               ILCOMP = true;
            } else {
               ILCOMP = SELECT( JE );
            }
            if ( ILCOMP ) {
               IEIG = IEIG - 1;

               if ( ABS1( S( JE, JE ) ) <= SAFMIN && ABS( DBLE( P( JE, JE ) ) ) <= SAFMIN ) {

                  // Singular matrix pencil -- return unit eigenvector

                  for (JR = 1; JR <= N; JR++) { // 150
                     VR( JR, IEIG ) = CZERO;
                  } // 150
                  VR( IEIG, IEIG ) = CONE;
                  GO TO 250;
               }

               // Non-singular eigenvalue:
               // Compute coefficients  a  and  b  in

               // ( a A - b B ) x  = 0

               TEMP = ONE / MAX( ABS1( S( JE, JE ) )*ASCALE, ABS( DBLE( P( JE, JE ) ) )*BSCALE, SAFMIN );
               SALPHA = ( TEMP*S( JE, JE ) )*ASCALE;
               SBETA = ( TEMP*DBLE( P( JE, JE ) ) )*BSCALE;
               ACOEFF = SBETA*ASCALE;
               BCOEFF = SALPHA*BSCALE;

               // Scale to avoid underflow

               LSA = ABS( SBETA ) >= SAFMIN && ABS( ACOEFF ) < SMALL;
               LSB = ABS1( SALPHA ) >= SAFMIN && ABS1( BCOEFF ) < SMALL;

               SCALE = ONE;
               if (LSA) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )                IF( LSB ) SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )* MIN( BNORM, BIG ) );
               if ( LSA || LSB ) {
                  SCALE = MIN( SCALE, ONE / ( SAFMIN*MAX( ONE, ABS( ACOEFF ), ABS1( BCOEFF ) ) ) );
                  if ( LSA ) {
                     ACOEFF = ASCALE*( SCALE*SBETA );
                  } else {
                     ACOEFF = SCALE*ACOEFF;
                  }
                  if ( LSB ) {
                     BCOEFF = BSCALE*( SCALE*SALPHA );
                  } else {
                     BCOEFF = SCALE*BCOEFF;
                  }
               }

               ACOEFA = ABS( ACOEFF );
               BCOEFA = ABS1( BCOEFF );
               XMAX = ONE;
               for (JR = 1; JR <= N; JR++) { // 160
                  WORK( JR ) = CZERO;
               } // 160
               WORK( JE ) = CONE;
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN );

               // Triangular solve of  (a A - b B) x = 0  (columnwise)

               // WORK(1:j-1) contains sums w,
               // WORK(j+1:JE) contains x

               for (JR = 1; JR <= JE - 1; JR++) { // 170
                  WORK( JR ) = ACOEFF*S( JR, JE ) - BCOEFF*P( JR, JE );
               } // 170
               WORK( JE ) = CONE;

               DO 210 J = JE - 1, 1, -1;

                  // Form x(j) := - w(j) / d
                  // with scaling and perturbation of the denominator

                  D = ACOEFF*S( J, J ) - BCOEFF*P( J, J );
                  IF( ABS1( D ) <= DMIN ) D = DCMPLX( DMIN );

                  if ( ABS1( D ) < ONE ) {
                     if ( ABS1( WORK( J ) ) >= BIGNUM*ABS1( D ) ) {
                        TEMP = ONE / ABS1( WORK( J ) );
                        for (JR = 1; JR <= JE; JR++) { // 180
                           WORK( JR ) = TEMP*WORK( JR );
                        } // 180
                     }
                  }

                  WORK( J ) = ZLADIV( -WORK( J ), D );

                  if ( J > 1 ) {

                     // w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling

                     if ( ABS1( WORK( J ) ) > ONE ) {
                        TEMP = ONE / ABS1( WORK( J ) );
                        if ( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ) >= BIGNUM*TEMP ) {
                           for (JR = 1; JR <= JE; JR++) { // 190
                              WORK( JR ) = TEMP*WORK( JR );
                           } // 190
                        }
                     }

                     CA = ACOEFF*WORK( J );
                     CB = BCOEFF*WORK( J );
                     for (JR = 1; JR <= J - 1; JR++) { // 200
                        WORK( JR ) = WORK( JR ) + CA*S( JR, J ) - CB*P( JR, J );
                     } // 200
                  }
               } // 210

               // Back transform eigenvector if HOWMNY='B'.

               if ( ILBACK ) {
                  zgemv('N', N, JE, CONE, VR, LDVR, WORK, 1, CZERO, WORK( N+1 ), 1 );
                  ISRC = 2;
                  IEND = N;
               } else {
                  ISRC = 1;
                  IEND = JE;
               }

               // Copy and scale eigenvector into column of VR

               XMAX = ZERO;
               for (JR = 1; JR <= IEND; JR++) { // 220
                  XMAX = MAX( XMAX, ABS1( WORK( ( ISRC-1 )*N+JR ) ) );
               } // 220

               if ( XMAX > SAFMIN ) {
                  TEMP = ONE / XMAX;
                  for (JR = 1; JR <= IEND; JR++) { // 230
                     VR( JR, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+JR );
                  } // 230
               } else {
                  IEND = 0;
               }

               for (JR = IEND + 1; JR <= N; JR++) { // 240
                  VR( JR, IEIG ) = CZERO;
               } // 240

            }
         } // 250
      }

      RETURN;

      // End of ZTGEVC

      }
