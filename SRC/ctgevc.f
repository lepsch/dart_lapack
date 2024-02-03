      SUBROUTINE CTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL, LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, SIDE;
      int                INFO, LDP, LDS, LDVL, LDVR, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      REAL               RWORK( * )
      COMPLEX            P( LDP, * ), S( LDS, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..


*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               COMPL, COMPR, ILALL, ILBACK, ILBBAD, ILCOMP, LSA, LSB;
      int                I, IBEG, IEIG, IEND, IHWMNY, IM, ISIDE, ISRC, J, JE, JR;
      REAL               ACOEFA, ACOEFF, ANORM, ASCALE, BCOEFA, BIG, BIGNUM, BNORM, BSCALE, DMIN, SAFMIN, SBETA, SCALE, SMALL, TEMP, ULP, XMAX;
      COMPLEX            BCOEFF, CA, CB, D, SALPHA, SUM, SUMA, SUMB, X
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH
      COMPLEX            CLADIV
      // EXTERNAL LSAME, SLAMCH, CLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, CONJG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               ABS1
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( REAL( X ) ) + ABS( AIMAG( X ) )
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
         CALL XERBLA( 'CTGEVC', -INFO )
         RETURN
      }

      // Count the number of eigenvectors

      if ( .NOT.ILALL ) {
         IM = 0
         DO 10 J = 1, N
            IF( SELECT( J ) ) IM = IM + 1
   10    CONTINUE
      } else {
         IM = N
      }

      // Check diagonal of B

      ILBBAD = .FALSE.
      DO 20 J = 1, N
         IF( AIMAG( P( J, J ) ).NE.ZERO ) ILBBAD = .TRUE.
   20 CONTINUE

      if ( ILBBAD ) {
         INFO = -7
      } else if ( COMPL .AND. LDVL.LT.N .OR. LDVL.LT.1 ) {
         INFO = -10
      } else if ( COMPR .AND. LDVR.LT.N .OR. LDVR.LT.1 ) {
         INFO = -12
      } else if ( MM.LT.IM ) {
         INFO = -13
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CTGEVC', -INFO )
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
      // part of A and B to check for possible overflow in the triangular
      // solver.

      ANORM = ABS1( S( 1, 1 ) )
      BNORM = ABS1( P( 1, 1 ) )
      RWORK( 1 ) = ZERO
      RWORK( N+1 ) = ZERO
      DO 40 J = 2, N
         RWORK( J ) = ZERO
         RWORK( N+J ) = ZERO
         DO 30 I = 1, J - 1
            RWORK( J ) = RWORK( J ) + ABS1( S( I, J ) )
            RWORK( N+J ) = RWORK( N+J ) + ABS1( P( I, J ) )
   30    CONTINUE
         ANORM = MAX( ANORM, RWORK( J )+ABS1( S( J, J ) ) )
         BNORM = MAX( BNORM, RWORK( N+J )+ABS1( P( J, J ) ) )
   40 CONTINUE

      ASCALE = ONE / MAX( ANORM, SAFMIN )
      BSCALE = ONE / MAX( BNORM, SAFMIN )

      // Left eigenvectors

      if ( COMPL ) {
         IEIG = 0

         // Main loop over eigenvalues

         DO 140 JE = 1, N
            if ( ILALL ) {
               ILCOMP = .TRUE.
            } else {
               ILCOMP = SELECT( JE )
            }
            if ( ILCOMP ) {
               IEIG = IEIG + 1

               if ( ABS1( S( JE, JE ) ).LE.SAFMIN .AND. ABS( REAL( P( JE, JE ) ) ).LE.SAFMIN ) {

                  // Singular matrix pencil -- return unit eigenvector

                  DO 50 JR = 1, N
                     VL( JR, IEIG ) = CZERO
   50             CONTINUE
                  VL( IEIG, IEIG ) = CONE
                  GO TO 140
               }

               // Non-singular eigenvalue:
               // Compute coefficients  a  and  b  in
                    // H
                  // y  ( a A - b B ) = 0

               TEMP = ONE / MAX( ABS1( S( JE, JE ) )*ASCALE, ABS( REAL( P( JE, JE ) ) )*BSCALE, SAFMIN )
               SALPHA = ( TEMP*S( JE, JE ) )*ASCALE
               SBETA = ( TEMP*REAL( P( JE, JE ) ) )*BSCALE
               ACOEFF = SBETA*ASCALE
               BCOEFF = SALPHA*BSCALE

               // Scale to avoid underflow

               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEFF ).LT.SMALL
               LSB = ABS1( SALPHA ).GE.SAFMIN .AND. ABS1( BCOEFF ).LT. SMALL

               SCALE = ONE
               IF( LSA ) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )                IF( LSB ) SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )* MIN( BNORM, BIG ) )
               if ( LSA .OR. LSB ) {
                  SCALE = MIN( SCALE, ONE / ( SAFMIN*MAX( ONE, ABS( ACOEFF ), ABS1( BCOEFF ) ) ) )
                  if ( LSA ) {
                     ACOEFF = ASCALE*( SCALE*SBETA )
                  } else {
                     ACOEFF = SCALE*ACOEFF
                  }
                  if ( LSB ) {
                     BCOEFF = BSCALE*( SCALE*SALPHA )
                  } else {
                     BCOEFF = SCALE*BCOEFF
                  }
               }

               ACOEFA = ABS( ACOEFF )
               BCOEFA = ABS1( BCOEFF )
               XMAX = ONE
               DO 60 JR = 1, N
                  WORK( JR ) = CZERO
   60          CONTINUE
               WORK( JE ) = CONE
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )

                                               // H
               // Triangular solve of  (a A - b B)  y = 0

                                       // H
               // (rowwise in  (a A - b B) , or columnwise in a A - b B)

               DO 100 J = JE + 1, N

                  // Compute
                        // j-1
                  // SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
                        // k=je
                  // (Scale if necessary)

                  TEMP = ONE / XMAX
                  if ( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ).GT.BIGNUM* TEMP ) {
                     DO 70 JR = JE, J - 1
                        WORK( JR ) = TEMP*WORK( JR )
   70                CONTINUE
                     XMAX = ONE
                  }
                  SUMA = CZERO
                  SUMB = CZERO

                  DO 80 JR = JE, J - 1
                     SUMA = SUMA + CONJG( S( JR, J ) )*WORK( JR )
                     SUMB = SUMB + CONJG( P( JR, J ) )*WORK( JR )
   80             CONTINUE
                  SUM = ACOEFF*SUMA - CONJG( BCOEFF )*SUMB

                  // Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) )

                  // with scaling and perturbation of the denominator

                  D = CONJG( ACOEFF*S( J, J )-BCOEFF*P( J, J ) )
                  IF( ABS1( D ).LE.DMIN ) D = CMPLX( DMIN )

                  if ( ABS1( D ).LT.ONE ) {
                     if ( ABS1( SUM ).GE.BIGNUM*ABS1( D ) ) {
                        TEMP = ONE / ABS1( SUM )
                        DO 90 JR = JE, J - 1
                           WORK( JR ) = TEMP*WORK( JR )
   90                   CONTINUE
                        XMAX = TEMP*XMAX
                        SUM = TEMP*SUM
                     }
                  }
                  WORK( J ) = CLADIV( -SUM, D )
                  XMAX = MAX( XMAX, ABS1( WORK( J ) ) )
  100          CONTINUE

               // Back transform eigenvector if HOWMNY='B'.

               if ( ILBACK ) {
                  CALL CGEMV( 'N', N, N+1-JE, CONE, VL( 1, JE ), LDVL, WORK( JE ), 1, CZERO, WORK( N+1 ), 1 )
                  ISRC = 2
                  IBEG = 1
               } else {
                  ISRC = 1
                  IBEG = JE
               }

               // Copy and scale eigenvector into column of VL

               XMAX = ZERO
               DO 110 JR = IBEG, N
                  XMAX = MAX( XMAX, ABS1( WORK( ( ISRC-1 )*N+JR ) ) )
  110          CONTINUE

               if ( XMAX.GT.SAFMIN ) {
                  TEMP = ONE / XMAX
                  DO 120 JR = IBEG, N
                     VL( JR, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+JR )
  120             CONTINUE
               } else {
                  IBEG = N + 1
               }

               DO 130 JR = 1, IBEG - 1
                  VL( JR, IEIG ) = CZERO
  130          CONTINUE

            }
  140    CONTINUE
      }

      // Right eigenvectors

      if ( COMPR ) {
         IEIG = IM + 1

         // Main loop over eigenvalues

         DO 250 JE = N, 1, -1
            if ( ILALL ) {
               ILCOMP = .TRUE.
            } else {
               ILCOMP = SELECT( JE )
            }
            if ( ILCOMP ) {
               IEIG = IEIG - 1

               if ( ABS1( S( JE, JE ) ).LE.SAFMIN .AND. ABS( REAL( P( JE, JE ) ) ).LE.SAFMIN ) {

                  // Singular matrix pencil -- return unit eigenvector

                  DO 150 JR = 1, N
                     VR( JR, IEIG ) = CZERO
  150             CONTINUE
                  VR( IEIG, IEIG ) = CONE
                  GO TO 250
               }

               // Non-singular eigenvalue:
               // Compute coefficients  a  and  b  in

               // ( a A - b B ) x  = 0

               TEMP = ONE / MAX( ABS1( S( JE, JE ) )*ASCALE, ABS( REAL( P( JE, JE ) ) )*BSCALE, SAFMIN )
               SALPHA = ( TEMP*S( JE, JE ) )*ASCALE
               SBETA = ( TEMP*REAL( P( JE, JE ) ) )*BSCALE
               ACOEFF = SBETA*ASCALE
               BCOEFF = SALPHA*BSCALE

               // Scale to avoid underflow

               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEFF ).LT.SMALL
               LSB = ABS1( SALPHA ).GE.SAFMIN .AND. ABS1( BCOEFF ).LT. SMALL

               SCALE = ONE
               IF( LSA ) SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )                IF( LSB ) SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )* MIN( BNORM, BIG ) )
               if ( LSA .OR. LSB ) {
                  SCALE = MIN( SCALE, ONE / ( SAFMIN*MAX( ONE, ABS( ACOEFF ), ABS1( BCOEFF ) ) ) )
                  if ( LSA ) {
                     ACOEFF = ASCALE*( SCALE*SBETA )
                  } else {
                     ACOEFF = SCALE*ACOEFF
                  }
                  if ( LSB ) {
                     BCOEFF = BSCALE*( SCALE*SALPHA )
                  } else {
                     BCOEFF = SCALE*BCOEFF
                  }
               }

               ACOEFA = ABS( ACOEFF )
               BCOEFA = ABS1( BCOEFF )
               XMAX = ONE
               DO 160 JR = 1, N
                  WORK( JR ) = CZERO
  160          CONTINUE
               WORK( JE ) = CONE
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )

               // Triangular solve of  (a A - b B) x = 0  (columnwise)

               // WORK(1:j-1) contains sums w,
               // WORK(j+1:JE) contains x

               DO 170 JR = 1, JE - 1
                  WORK( JR ) = ACOEFF*S( JR, JE ) - BCOEFF*P( JR, JE )
  170          CONTINUE
               WORK( JE ) = CONE

               DO 210 J = JE - 1, 1, -1

                  // Form x(j) := - w(j) / d
                  // with scaling and perturbation of the denominator

                  D = ACOEFF*S( J, J ) - BCOEFF*P( J, J )
                  IF( ABS1( D ).LE.DMIN ) D = CMPLX( DMIN )

                  if ( ABS1( D ).LT.ONE ) {
                     if ( ABS1( WORK( J ) ).GE.BIGNUM*ABS1( D ) ) {
                        TEMP = ONE / ABS1( WORK( J ) )
                        DO 180 JR = 1, JE
                           WORK( JR ) = TEMP*WORK( JR )
  180                   CONTINUE
                     }
                  }

                  WORK( J ) = CLADIV( -WORK( J ), D )

                  if ( J.GT.1 ) {

                     // w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling

                     if ( ABS1( WORK( J ) ).GT.ONE ) {
                        TEMP = ONE / ABS1( WORK( J ) )
                        if ( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ).GE. BIGNUM*TEMP ) {
                           DO 190 JR = 1, JE
                              WORK( JR ) = TEMP*WORK( JR )
  190                      CONTINUE
                        }
                     }

                     CA = ACOEFF*WORK( J )
                     CB = BCOEFF*WORK( J )
                     DO 200 JR = 1, J - 1
                        WORK( JR ) = WORK( JR ) + CA*S( JR, J ) - CB*P( JR, J )
  200                CONTINUE
                  }
  210          CONTINUE

               // Back transform eigenvector if HOWMNY='B'.

               if ( ILBACK ) {
                  CALL CGEMV( 'N', N, JE, CONE, VR, LDVR, WORK, 1, CZERO, WORK( N+1 ), 1 )
                  ISRC = 2
                  IEND = N
               } else {
                  ISRC = 1
                  IEND = JE
               }

               // Copy and scale eigenvector into column of VR

               XMAX = ZERO
               DO 220 JR = 1, IEND
                  XMAX = MAX( XMAX, ABS1( WORK( ( ISRC-1 )*N+JR ) ) )
  220          CONTINUE

               if ( XMAX.GT.SAFMIN ) {
                  TEMP = ONE / XMAX
                  DO 230 JR = 1, IEND
                     VR( JR, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+JR )
  230             CONTINUE
               } else {
                  IEND = 0
               }

               DO 240 JR = IEND + 1, N
                  VR( JR, IEIG ) = CZERO
  240          CONTINUE

            }
  250    CONTINUE
      }

      RETURN

      // End of CTGEVC

      }
