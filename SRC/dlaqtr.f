      SUBROUTINE DLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LREAL, LTRAN;
      int                INFO, LDT, N;
      double             SCALE, W;
      // ..
      // .. Array Arguments ..
      double             B( * ), T( LDT, * ), WORK( * ), X( * );
      // ..

* =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IERR, J, J1, J2, JNEXT, K, N1, N2;
      double             BIGNUM, EPS, REC, SCALOC, SI, SMIN, SMINW, SMLNUM, SR, TJJ, TMP, XJ, XMAX, XNORM, Z;
      // ..
      // .. Local Arrays ..
      double             D( 2, 2 ), V( 2, 2 );
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DASUM, DDOT, DLAMCH, DLANGE;
      // EXTERNAL IDAMAX, DASUM, DDOT, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLADIV, DLALN2, DSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Do not test the input parameters for errors

      NOTRAN = .NOT.LTRAN
      INFO = 0

      // Quick return if possible

      if (N == 0) RETURN;

      // Set constants to control overflow

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM

      XNORM = DLANGE( 'M', N, N, T, LDT, D )
      if (.NOT.LREAL) XNORM = MAX( XNORM, ABS( W ), DLANGE( 'M', N, 1, B, N, D ) );
      SMIN = MAX( SMLNUM, EPS*XNORM )

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      WORK( 1 ) = ZERO
      for (J = 2; J <= N; J++) { // 10
         WORK( J ) = DASUM( J-1, T( 1, J ), 1 )
      } // 10

      if ( .NOT.LREAL ) {
         for (I = 2; I <= N; I++) { // 20
            WORK( I ) = WORK( I ) + ABS( B( I ) )
         } // 20
      }

      N2 = 2*N
      N1 = N
      if (.NOT.LREAL) N1 = N2;
      K = IDAMAX( N1, X, 1 )
      XMAX = ABS( X( K ) )
      SCALE = ONE

      if ( XMAX.GT.BIGNUM ) {
         SCALE = BIGNUM / XMAX
         dscal(N1, SCALE, X, 1 );
         XMAX = BIGNUM
      }

      if ( LREAL ) {

         if ( NOTRAN ) {

            // Solve T*p = scale*c

            JNEXT = N
            DO 30 J = N, 1, -1
               if (J.GT.JNEXT) GO TO 30;
               J1 = J
               J2 = J
               JNEXT = J - 1
               if ( J.GT.1 ) {
                  if ( T( J, J-1 ) != ZERO ) {
                     J1 = J - 1
                     JNEXT = J - 2
                  }
               }

               if ( J1 == J2 ) {

                  // Meet 1 by 1 diagonal block

                  // Scale to avoid overflow when computing
                      // x(j) = b(j)/T(j,j)

                  XJ = ABS( X( J1 ) )
                  TJJ = ABS( T( J1, J1 ) )
                  TMP = T( J1, J1 )
                  if ( TJJ.LT.SMIN ) {
                     TMP = SMIN
                     TJJ = SMIN
                     INFO = 1
                  }

                  if (XJ == ZERO) GO TO 30;

                  if ( TJJ.LT.ONE ) {
                     if ( XJ.GT.BIGNUM*TJJ ) {
                        REC = ONE / XJ
                        dscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }
                  X( J1 ) = X( J1 ) / TMP
                  XJ = ABS( X( J1 ) )

                  // Scale x if necessary to avoid overflow when adding a
                  // multiple of column j1 of T.

                  if ( XJ.GT.ONE ) {
                     REC = ONE / XJ
                     if ( WORK( J1 ).GT.( BIGNUM-XMAX )*REC ) {
                        dscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                     }
                  }
                  if ( J1.GT.1 ) {
                     daxpy(J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 );
                     K = IDAMAX( J1-1, X, 1 )
                     XMAX = ABS( X( K ) )
                  }

               } else {

                  // Meet 2 by 2 diagonal block

                  // Call 2 by 2 linear system solve, to take
                  // care of possible overflow by scaling factor.

                  D( 1, 1 ) = X( J1 )
                  D( 2, 1 ) = X( J2 )
                  dlaln2( false , 2, 1, SMIN, ONE, T( J1, J1 ), LDT, ONE, ONE, D, 2, ZERO, ZERO, V, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 2;

                  if ( SCALOC != ONE ) {
                     dscal(N, SCALOC, X, 1 );
                     SCALE = SCALE*SCALOC
                  }
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )

                  // Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
                  // to avoid overflow in updating right-hand side.

                  XJ = MAX( ABS( V( 1, 1 ) ), ABS( V( 2, 1 ) ) )
                  if ( XJ.GT.ONE ) {
                     REC = ONE / XJ
                     if ( MAX( WORK( J1 ), WORK( J2 ) ).GT. ( BIGNUM-XMAX )*REC ) {
                        dscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                     }
                  }

                  // Update right-hand side

                  if ( J1.GT.1 ) {
                     daxpy(J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 );
                     daxpy(J1-1, -X( J2 ), T( 1, J2 ), 1, X, 1 );
                     K = IDAMAX( J1-1, X, 1 )
                     XMAX = ABS( X( K ) )
                  }

               }

            } // 30

         } else {

            // Solve T**T*p = scale*c

            JNEXT = 1
            for (J = 1; J <= N; J++) { // 40
               if (J.LT.JNEXT) GO TO 40;
               J1 = J
               J2 = J
               JNEXT = J + 1
               if ( J.LT.N ) {
                  if ( T( J+1, J ) != ZERO ) {
                     J2 = J + 1
                     JNEXT = J + 2
                  }
               }

               if ( J1 == J2 ) {

                  // 1 by 1 diagonal block

                  // Scale if necessary to avoid overflow in forming the
                  // right-hand side element by inner product.

                  XJ = ABS( X( J1 ) )
                  if ( XMAX.GT.ONE ) {
                     REC = ONE / XMAX
                     if ( WORK( J1 ).GT.( BIGNUM-XJ )*REC ) {
                        dscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }

                  X( J1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, 1 )

                  XJ = ABS( X( J1 ) )
                  TJJ = ABS( T( J1, J1 ) )
                  TMP = T( J1, J1 )
                  if ( TJJ.LT.SMIN ) {
                     TMP = SMIN
                     TJJ = SMIN
                     INFO = 1
                  }

                  if ( TJJ.LT.ONE ) {
                     if ( XJ.GT.BIGNUM*TJJ ) {
                        REC = ONE / XJ
                        dscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }
                  X( J1 ) = X( J1 ) / TMP
                  XMAX = MAX( XMAX, ABS( X( J1 ) ) )

               } else {

                  // 2 by 2 diagonal block

                  // Scale if necessary to avoid overflow in forming the
                  // right-hand side elements by inner product.

                  XJ = MAX( ABS( X( J1 ) ), ABS( X( J2 ) ) )
                  if ( XMAX.GT.ONE ) {
                     REC = ONE / XMAX
                     if ( MAX( WORK( J2 ), WORK( J1 ) ).GT.( BIGNUM-XJ )* REC ) {
                        dscal(N, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }

                  D( 1, 1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, 1 )                   D( 2, 1 ) = X( J2 ) - DDOT( J1-1, T( 1, J2 ), 1, X, 1 )

                  dlaln2( true , 2, 1, SMIN, ONE, T( J1, J1 ), LDT, ONE, ONE, D, 2, ZERO, ZERO, V, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 2;

                  if ( SCALOC != ONE ) {
                     dscal(N, SCALOC, X, 1 );
                     SCALE = SCALE*SCALOC
                  }
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  XMAX = MAX( ABS( X( J1 ) ), ABS( X( J2 ) ), XMAX )

               }
            } // 40
         }

      } else {

         SMINW = MAX( EPS*ABS( W ), SMIN )
         if ( NOTRAN ) {

            // Solve (T + iB)*(p+iq) = c+id

            JNEXT = N
            DO 70 J = N, 1, -1
               if (J.GT.JNEXT) GO TO 70;
               J1 = J
               J2 = J
               JNEXT = J - 1
               if ( J.GT.1 ) {
                  if ( T( J, J-1 ) != ZERO ) {
                     J1 = J - 1
                     JNEXT = J - 2
                  }
               }

               if ( J1 == J2 ) {

                  // 1 by 1 diagonal block

                  // Scale if necessary to avoid overflow in division

                  Z = W
                  if (J1 == 1) Z = B( 1 );
                  XJ = ABS( X( J1 ) ) + ABS( X( N+J1 ) )
                  TJJ = ABS( T( J1, J1 ) ) + ABS( Z )
                  TMP = T( J1, J1 )
                  if ( TJJ.LT.SMINW ) {
                     TMP = SMINW
                     TJJ = SMINW
                     INFO = 1
                  }

                  if (XJ == ZERO) GO TO 70;

                  if ( TJJ.LT.ONE ) {
                     if ( XJ.GT.BIGNUM*TJJ ) {
                        REC = ONE / XJ
                        dscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }
                  dladiv(X( J1 ), X( N+J1 ), TMP, Z, SR, SI );
                  X( J1 ) = SR
                  X( N+J1 ) = SI
                  XJ = ABS( X( J1 ) ) + ABS( X( N+J1 ) )

                  // Scale x if necessary to avoid overflow when adding a
                  // multiple of column j1 of T.

                  if ( XJ.GT.ONE ) {
                     REC = ONE / XJ
                     if ( WORK( J1 ).GT.( BIGNUM-XMAX )*REC ) {
                        dscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC
                     }
                  }

                  if ( J1.GT.1 ) {
                     daxpy(J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 );
                     daxpy(J1-1, -X( N+J1 ), T( 1, J1 ), 1, X( N+1 ), 1 );

                     X( 1 ) = X( 1 ) + B( J1 )*X( N+J1 )
                     X( N+1 ) = X( N+1 ) - B( J1 )*X( J1 )

                     XMAX = ZERO
                     for (K = 1; K <= J1 - 1; K++) { // 50
                        XMAX = MAX( XMAX, ABS( X( K ) )+ ABS( X( K+N ) ) )
                     } // 50
                  }

               } else {

                  // Meet 2 by 2 diagonal block

                  D( 1, 1 ) = X( J1 )
                  D( 2, 1 ) = X( J2 )
                  D( 1, 2 ) = X( N+J1 )
                  D( 2, 2 ) = X( N+J2 )
                  dlaln2( false , 2, 2, SMINW, ONE, T( J1, J1 ), LDT, ONE, ONE, D, 2, ZERO, -W, V, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 2;

                  if ( SCALOC != ONE ) {
                     dscal(2*N, SCALOC, X, 1 );
                     SCALE = SCALOC*SCALE
                  }
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  X( N+J1 ) = V( 1, 2 )
                  X( N+J2 ) = V( 2, 2 )

                  // Scale X(J1), .... to avoid overflow in
                  // updating right hand side.

                  XJ = MAX( ABS( V( 1, 1 ) )+ABS( V( 1, 2 ) ), ABS( V( 2, 1 ) )+ABS( V( 2, 2 ) ) )
                  if ( XJ.GT.ONE ) {
                     REC = ONE / XJ
                     if ( MAX( WORK( J1 ), WORK( J2 ) ).GT. ( BIGNUM-XMAX )*REC ) {
                        dscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC
                     }
                  }

                  // Update the right-hand side.

                  if ( J1.GT.1 ) {
                     daxpy(J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 );
                     daxpy(J1-1, -X( J2 ), T( 1, J2 ), 1, X, 1 );

                     daxpy(J1-1, -X( N+J1 ), T( 1, J1 ), 1, X( N+1 ), 1 );
                     daxpy(J1-1, -X( N+J2 ), T( 1, J2 ), 1, X( N+1 ), 1 );

                     X( 1 ) = X( 1 ) + B( J1 )*X( N+J1 ) + B( J2 )*X( N+J2 )                      X( N+1 ) = X( N+1 ) - B( J1 )*X( J1 ) - B( J2 )*X( J2 )

                     XMAX = ZERO
                     for (K = 1; K <= J1 - 1; K++) { // 60
                        XMAX = MAX( ABS( X( K ) )+ABS( X( K+N ) ), XMAX )
                     } // 60
                  }

               }
            } // 70

         } else {

            // Solve (T + iB)**T*(p+iq) = c+id

            JNEXT = 1
            for (J = 1; J <= N; J++) { // 80
               if (J.LT.JNEXT) GO TO 80;
               J1 = J
               J2 = J
               JNEXT = J + 1
               if ( J.LT.N ) {
                  if ( T( J+1, J ) != ZERO ) {
                     J2 = J + 1
                     JNEXT = J + 2
                  }
               }

               if ( J1 == J2 ) {

                  // 1 by 1 diagonal block

                  // Scale if necessary to avoid overflow in forming the
                  // right-hand side element by inner product.

                  XJ = ABS( X( J1 ) ) + ABS( X( J1+N ) )
                  if ( XMAX.GT.ONE ) {
                     REC = ONE / XMAX
                     if ( WORK( J1 ).GT.( BIGNUM-XJ )*REC ) {
                        dscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }

                  X( J1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, 1 )
                  X( N+J1 ) = X( N+J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X( N+1 ), 1 )
                  if ( J1.GT.1 ) {
                     X( J1 ) = X( J1 ) - B( J1 )*X( N+1 )
                     X( N+J1 ) = X( N+J1 ) + B( J1 )*X( 1 )
                  }
                  XJ = ABS( X( J1 ) ) + ABS( X( J1+N ) )

                  Z = W
                  if (J1 == 1) Z = B( 1 );

                  // Scale if necessary to avoid overflow in
                  // complex division

                  TJJ = ABS( T( J1, J1 ) ) + ABS( Z )
                  TMP = T( J1, J1 )
                  if ( TJJ.LT.SMINW ) {
                     TMP = SMINW
                     TJJ = SMINW
                     INFO = 1
                  }

                  if ( TJJ.LT.ONE ) {
                     if ( XJ.GT.BIGNUM*TJJ ) {
                        REC = ONE / XJ
                        dscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }
                  dladiv(X( J1 ), X( N+J1 ), TMP, -Z, SR, SI );
                  X( J1 ) = SR
                  X( J1+N ) = SI
                  XMAX = MAX( ABS( X( J1 ) )+ABS( X( J1+N ) ), XMAX )

               } else {

                  // 2 by 2 diagonal block

                  // Scale if necessary to avoid overflow in forming the
                  // right-hand side element by inner product.

                  XJ = MAX( ABS( X( J1 ) )+ABS( X( N+J1 ) ), ABS( X( J2 ) )+ABS( X( N+J2 ) ) )
                  if ( XMAX.GT.ONE ) {
                     REC = ONE / XMAX
                     if ( MAX( WORK( J1 ), WORK( J2 ) ).GT. ( BIGNUM-XJ ) / XMAX ) {
                        dscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     }
                  }

                  D( 1, 1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, 1 )                   D( 2, 1 ) = X( J2 ) - DDOT( J1-1, T( 1, J2 ), 1, X, 1 )                   D( 1, 2 ) = X( N+J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X( N+1 ), 1 )                   D( 2, 2 ) = X( N+J2 ) - DDOT( J1-1, T( 1, J2 ), 1, X( N+1 ), 1 )
                  D( 1, 1 ) = D( 1, 1 ) - B( J1 )*X( N+1 )
                  D( 2, 1 ) = D( 2, 1 ) - B( J2 )*X( N+1 )
                  D( 1, 2 ) = D( 1, 2 ) + B( J1 )*X( 1 )
                  D( 2, 2 ) = D( 2, 2 ) + B( J2 )*X( 1 )

                  dlaln2( true , 2, 2, SMINW, ONE, T( J1, J1 ), LDT, ONE, ONE, D, 2, ZERO, W, V, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 2;

                  if ( SCALOC != ONE ) {
                     dscal(N2, SCALOC, X, 1 );
                     SCALE = SCALOC*SCALE
                  }
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  X( N+J1 ) = V( 1, 2 )
                  X( N+J2 ) = V( 2, 2 )
                  XMAX = MAX( ABS( X( J1 ) )+ABS( X( N+J1 ) ), ABS( X( J2 ) )+ABS( X( N+J2 ) ), XMAX )

               }

            } // 80

         }

      }

      RETURN

      // End of DLAQTR

      }
