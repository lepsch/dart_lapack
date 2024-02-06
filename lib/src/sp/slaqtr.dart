      void slaqtr(LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               LREAL, LTRAN;
      int                INFO, LDT, N;
      double               SCALE, W;
      // ..
      // .. Array Arguments ..
      double               B( * ), T( LDT, * ), WORK( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      int                I, IERR, J, J1, J2, JNEXT, K, N1, N2;
      double               BIGNUM, EPS, REC, SCALOC, SI, SMIN, SMINW, SMLNUM, SR, TJJ, TMP, XJ, XMAX, XNORM, Z;
      // ..
      // .. Local Arrays ..
      double               D( 2, 2 ), V( 2, 2 );
      // ..
      // .. External Functions ..
      //- int                ISAMAX;
      //- REAL               SASUM, SDOT, SLAMCH, SLANGE;
      // EXTERNAL ISAMAX, SASUM, SDOT, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SLADIV, SLALN2, SSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Do not test the input parameters for errors

      NOTRAN = !LTRAN;
      INFO = 0;

      // Quick return if possible

      if (N == 0) return;

      // Set constants to control overflow

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      XNORM = SLANGE( 'M', N, N, T, LDT, D );
      if ( !LREAL) XNORM = max( XNORM, ( W ).abs(), SLANGE( 'M', N, 1, B, N, D ) );
      SMIN = max( SMLNUM, EPS*XNORM );

      // Compute 1-norm of each column of strictly upper triangular
      // part of T to control overflow in triangular solver.

      WORK[1] = ZERO;
      for (J = 2; J <= N; J++) { // 10
         WORK[J] = SASUM( J-1, T( 1, J ), 1 );
      } // 10

      if ( !LREAL ) {
         for (I = 2; I <= N; I++) { // 20
            WORK[I] = WORK( I ) + ( B( I ) ).abs();
         } // 20
      }

      N2 = 2*N;
      N1 = N;
      if ( !LREAL) N1 = N2;
      K = ISAMAX( N1, X, 1 );
      XMAX = ( X( K ) ).abs();
      SCALE = ONE;

      if ( XMAX > BIGNUM ) {
         SCALE = BIGNUM / XMAX;
         sscal(N1, SCALE, X, 1 );
         XMAX = BIGNUM;
      }

      if ( LREAL ) {

         if ( NOTRAN ) {

            // Solve T*p = scale*c

            JNEXT = N;
            for (J = N; J >= 1; J--) { // 30
               if (J > JNEXT) GO TO 30;
               J1 = J;
               J2 = J;
               JNEXT = J - 1;
               if ( J > 1 ) {
                  if ( T( J, J-1 ) != ZERO ) {
                     J1 = J - 1;
                     JNEXT = J - 2;
                  }
               }

               if ( J1 == J2 ) {

                  // Meet 1 by 1 diagonal block

                  // Scale to avoid overflow when computing
                      // x(j) = b(j)/T(j,j)

                  XJ = ( X( J1 ) ).abs();
                  TJJ = ( T( J1, J1 ) ).abs();
                  TMP = T( J1, J1 );
                  if ( TJJ < SMIN ) {
                     TMP = SMIN;
                     TJJ = SMIN;
                     INFO = 1;
                  }

                  if (XJ == ZERO) GO TO 30;

                  if ( TJJ < ONE ) {
                     if ( XJ > BIGNUM*TJJ ) {
                        REC = ONE / XJ;
                        sscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }
                  X[J1] = X( J1 ) / TMP;
                  XJ = ( X( J1 ) ).abs();

                  // Scale x if necessary to avoid overflow when adding a
                  // multiple of column j1 of T.

                  if ( XJ > ONE ) {
                     REC = ONE / XJ;
                     if ( WORK( J1 ) > ( BIGNUM-XMAX )*REC ) {
                        sscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                     }
                  }
                  if ( J1 > 1 ) {
                     saxpy(J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 );
                     K = ISAMAX( J1-1, X, 1 );
                     XMAX = ( X( K ) ).abs();
                  }

               } else {

                  // Meet 2 by 2 diagonal block

                  // Call 2 by 2 linear system solve, to take
                  // care of possible overflow by scaling factor.

                  D[1][1] = X( J1 );
                  D[2][1] = X( J2 );
                  slaln2( false , 2, 1, SMIN, ONE, T( J1, J1 ), LDT, ONE, ONE, D, 2, ZERO, ZERO, V, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 2;

                  if ( SCALOC != ONE ) {
                     sscal(N, SCALOC, X, 1 );
                     SCALE = SCALE*SCALOC;
                  }
                  X[J1] = V( 1, 1 );
                  X[J2] = V( 2, 1 );

                  // Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
                  // to avoid overflow in updating right-hand side.

                  XJ = max( ( V( 1, 1 ) ).abs(), ( V( 2, 1 ) ) ).abs();
                  if ( XJ > ONE ) {
                     REC = ONE / XJ;
                     if ( max( WORK( J1 ), WORK( J2 ) ) > ( BIGNUM-XMAX )*REC ) {
                        sscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                     }
                  }

                  // Update right-hand side

                  if ( J1 > 1 ) {
                     saxpy(J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 );
                     saxpy(J1-1, -X( J2 ), T( 1, J2 ), 1, X, 1 );
                     K = ISAMAX( J1-1, X, 1 );
                     XMAX = ( X( K ) ).abs();
                  }

               }

            } // 30

         } else {

            // Solve T**T*p = scale*c

            JNEXT = 1;
            for (J = 1; J <= N; J++) { // 40
               if (J < JNEXT) GO TO 40;
               J1 = J;
               J2 = J;
               JNEXT = J + 1;
               if ( J < N ) {
                  if ( T( J+1, J ) != ZERO ) {
                     J2 = J + 1;
                     JNEXT = J + 2;
                  }
               }

               if ( J1 == J2 ) {

                  // 1 by 1 diagonal block

                  // Scale if necessary to avoid overflow in forming the
                  // right-hand side element by inner product.

                  XJ = ( X( J1 ) ).abs();
                  if ( XMAX > ONE ) {
                     REC = ONE / XMAX;
                     if ( WORK( J1 ) > ( BIGNUM-XJ )*REC ) {
                        sscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }

                  X[J1] = X( J1 ) - SDOT( J1-1, T( 1, J1 ), 1, X, 1 );

                  XJ = ( X( J1 ) ).abs();
                  TJJ = ( T( J1, J1 ) ).abs();
                  TMP = T( J1, J1 );
                  if ( TJJ < SMIN ) {
                     TMP = SMIN;
                     TJJ = SMIN;
                     INFO = 1;
                  }

                  if ( TJJ < ONE ) {
                     if ( XJ > BIGNUM*TJJ ) {
                        REC = ONE / XJ;
                        sscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }
                  X[J1] = X( J1 ) / TMP;
                  XMAX = max( XMAX, ( X( J1 ) ) ).abs();

               } else {

                  // 2 by 2 diagonal block

                  // Scale if necessary to avoid overflow in forming the
                  // right-hand side elements by inner product.

                  XJ = max( ( X( J1 ) ).abs(), ( X( J2 ) ) ).abs();
                  if ( XMAX > ONE ) {
                     REC = ONE / XMAX;
                     if ( max( WORK( J2 ), WORK( J1 ) ) > ( BIGNUM-XJ )* REC ) {
                        sscal(N, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }

                  D[1][1] = X( J1 ) - SDOT( J1-1, T( 1, J1 ), 1, X, 1 )                   D( 2, 1 ) = X( J2 ) - SDOT( J1-1, T( 1, J2 ), 1, X, 1 );

                  slaln2( true , 2, 1, SMIN, ONE, T( J1, J1 ), LDT, ONE, ONE, D, 2, ZERO, ZERO, V, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 2;

                  if ( SCALOC != ONE ) {
                     sscal(N, SCALOC, X, 1 );
                     SCALE = SCALE*SCALOC;
                  }
                  X[J1] = V( 1, 1 );
                  X[J2] = V( 2, 1 );
                  XMAX = max( ( X( J1 ) ).abs(), ( X( J2 ) ).abs(), XMAX );

               }
            } // 40
         }

      } else {

         SMINW = max( EPS*( W ).abs(), SMIN );
         if ( NOTRAN ) {

            // Solve (T + iB)*(p+iq) = c+id

            JNEXT = N;
            for (J = N; J >= 1; J--) { // 70
               if (J > JNEXT) GO TO 70;
               J1 = J;
               J2 = J;
               JNEXT = J - 1;
               if ( J > 1 ) {
                  if ( T( J, J-1 ) != ZERO ) {
                     J1 = J - 1;
                     JNEXT = J - 2;
                  }
               }

               if ( J1 == J2 ) {

                  // 1 by 1 diagonal block

                  // Scale if necessary to avoid overflow in division

                  Z = W;
                  if (J1 == 1) Z = B( 1 );
                  XJ = ( X( J1 ) ).abs() + ( X( N+J1 ) ).abs();
                  TJJ = ( T( J1, J1 ) ).abs() + ( Z ).abs();
                  TMP = T( J1, J1 );
                  if ( TJJ < SMINW ) {
                     TMP = SMINW;
                     TJJ = SMINW;
                     INFO = 1;
                  }

                  if (XJ == ZERO) GO TO 70;

                  if ( TJJ < ONE ) {
                     if ( XJ > BIGNUM*TJJ ) {
                        REC = ONE / XJ;
                        sscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }
                  sladiv(X( J1 ), X( N+J1 ), TMP, Z, SR, SI );
                  X[J1] = SR;
                  X[N+J1] = SI;
                  XJ = ( X( J1 ) ).abs() + ( X( N+J1 ) ).abs();

                  // Scale x if necessary to avoid overflow when adding a
                  // multiple of column j1 of T.

                  if ( XJ > ONE ) {
                     REC = ONE / XJ;
                     if ( WORK( J1 ) > ( BIGNUM-XMAX )*REC ) {
                        sscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC;
                     }
                  }

                  if ( J1 > 1 ) {
                     saxpy(J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 );
                     saxpy(J1-1, -X( N+J1 ), T( 1, J1 ), 1, X( N+1 ), 1 );

                     X[1] = X( 1 ) + B( J1 )*X( N+J1 );
                     X[N+1] = X( N+1 ) - B( J1 )*X( J1 );

                     XMAX = ZERO;
                     for (K = 1; K <= J1 - 1; K++) { // 50
                        XMAX = max( XMAX, ( X( K ) ).abs()+ ( X( K+N ) ) ).abs();
                     } // 50
                  }

               } else {

                  // Meet 2 by 2 diagonal block

                  D[1][1] = X( J1 );
                  D[2][1] = X( J2 );
                  D[1][2] = X( N+J1 );
                  D[2][2] = X( N+J2 );
                  slaln2( false , 2, 2, SMINW, ONE, T( J1, J1 ), LDT, ONE, ONE, D, 2, ZERO, -W, V, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 2;

                  if ( SCALOC != ONE ) {
                     sscal(2*N, SCALOC, X, 1 );
                     SCALE = SCALOC*SCALE;
                  }
                  X[J1] = V( 1, 1 );
                  X[J2] = V( 2, 1 );
                  X[N+J1] = V( 1, 2 );
                  X[N+J2] = V( 2, 2 );

                  // Scale X(J1), .... to avoid overflow in
                  // updating right hand side.

                  XJ = max( ( V( 1, 1 ) ).abs()+( V( 1, 2 ) ).abs(), ( V( 2, 1 ) ).abs()+( V( 2, 2 ) ) ).abs();
                  if ( XJ > ONE ) {
                     REC = ONE / XJ;
                     if ( max( WORK( J1 ), WORK( J2 ) ) > ( BIGNUM-XMAX )*REC ) {
                        sscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC;
                     }
                  }

                  // Update the right-hand side.

                  if ( J1 > 1 ) {
                     saxpy(J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 );
                     saxpy(J1-1, -X( J2 ), T( 1, J2 ), 1, X, 1 );

                     saxpy(J1-1, -X( N+J1 ), T( 1, J1 ), 1, X( N+1 ), 1 );
                     saxpy(J1-1, -X( N+J2 ), T( 1, J2 ), 1, X( N+1 ), 1 );

                     X[1] = X( 1 ) + B( J1 )*X( N+J1 ) + B( J2 )*X( N+J2 )                      X( N+1 ) = X( N+1 ) - B( J1 )*X( J1 ) - B( J2 )*X( J2 );

                     XMAX = ZERO;
                     for (K = 1; K <= J1 - 1; K++) { // 60
                        XMAX = max( ( X( K ) ).abs()+( X( K+N ) ).abs(), XMAX );
                     } // 60
                  }

               }
            } // 70

         } else {

            // Solve (T + iB)**T*(p+iq) = c+id

            JNEXT = 1;
            for (J = 1; J <= N; J++) { // 80
               if (J < JNEXT) GO TO 80;
               J1 = J;
               J2 = J;
               JNEXT = J + 1;
               if ( J < N ) {
                  if ( T( J+1, J ) != ZERO ) {
                     J2 = J + 1;
                     JNEXT = J + 2;
                  }
               }

               if ( J1 == J2 ) {

                  // 1 by 1 diagonal block

                  // Scale if necessary to avoid overflow in forming the
                  // right-hand side element by inner product.

                  XJ = ( X( J1 ) ).abs() + ( X( J1+N ) ).abs();
                  if ( XMAX > ONE ) {
                     REC = ONE / XMAX;
                     if ( WORK( J1 ) > ( BIGNUM-XJ )*REC ) {
                        sscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }

                  X[J1] = X( J1 ) - SDOT( J1-1, T( 1, J1 ), 1, X, 1 );
                  X[N+J1] = X( N+J1 ) - SDOT( J1-1, T( 1, J1 ), 1, X( N+1 ), 1 );
                  if ( J1 > 1 ) {
                     X[J1] = X( J1 ) - B( J1 )*X( N+1 );
                     X[N+J1] = X( N+J1 ) + B( J1 )*X( 1 );
                  }
                  XJ = ( X( J1 ) ).abs() + ( X( J1+N ) ).abs();

                  Z = W;
                  if (J1 == 1) Z = B( 1 );

                  // Scale if necessary to avoid overflow in
                  // complex division

                  TJJ = ( T( J1, J1 ) ).abs() + ( Z ).abs();
                  TMP = T( J1, J1 );
                  if ( TJJ < SMINW ) {
                     TMP = SMINW;
                     TJJ = SMINW;
                     INFO = 1;
                  }

                  if ( TJJ < ONE ) {
                     if ( XJ > BIGNUM*TJJ ) {
                        REC = ONE / XJ;
                        sscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }
                  sladiv(X( J1 ), X( N+J1 ), TMP, -Z, SR, SI );
                  X[J1] = SR;
                  X[J1+N] = SI;
                  XMAX = max( ( X( J1 ) ).abs()+( X( J1+N ) ).abs(), XMAX );

               } else {

                  // 2 by 2 diagonal block

                  // Scale if necessary to avoid overflow in forming the
                  // right-hand side element by inner product.

                  XJ = max( ( X( J1 ) ).abs()+( X( N+J1 ) ).abs(), ( X( J2 ) ).abs()+( X( N+J2 ) ) ).abs();
                  if ( XMAX > ONE ) {
                     REC = ONE / XMAX;
                     if ( max( WORK( J1 ), WORK( J2 ) ) > ( BIGNUM-XJ ) / XMAX ) {
                        sscal(N2, REC, X, 1 );
                        SCALE = SCALE*REC;
                        XMAX = XMAX*REC;
                     }
                  }

                  D[1][1] = X( J1 ) - SDOT( J1-1, T( 1, J1 ), 1, X, 1 )                   D( 2, 1 ) = X( J2 ) - SDOT( J1-1, T( 1, J2 ), 1, X, 1 )                   D( 1, 2 ) = X( N+J1 ) - SDOT( J1-1, T( 1, J1 ), 1, X( N+1 ), 1 )                   D( 2, 2 ) = X( N+J2 ) - SDOT( J1-1, T( 1, J2 ), 1, X( N+1 ), 1 );
                  D[1][1] = D( 1, 1 ) - B( J1 )*X( N+1 );
                  D[2][1] = D( 2, 1 ) - B( J2 )*X( N+1 );
                  D[1][2] = D( 1, 2 ) + B( J1 )*X( 1 );
                  D[2][2] = D( 2, 2 ) + B( J2 )*X( 1 );

                  slaln2( true , 2, 2, SMINW, ONE, T( J1, J1 ), LDT, ONE, ONE, D, 2, ZERO, W, V, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 2;

                  if ( SCALOC != ONE ) {
                     sscal(N2, SCALOC, X, 1 );
                     SCALE = SCALOC*SCALE;
                  }
                  X[J1] = V( 1, 1 );
                  X[J2] = V( 2, 1 );
                  X[N+J1] = V( 1, 2 );
                  X[N+J2] = V( 2, 2 );
                  XMAX = max( ( X( J1 ) ).abs()+( X( N+J1 ) ).abs(), ( X( J2 ) ).abs()+( X( N+J2 ) ).abs(), XMAX );

               }

            } // 80

         }

      }

      return;
      }
