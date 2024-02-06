      void ztprfs(UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      double             BERR( * ), FERR( * ), RWORK( * );
      Complex         AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               NOTRAN, NOUNIT, UPPER;
      String             TRANSN, TRANST;
      int                I, J, K, KASE, KC, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      Complex         ZDUM;
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZLACN2, ZTPMV, ZTPSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH;
      // EXTERNAL lsame, DLAMCH
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      NOTRAN = lsame( TRANS, 'N' );
      NOUNIT = lsame( DIAG, 'N' );

      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('ZTPRFS', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 || NRHS == 0 ) {
         for (J = 1; J <= NRHS; J++) { // 10
            FERR[J] = ZERO;
            BERR[J] = ZERO;
         } // 10
         return;
      }

      if ( NOTRAN ) {
         TRANSN = 'N';
         TRANST = 'C';
      } else {
         TRANSN = 'C';
         TRANST = 'N';
      }

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = N + 1;
      EPS = dlamch( 'Epsilon' );
      SAFMIN = dlamch( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 250

         // Compute residual R = B - op(A) * X,
         // where op(A) = A, A**T, or A**H, depending on TRANS.

         zcopy(N, X( 1, J ), 1, WORK, 1 );
         ztpmv(UPLO, TRANS, DIAG, N, AP, WORK, 1 );
         zaxpy(N, -ONE, B( 1, J ), 1, WORK, 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 20
            RWORK[I] = CABS1( B( I, J ) );
         } // 20

         if ( NOTRAN ) {

            // Compute abs(A)*abs(X) + abs(B).

            if ( UPPER ) {
               KC = 1;
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 40
                     XK = CABS1( X( K, J ) );
                     for (I = 1; I <= K; I++) { // 30
                        RWORK[I] = RWORK( I ) + CABS1( AP( KC+I-1 ) )*XK;
                     } // 30
                     KC = KC + K;
                  } // 40
               } else {
                  for (K = 1; K <= N; K++) { // 60
                     XK = CABS1( X( K, J ) );
                     for (I = 1; I <= K - 1; I++) { // 50
                        RWORK[I] = RWORK( I ) + CABS1( AP( KC+I-1 ) )*XK;
                     } // 50
                     RWORK[K] = RWORK( K ) + XK;
                     KC = KC + K;
                  } // 60
               }
            } else {
               KC = 1;
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 80
                     XK = CABS1( X( K, J ) );
                     for (I = K; I <= N; I++) { // 70
                        RWORK[I] = RWORK( I ) + CABS1( AP( KC+I-K ) )*XK;
                     } // 70
                     KC = KC + N - K + 1;
                  } // 80
               } else {
                  for (K = 1; K <= N; K++) { // 100
                     XK = CABS1( X( K, J ) );
                     for (I = K + 1; I <= N; I++) { // 90
                        RWORK[I] = RWORK( I ) + CABS1( AP( KC+I-K ) )*XK;
                     } // 90
                     RWORK[K] = RWORK( K ) + XK;
                     KC = KC + N - K + 1;
                  } // 100
               }
            }
         } else {

            // Compute abs(A**H)*abs(X) + abs(B).

            if ( UPPER ) {
               KC = 1;
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 120
                     S = ZERO;
                     for (I = 1; I <= K; I++) { // 110
                        S = S + CABS1( AP( KC+I-1 ) )*CABS1( X( I, J ) );
                     } // 110
                     RWORK[K] = RWORK( K ) + S;
                     KC = KC + K;
                  } // 120
               } else {
                  for (K = 1; K <= N; K++) { // 140
                     S = CABS1( X( K, J ) );
                     for (I = 1; I <= K - 1; I++) { // 130
                        S = S + CABS1( AP( KC+I-1 ) )*CABS1( X( I, J ) );
                     } // 130
                     RWORK[K] = RWORK( K ) + S;
                     KC = KC + K;
                  } // 140
               }
            } else {
               KC = 1;
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 160
                     S = ZERO;
                     for (I = K; I <= N; I++) { // 150
                        S = S + CABS1( AP( KC+I-K ) )*CABS1( X( I, J ) );
                     } // 150
                     RWORK[K] = RWORK( K ) + S;
                     KC = KC + N - K + 1;
                  } // 160
               } else {
                  for (K = 1; K <= N; K++) { // 180
                     S = CABS1( X( K, J ) );
                     for (I = K + 1; I <= N; I++) { // 170
                        S = S + CABS1( AP( KC+I-K ) )*CABS1( X( I, J ) );
                     } // 170
                     RWORK[K] = RWORK( K ) + S;
                     KC = KC + N - K + 1;
                  } // 180
               }
            }
         }
         S = ZERO;
         for (I = 1; I <= N; I++) { // 190
            if ( RWORK( I ) > SAFE2 ) {
               S = max( S, CABS1( WORK( I ) ) / RWORK( I ) );
            } else {
               S = max( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) );
            }
         } // 190
         BERR[J] = S;

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) <= FERR =
         // norm( abs(inv(op(A)))*
            // ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)

         // where
           // norm(Z) is the magnitude of the largest component of Z
           // inv(op(A)) is the inverse of op(A)
           // abs(Z) is the componentwise absolute value of the matrix or
              // vector Z
           // NZ is the maximum number of nonzeros in any row of A, plus 1
           // EPS is machine epsilon

         // The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(op(A))*abs(X) + abs(B) is less than SAFE2.

         // Use ZLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 200
            if ( RWORK( I ) > SAFE2 ) {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I );
            } else {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1;
            }
         } // 200

         KASE = 0;
         } // 210
         zlacn2(N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(op(A)**H).

               ztpsv(UPLO, TRANST, DIAG, N, AP, WORK, 1 );
               for (I = 1; I <= N; I++) { // 220
                  WORK[I] = RWORK( I )*WORK( I );
               } // 220
            } else {

               // Multiply by inv(op(A))*diag(W).

               for (I = 1; I <= N; I++) { // 230
                  WORK[I] = RWORK( I )*WORK( I );
               } // 230
               ztpsv(UPLO, TRANSN, DIAG, N, AP, WORK, 1 );
            }
            GO TO 210;
         }

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 240
            LSTRES = max( LSTRES, CABS1( X( I, J ) ) );
         } // 240
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 250

      return;
      }
