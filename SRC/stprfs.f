      SUBROUTINE STPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AP( * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E+0 ;
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      String             TRANST;
      int                I, J, K, KASE, KC, NZ;
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SLACN2, STPMV, STPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )

      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN && .NOT.LSAME( TRANS, 'T' ) && .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT && .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( NRHS.LT.0 ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LDX.LT.MAX( 1, N ) ) {
         INFO = -10
      }
      if ( INFO != 0 ) {
         xerbla('STPRFS', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N == 0 .OR. NRHS == 0 ) {
         for (J = 1; J <= NRHS; J++) { // 10
            FERR( J ) = ZERO
            BERR( J ) = ZERO
         } // 10
         RETURN
      }

      if ( NOTRAN ) {
         TRANST = 'T'
      } else {
         TRANST = 'N'
      }

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = N + 1
      EPS = SLAMCH( 'Epsilon' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 250

         // Compute residual R = B - op(A) * X,
         // where op(A) = A or A**T, depending on TRANS.

         scopy(N, X( 1, J ), 1, WORK( N+1 ), 1 );
         stpmv(UPLO, TRANS, DIAG, N, AP, WORK( N+1 ), 1 );
         saxpy(N, -ONE, B( 1, J ), 1, WORK( N+1 ), 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 20
            WORK( I ) = ABS( B( I, J ) )
         } // 20

         if ( NOTRAN ) {

            // Compute abs(A)*abs(X) + abs(B).

            if ( UPPER ) {
               KC = 1
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 40
                     XK = ABS( X( K, J ) )
                     for (I = 1; I <= K; I++) { // 30
                        WORK( I ) = WORK( I ) + ABS( AP( KC+I-1 ) )*XK
                     } // 30
                     KC = KC + K
                  } // 40
               } else {
                  for (K = 1; K <= N; K++) { // 60
                     XK = ABS( X( K, J ) )
                     for (I = 1; I <= K - 1; I++) { // 50
                        WORK( I ) = WORK( I ) + ABS( AP( KC+I-1 ) )*XK
                     } // 50
                     WORK( K ) = WORK( K ) + XK
                     KC = KC + K
                  } // 60
               }
            } else {
               KC = 1
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 80
                     XK = ABS( X( K, J ) )
                     for (I = K; I <= N; I++) { // 70
                        WORK( I ) = WORK( I ) + ABS( AP( KC+I-K ) )*XK
                     } // 70
                     KC = KC + N - K + 1
                  } // 80
               } else {
                  for (K = 1; K <= N; K++) { // 100
                     XK = ABS( X( K, J ) )
                     for (I = K + 1; I <= N; I++) { // 90
                        WORK( I ) = WORK( I ) + ABS( AP( KC+I-K ) )*XK
                     } // 90
                     WORK( K ) = WORK( K ) + XK
                     KC = KC + N - K + 1
                  } // 100
               }
            }
         } else {

            // Compute abs(A**T)*abs(X) + abs(B).

            if ( UPPER ) {
               KC = 1
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 120
                     S = ZERO
                     for (I = 1; I <= K; I++) { // 110
                        S = S + ABS( AP( KC+I-1 ) )*ABS( X( I, J ) )
                     } // 110
                     WORK( K ) = WORK( K ) + S
                     KC = KC + K
                  } // 120
               } else {
                  for (K = 1; K <= N; K++) { // 140
                     S = ABS( X( K, J ) )
                     for (I = 1; I <= K - 1; I++) { // 130
                        S = S + ABS( AP( KC+I-1 ) )*ABS( X( I, J ) )
                     } // 130
                     WORK( K ) = WORK( K ) + S
                     KC = KC + K
                  } // 140
               }
            } else {
               KC = 1
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 160
                     S = ZERO
                     for (I = K; I <= N; I++) { // 150
                        S = S + ABS( AP( KC+I-K ) )*ABS( X( I, J ) )
                     } // 150
                     WORK( K ) = WORK( K ) + S
                     KC = KC + N - K + 1
                  } // 160
               } else {
                  for (K = 1; K <= N; K++) { // 180
                     S = ABS( X( K, J ) )
                     for (I = K + 1; I <= N; I++) { // 170
                        S = S + ABS( AP( KC+I-K ) )*ABS( X( I, J ) )
                     } // 170
                     WORK( K ) = WORK( K ) + S
                     KC = KC + N - K + 1
                  } // 180
               }
            }
         }
         S = ZERO
         for (I = 1; I <= N; I++) { // 190
            if ( WORK( I ).GT.SAFE2 ) {
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            } else {
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) / ( WORK( I )+SAFE1 ) )
            }
         } // 190
         BERR( J ) = S

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) .le. FERR =
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

         // Use SLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 200
            if ( WORK( I ).GT.SAFE2 ) {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            } else {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            }
         } // 200

         KASE = 0
         } // 210
         slacn2(N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(op(A)**T).

               stpsv(UPLO, TRANST, DIAG, N, AP, WORK( N+1 ), 1 );
               for (I = 1; I <= N; I++) { // 220
                  WORK( N+I ) = WORK( I )*WORK( N+I )
               } // 220
            } else {

               // Multiply by inv(op(A))*diag(W).

               for (I = 1; I <= N; I++) { // 230
                  WORK( N+I ) = WORK( I )*WORK( N+I )
               } // 230
               stpsv(UPLO, TRANS, DIAG, N, AP, WORK( N+1 ), 1 );
            }
            GO TO 210
         }

         // Normalize error.

         LSTRES = ZERO
         for (I = 1; I <= N; I++) { // 240
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
         } // 240
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 250

      RETURN

      // End of STPRFS

      }
