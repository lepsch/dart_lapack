      void spbrfs(UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      REAL               ZERO;
      const              ZERO = 0.0 ;
      REAL               ONE;
      const              ONE = 1.0 ;
      REAL               TWO;
      const              TWO = 2.0 ;
      REAL               THREE;
      const              THREE = 3.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                COUNT, I, J, K, KASE, L, NZ;
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SLACN2, SPBTRS, SSBMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, SLAMCH
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDAB < KD+1 ) {
         INFO = -6;
      } else if ( LDAFB < KD+1 ) {
         INFO = -8;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -10;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -12;
      }
      if ( INFO != 0 ) {
         xerbla('SPBRFS', -INFO );
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

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = min( N+1, 2*KD+2 );
      EPS = SLAMCH( 'Epsilon' );
      SAFMIN = SLAMCH( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 140

         COUNT = 1;
         LSTRES = THREE;
         } // 20

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - A * X

         scopy(N, B( 1, J ), 1, WORK( N+1 ), 1 );
         ssbmv(UPLO, N, KD, -ONE, AB, LDAB, X( 1, J ), 1, ONE, WORK( N+1 ), 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 30
            WORK[I] = ( B( I, J ) ).abs();
         } // 30

         // Compute abs(A)*abs(X) + abs(B).

         if ( UPPER ) {
            for (K = 1; K <= N; K++) { // 50
               S = ZERO;
               XK = ( X( K, J ) ).abs();
               L = KD + 1 - K;
               for (I = max( 1, K-KD ); I <= K - 1; I++) { // 40
                  WORK[I] = WORK( I ) + ( AB( L+I, K ) ).abs()*XK;
                  S = S + ( AB( L+I, K ) ).abs()*( X( I, J ) ).abs();
               } // 40
               WORK[K] = WORK( K ) + ( AB( KD+1, K ) ).abs()*XK + S;
            } // 50
         } else {
            for (K = 1; K <= N; K++) { // 70
               S = ZERO;
               XK = ( X( K, J ) ).abs();
               WORK[K] = WORK( K ) + ( AB( 1, K ) ).abs()*XK;
               L = 1 - K;
               for (I = K + 1; I <= min( N, K+KD ); I++) { // 60
                  WORK[I] = WORK( I ) + ( AB( L+I, K ) ).abs()*XK;
                  S = S + ( AB( L+I, K ) ).abs()*( X( I, J ) ).abs();
               } // 60
               WORK[K] = WORK( K ) + S;
            } // 70
         }
         S = ZERO;
         for (I = 1; I <= N; I++) { // 80
            if ( WORK( I ) > SAFE2 ) {
               S = max( S, ( WORK( N+I ) ).abs() / WORK( I ) );
            } else {
               S = max( S, ( ( WORK( N+I ) ).abs()+SAFE1 ) / ( WORK( I )+SAFE1 ) );
            }
         } // 80
         BERR[J] = S;

         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.

         if ( BERR( J ) > EPS && TWO*BERR( J ) <= LSTRES && COUNT <= ITMAX ) {

            // Update solution and try again.

            spbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK( N+1 ), N, INFO );
            saxpy(N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 );
            LSTRES = BERR( J );
            COUNT = COUNT + 1;
            GO TO 20;
         }

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) <= FERR =
         // norm( abs(inv(A))*
            // ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)

         // where
           // norm(Z) is the magnitude of the largest component of Z
           // inv(A) is the inverse of A
           // abs(Z) is the componentwise absolute value of the matrix or
              // vector Z
           // NZ is the maximum number of nonzeros in any row of A, plus 1
           // EPS is machine epsilon

         // The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(A)*abs(X) + abs(B) is less than SAFE2.

         // Use SLACN2 to estimate the infinity-norm of the matrix
            // inv(A) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 90
            if ( WORK( I ) > SAFE2 ) {
               WORK[I] = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I );
            } else {
               WORK[I] = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I ) + SAFE1;
            }
         } // 90

         KASE = 0;
         } // 100
         slacn2(N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(A**T).

               spbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK( N+1 ), N, INFO );
               for (I = 1; I <= N; I++) { // 110
                  WORK[N+I] = WORK( N+I )*WORK( I );
               } // 110
            } else if ( KASE == 2 ) {

               // Multiply by inv(A)*diag(W).

               for (I = 1; I <= N; I++) { // 120
                  WORK[N+I] = WORK( N+I )*WORK( I );
               } // 120
               spbtrs(UPLO, N, KD, 1, AFB, LDAFB, WORK( N+1 ), N, INFO );
            }
            GO TO 100;
         }

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 130
            LSTRES = max( LSTRES, ( X( I, J ) ) ).abs();
         } // 130
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 140

      return;
      }
