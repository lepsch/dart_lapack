      void cporfs(UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      double               BERR( * ), FERR( * ), RWORK( * );
      Complex            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

// ====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double               ZERO;
      const              ZERO = 0.0 ;
      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      double               TWO;
      const              TWO = 2.0 ;
      double               THREE;
      const              THREE = 3.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                COUNT, I, J, K, KASE, NZ;
      double               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      Complex            ZDUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CHEMV, CLACN2, CPOTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, SLAMCH
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('CPORFS', -INFO );
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

      NZ = N + 1;
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

         ccopy(N, B( 1, J ), 1, WORK, 1 );
         chemv(UPLO, N, -ONE, A, LDA, X( 1, J ), 1, ONE, WORK, 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 30
            RWORK[I] = CABS1( B( I, J ) );
         } // 30

         // Compute abs(A)*abs(X) + abs(B).

         if ( UPPER ) {
            for (K = 1; K <= N; K++) { // 50
               S = ZERO;
               XK = CABS1( X( K, J ) );
               for (I = 1; I <= K - 1; I++) { // 40
                  RWORK[I] = RWORK( I ) + CABS1( A( I, K ) )*XK;
                  S = S + CABS1( A( I, K ) )*CABS1( X( I, J ) );
               } // 40
               RWORK[K] = RWORK( K ) + ABS( double( A( K, K ) ) )*XK + S;
            } // 50
         } else {
            for (K = 1; K <= N; K++) { // 70
               S = ZERO;
               XK = CABS1( X( K, J ) );
               RWORK[K] = RWORK( K ) + ABS( double( A( K, K ) ) )*XK;
               for (I = K + 1; I <= N; I++) { // 60
                  RWORK[I] = RWORK( I ) + CABS1( A( I, K ) )*XK;
                  S = S + CABS1( A( I, K ) )*CABS1( X( I, J ) );
               } // 60
               RWORK[K] = RWORK( K ) + S;
            } // 70
         }
         S = ZERO;
         for (I = 1; I <= N; I++) { // 80
            if ( RWORK( I ) > SAFE2 ) {
               S = max( S, CABS1( WORK( I ) ) / RWORK( I ) );
            } else {
               S = max( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) );
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

            cpotrs(UPLO, N, 1, AF, LDAF, WORK, N, INFO );
            caxpy(N, ONE, WORK, 1, X( 1, J ), 1 );
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

         // Use CLACN2 to estimate the infinity-norm of the matrix
            // inv(A) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 90
            if ( RWORK( I ) > SAFE2 ) {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I );
            } else {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1;
            }
         } // 90

         KASE = 0;
         } // 100
         clacn2(N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(A**H).

               cpotrs(UPLO, N, 1, AF, LDAF, WORK, N, INFO );
               for (I = 1; I <= N; I++) { // 110
                  WORK[I] = RWORK( I )*WORK( I );
               } // 110
            } else if ( KASE == 2 ) {

               // Multiply by inv(A)*diag(W).

               for (I = 1; I <= N; I++) { // 120
                  WORK[I] = RWORK( I )*WORK( I );
               } // 120
               cpotrs(UPLO, N, 1, AF, LDAF, WORK, N, INFO );
            }
            GO TO 100;
         }

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 130
            LSTRES = max( LSTRES, CABS1( X( I, J ) ) );
         } // 130
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 140

      return;
      }
