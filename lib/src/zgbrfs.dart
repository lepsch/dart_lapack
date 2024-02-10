      void zgbrfs(TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS;
      int                IPIV( * );
      double             BERR( * ), FERR( * ), RWORK( * );
      Complex         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO;
      const              ZERO = 0.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      double             TWO;
      const              TWO = 2.0 ;
      double             THREE;
      const              THREE = 3.0 ;
      bool               NOTRAN;
      String             TRANSN, TRANST;
      int                COUNT, I, J, K, KASE, KK, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      Complex         ZDUM;
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZGBMV, ZGBTRS, ZLACN2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
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
      NOTRAN = lsame( TRANS, 'N' );
      if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KL < 0 ) {
         INFO = -3;
      } else if ( KU < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDAB < KL+KU+1 ) {
         INFO = -7;
      } else if ( LDAFB < 2*KL+KU+1 ) {
         INFO = -9;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -12;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -14;
      }
      if ( INFO != 0 ) {
         xerbla('ZGBRFS', -INFO );
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

      NZ = min( KL+KU+2, N+1 );
      EPS = dlamch( 'Epsilon' );
      SAFMIN = dlamch( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 140

         COUNT = 1;
         LSTRES = THREE;
         } // 20

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - op(A) * X,
         // where op(A) = A, A**T, or A**H, depending on TRANS.

         zcopy(N, B( 1, J ), 1, WORK, 1 );
         zgbmv(TRANS, N, N, KL, KU, -CONE, AB, LDAB, X( 1, J ), 1, CONE, WORK, 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 30
            RWORK[I] = CABS1( B( I, J ) );
         } // 30

         // Compute abs(op(A))*abs(X) + abs(B).

         if ( NOTRAN ) {
            for (K = 1; K <= N; K++) { // 50
               KK = KU + 1 - K;
               XK = CABS1( X( K, J ) );
               for (I = max( 1, K-KU ); I <= min( N, K+KL ); I++) { // 40
                  RWORK[I] = RWORK( I ) + CABS1( AB( KK+I, K ) )*XK;
               } // 40
            } // 50
         } else {
            for (K = 1; K <= N; K++) { // 70
               S = ZERO;
               KK = KU + 1 - K;
               for (I = max( 1, K-KU ); I <= min( N, K+KL ); I++) { // 60
                  S = S + CABS1( AB( KK+I, K ) )*CABS1( X( I, J ) );
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
         //    1) The residual BERR(J) is larger than machine epsilon, and
         //    2) BERR(J) decreased by at least a factor of 2 during the
         //       last iteration, and
         //    3) At most ITMAX iterations tried.

         if ( BERR( J ) > EPS && TWO*BERR( J ) <= LSTRES && COUNT <= ITMAX ) {

            // Update solution and try again.

            zgbtrs(TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            zaxpy(N, CONE, WORK, 1, X( 1, J ), 1 );
            LSTRES = BERR( J );
            COUNT = COUNT + 1;
            GO TO 20;
         }

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) <= FERR =
         // norm( abs(inv(op(A)))*
         //    ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)

         // where
         //   norm(Z) is the magnitude of the largest component of Z
         //   inv(op(A)) is the inverse of op(A)
         //   abs(Z) is the componentwise absolute value of the matrix or
         //      vector Z
         //   NZ is the maximum number of nonzeros in any row of A, plus 1
         //   EPS is machine epsilon

         // The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(op(A))*abs(X) + abs(B) is less than SAFE2.

         // Use ZLACN2 to estimate the infinity-norm of the matrix
         //    inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 90
            if ( RWORK( I ) > SAFE2 ) {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I );
            } else {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1;
            }
         } // 90

         KASE = 0;
         } // 100
         zlacn2(N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(op(A)**H).

               zgbtrs(TRANST, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
               for (I = 1; I <= N; I++) { // 110
                  WORK[I] = RWORK( I )*WORK( I );
               } // 110
            } else {

               // Multiply by inv(op(A))*diag(W).

               for (I = 1; I <= N; I++) { // 120
                  WORK[I] = RWORK( I )*WORK( I );
               } // 120
               zgbtrs(TRANSN, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
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

      }
