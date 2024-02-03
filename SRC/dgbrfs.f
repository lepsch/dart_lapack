      SUBROUTINE DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      double             ONE;
      const              ONE = 1.0D+0 ;
      double             TWO;
      const              TWO = 2.0D+0 ;
      double             THREE;
      const              THREE = 3.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      String             TRANST;
      int                COUNT, I, J, K, KASE, KK, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGBMV, DGBTRS, DLACN2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      if ( .NOT.NOTRAN && .NOT.LSAME( TRANS, 'T' ) && .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( KL < 0 ) {
         INFO = -3
      } else if ( KU < 0 ) {
         INFO = -4
      } else if ( NRHS < 0 ) {
         INFO = -5
      } else if ( LDAB < KL+KU+1 ) {
         INFO = -7
      } else if ( LDAFB < 2*KL+KU+1 ) {
         INFO = -9
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -12
      } else if ( LDX < MAX( 1, N ) ) {
         INFO = -14
      }
      if ( INFO != 0 ) {
         xerbla('DGBRFS', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N == 0 || NRHS == 0 ) {
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

      NZ = MIN( KL+KU+2, N+1 )
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 140

         COUNT = 1
         LSTRES = THREE
         } // 20

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - op(A) * X,
         // where op(A) = A, A**T, or A**H, depending on TRANS.

         dcopy(N, B( 1, J ), 1, WORK( N+1 ), 1 );
         dgbmv(TRANS, N, N, KL, KU, -ONE, AB, LDAB, X( 1, J ), 1, ONE, WORK( N+1 ), 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 30
            WORK( I ) = ABS( B( I, J ) )
         } // 30

         // Compute abs(op(A))*abs(X) + abs(B).

         if ( NOTRAN ) {
            for (K = 1; K <= N; K++) { // 50
               KK = KU + 1 - K
               XK = ABS( X( K, J ) )
               DO 40 I = MAX( 1, K-KU ), MIN( N, K+KL )
                  WORK( I ) = WORK( I ) + ABS( AB( KK+I, K ) )*XK
               } // 40
            } // 50
         } else {
            for (K = 1; K <= N; K++) { // 70
               S = ZERO
               KK = KU + 1 - K
               DO 60 I = MAX( 1, K-KU ), MIN( N, K+KL )
                  S = S + ABS( AB( KK+I, K ) )*ABS( X( I, J ) )
               } // 60
               WORK( K ) = WORK( K ) + S
            } // 70
         }
         S = ZERO
         for (I = 1; I <= N; I++) { // 80
            if ( WORK( I ) > SAFE2 ) {
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            } else {
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) / ( WORK( I )+SAFE1 ) )
            }
         } // 80
         BERR( J ) = S

         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.

         if ( BERR( J ) > EPS && TWO*BERR( J ).LE.LSTRES && COUNT.LE.ITMAX ) {

            // Update solution and try again.

            dgbtrs(TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK( N+1 ), N, INFO );
            daxpy(N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 );
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         }

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

         // Use DLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 90
            if ( WORK( I ) > SAFE2 ) {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            } else {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            }
         } // 90

         KASE = 0
         } // 100
         dlacn2(N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(op(A)**T).

               dgbtrs(TRANST, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK( N+1 ), N, INFO );
               for (I = 1; I <= N; I++) { // 110
                  WORK( N+I ) = WORK( N+I )*WORK( I )
               } // 110
            } else {

               // Multiply by inv(op(A))*diag(W).

               for (I = 1; I <= N; I++) { // 120
                  WORK( N+I ) = WORK( N+I )*WORK( I )
               } // 120
               dgbtrs(TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK( N+1 ), N, INFO );
            }
            GO TO 100
         }

         // Normalize error.

         LSTRES = ZERO
         for (I = 1; I <= N; I++) { // 130
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
         } // 130
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 140

      RETURN

      // End of DGBRFS

      }
