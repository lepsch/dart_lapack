      SUBROUTINE DTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AP( * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      double             ONE;
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      String             TRANST;
      int                I, J, K, KASE, KC, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DLACN2, DTPMV, DTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )

      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) {
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
      if ( INFO.NE.0 ) {
         xerbla('DTPRFS', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 .OR. NRHS.EQ.0 ) {
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      }

      if ( NOTRAN ) {
         TRANST = 'T'
      } else {
         TRANST = 'N'
      }

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = N + 1
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS

      // Do for each right hand side

      DO 250 J = 1, NRHS

         // Compute residual R = B - op(A) * X,
         // where op(A) = A or A**T, depending on TRANS.

         dcopy(N, X( 1, J ), 1, WORK( N+1 ), 1 );
         dtpmv(UPLO, TRANS, DIAG, N, AP, WORK( N+1 ), 1 );
         daxpy(N, -ONE, B( 1, J ), 1, WORK( N+1 ), 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         DO 20 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   20    CONTINUE

         if ( NOTRAN ) {

            // Compute abs(A)*abs(X) + abs(B).

            if ( UPPER ) {
               KC = 1
               if ( NOUNIT ) {
                  DO 40 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 30 I = 1, K
                        WORK( I ) = WORK( I ) + ABS( AP( KC+I-1 ) )*XK
   30                CONTINUE
                     KC = KC + K
   40             CONTINUE
               } else {
                  DO 60 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 50 I = 1, K - 1
                        WORK( I ) = WORK( I ) + ABS( AP( KC+I-1 ) )*XK
   50                CONTINUE
                     WORK( K ) = WORK( K ) + XK
                     KC = KC + K
   60             CONTINUE
               }
            } else {
               KC = 1
               if ( NOUNIT ) {
                  DO 80 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 70 I = K, N
                        WORK( I ) = WORK( I ) + ABS( AP( KC+I-K ) )*XK
   70                CONTINUE
                     KC = KC + N - K + 1
   80             CONTINUE
               } else {
                  DO 100 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 90 I = K + 1, N
                        WORK( I ) = WORK( I ) + ABS( AP( KC+I-K ) )*XK
   90                CONTINUE
                     WORK( K ) = WORK( K ) + XK
                     KC = KC + N - K + 1
  100             CONTINUE
               }
            }
         } else {

            // Compute abs(A**T)*abs(X) + abs(B).

            if ( UPPER ) {
               KC = 1
               if ( NOUNIT ) {
                  DO 120 K = 1, N
                     S = ZERO
                     DO 110 I = 1, K
                        S = S + ABS( AP( KC+I-1 ) )*ABS( X( I, J ) )
  110                CONTINUE
                     WORK( K ) = WORK( K ) + S
                     KC = KC + K
  120             CONTINUE
               } else {
                  DO 140 K = 1, N
                     S = ABS( X( K, J ) )
                     DO 130 I = 1, K - 1
                        S = S + ABS( AP( KC+I-1 ) )*ABS( X( I, J ) )
  130                CONTINUE
                     WORK( K ) = WORK( K ) + S
                     KC = KC + K
  140             CONTINUE
               }
            } else {
               KC = 1
               if ( NOUNIT ) {
                  DO 160 K = 1, N
                     S = ZERO
                     DO 150 I = K, N
                        S = S + ABS( AP( KC+I-K ) )*ABS( X( I, J ) )
  150                CONTINUE
                     WORK( K ) = WORK( K ) + S
                     KC = KC + N - K + 1
  160             CONTINUE
               } else {
                  DO 180 K = 1, N
                     S = ABS( X( K, J ) )
                     DO 170 I = K + 1, N
                        S = S + ABS( AP( KC+I-K ) )*ABS( X( I, J ) )
  170                CONTINUE
                     WORK( K ) = WORK( K ) + S
                     KC = KC + N - K + 1
  180             CONTINUE
               }
            }
         }
         S = ZERO
         DO 190 I = 1, N
            if ( WORK( I ).GT.SAFE2 ) {
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            } else {
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) / ( WORK( I )+SAFE1 ) )
            }
  190    CONTINUE
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

         // Use DLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         DO 200 I = 1, N
            if ( WORK( I ).GT.SAFE2 ) {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            } else {
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            }
  200    CONTINUE

         KASE = 0
  210    CONTINUE
         dlacn2(N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE );
         if ( KASE.NE.0 ) {
            if ( KASE.EQ.1 ) {

               // Multiply by diag(W)*inv(op(A)**T).

               dtpsv(UPLO, TRANST, DIAG, N, AP, WORK( N+1 ), 1 );
               DO 220 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  220          CONTINUE
            } else {

               // Multiply by inv(op(A))*diag(W).

               DO 230 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  230          CONTINUE
               dtpsv(UPLO, TRANS, DIAG, N, AP, WORK( N+1 ), 1 );
            }
            GO TO 210
         }

         // Normalize error.

         LSTRES = ZERO
         DO 240 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  240    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES

  250 CONTINUE

      RETURN

      // End of DTPRFS

      }
