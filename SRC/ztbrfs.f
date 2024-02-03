      SUBROUTINE ZTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, KD, LDAB, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), FERR( * ), RWORK( * );
      COMPLEX*16         AB( LDAB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      String             TRANSN, TRANST;
      int                I, J, K, KASE, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      COMPLEX*16         ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZLACN2, ZTBMV, ZTBSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
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
      } else if ( KD.LT.0 ) {
         INFO = -5
      } else if ( NRHS.LT.0 ) {
         INFO = -6
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -8
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -10
      } else if ( LDX.LT.MAX( 1, N ) ) {
         INFO = -12
      }
      if ( INFO.NE.0 ) {
         xerbla('ZTBRFS', -INFO );
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
         TRANSN = 'N'
         TRANST = 'C'
      } else {
         TRANSN = 'C'
         TRANST = 'N'
      }

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = KD + 2
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS

      // Do for each right hand side

      DO 250 J = 1, NRHS

         // Compute residual R = B - op(A) * X,
         // where op(A) = A, A**T, or A**H, depending on TRANS.

         zcopy(N, X( 1, J ), 1, WORK, 1 );
         ztbmv(UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 );
         zaxpy(N, -ONE, B( 1, J ), 1, WORK, 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         DO 20 I = 1, N
            RWORK( I ) = CABS1( B( I, J ) )
   20    CONTINUE

         if ( NOTRAN ) {

            // Compute abs(A)*abs(X) + abs(B).

            if ( UPPER ) {
               if ( NOUNIT ) {
                  DO 40 K = 1, N
                     XK = CABS1( X( K, J ) )
                     DO 30 I = MAX( 1, K-KD ), K
                        RWORK( I ) = RWORK( I ) + CABS1( AB( KD+1+I-K, K ) )*XK
   30                CONTINUE
   40             CONTINUE
               } else {
                  DO 60 K = 1, N
                     XK = CABS1( X( K, J ) )
                     DO 50 I = MAX( 1, K-KD ), K - 1
                        RWORK( I ) = RWORK( I ) + CABS1( AB( KD+1+I-K, K ) )*XK
   50                CONTINUE
                     RWORK( K ) = RWORK( K ) + XK
   60             CONTINUE
               }
            } else {
               if ( NOUNIT ) {
                  DO 80 K = 1, N
                     XK = CABS1( X( K, J ) )
                     DO 70 I = K, MIN( N, K+KD )
                        RWORK( I ) = RWORK( I ) + CABS1( AB( 1+I-K, K ) )*XK
   70                CONTINUE
   80             CONTINUE
               } else {
                  DO 100 K = 1, N
                     XK = CABS1( X( K, J ) )
                     DO 90 I = K + 1, MIN( N, K+KD )
                        RWORK( I ) = RWORK( I ) + CABS1( AB( 1+I-K, K ) )*XK
   90                CONTINUE
                     RWORK( K ) = RWORK( K ) + XK
  100             CONTINUE
               }
            }
         } else {

            // Compute abs(A**H)*abs(X) + abs(B).

            if ( UPPER ) {
               if ( NOUNIT ) {
                  DO 120 K = 1, N
                     S = ZERO
                     DO 110 I = MAX( 1, K-KD ), K
                        S = S + CABS1( AB( KD+1+I-K, K ) )* CABS1( X( I, J ) )
  110                CONTINUE
                     RWORK( K ) = RWORK( K ) + S
  120             CONTINUE
               } else {
                  DO 140 K = 1, N
                     S = CABS1( X( K, J ) )
                     DO 130 I = MAX( 1, K-KD ), K - 1
                        S = S + CABS1( AB( KD+1+I-K, K ) )* CABS1( X( I, J ) )
  130                CONTINUE
                     RWORK( K ) = RWORK( K ) + S
  140             CONTINUE
               }
            } else {
               if ( NOUNIT ) {
                  DO 160 K = 1, N
                     S = ZERO
                     DO 150 I = K, MIN( N, K+KD )
                        S = S + CABS1( AB( 1+I-K, K ) )* CABS1( X( I, J ) )
  150                CONTINUE
                     RWORK( K ) = RWORK( K ) + S
  160             CONTINUE
               } else {
                  DO 180 K = 1, N
                     S = CABS1( X( K, J ) )
                     DO 170 I = K + 1, MIN( N, K+KD )
                        S = S + CABS1( AB( 1+I-K, K ) )* CABS1( X( I, J ) )
  170                CONTINUE
                     RWORK( K ) = RWORK( K ) + S
  180             CONTINUE
               }
            }
         }
         S = ZERO
         DO 190 I = 1, N
            if ( RWORK( I ).GT.SAFE2 ) {
               S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) )
            } else {
               S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) )
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

         // Use ZLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         DO 200 I = 1, N
            if ( RWORK( I ).GT.SAFE2 ) {
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I )
            } else {
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1
            }
  200    CONTINUE

         KASE = 0
  210    CONTINUE
         zlacn2(N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE );
         if ( KASE.NE.0 ) {
            if ( KASE.EQ.1 ) {

               // Multiply by diag(W)*inv(op(A)**H).

               ztbsv(UPLO, TRANST, DIAG, N, KD, AB, LDAB, WORK, 1 );
               DO 220 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
  220          CONTINUE
            } else {

               // Multiply by inv(op(A))*diag(W).

               DO 230 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
  230          CONTINUE
               ztbsv(UPLO, TRANSN, DIAG, N, KD, AB, LDAB, WORK, 1 );
            }
            GO TO 210
         }

         // Normalize error.

         LSTRES = ZERO
         DO 240 I = 1, N
            LSTRES = MAX( LSTRES, CABS1( X( I, J ) ) )
  240    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES

  250 CONTINUE

      RETURN

      // End of ZTBRFS

      }
