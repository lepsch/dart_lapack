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
      PARAMETER          ( ZERO = 0.0D+0 )
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
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

      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( KD.LT.0 ) THEN
         INFO = -5
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTBRFS', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF

      IF( NOTRAN ) THEN
         TRANSN = 'N'
         TRANST = 'C'
      ELSE
         TRANSN = 'C'
         TRANST = 'N'
      END IF

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

         CALL ZCOPY( N, X( 1, J ), 1, WORK, 1 )
         CALL ZTBMV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 )
         CALL ZAXPY( N, -ONE, B( 1, J ), 1, WORK, 1 )

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
        t // han SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         DO 20 I = 1, N
            RWORK( I ) = CABS1( B( I, J ) )
   20    CONTINUE

         IF( NOTRAN ) THEN

            // Compute abs(A)*abs(X) + abs(B).

            IF( UPPER ) THEN
               IF( NOUNIT ) THEN
                  DO 40 K = 1, N
                     XK = CABS1( X( K, J ) )
                     DO 30 I = MAX( 1, K-KD ), K
                        RWORK( I ) = RWORK( I ) + CABS1( AB( KD+1+I-K, K ) )*XK
   30                CONTINUE
   40             CONTINUE
               ELSE
                  DO 60 K = 1, N
                     XK = CABS1( X( K, J ) )
                     DO 50 I = MAX( 1, K-KD ), K - 1
                        RWORK( I ) = RWORK( I ) + CABS1( AB( KD+1+I-K, K ) )*XK
   50                CONTINUE
                     RWORK( K ) = RWORK( K ) + XK
   60             CONTINUE
               END IF
            ELSE
               IF( NOUNIT ) THEN
                  DO 80 K = 1, N
                     XK = CABS1( X( K, J ) )
                     DO 70 I = K, MIN( N, K+KD )
                        RWORK( I ) = RWORK( I ) + CABS1( AB( 1+I-K, K ) )*XK
   70                CONTINUE
   80             CONTINUE
               ELSE
                  DO 100 K = 1, N
                     XK = CABS1( X( K, J ) )
                     DO 90 I = K + 1, MIN( N, K+KD )
                        RWORK( I ) = RWORK( I ) + CABS1( AB( 1+I-K, K ) )*XK
   90                CONTINUE
                     RWORK( K ) = RWORK( K ) + XK
  100             CONTINUE
               END IF
            END IF
         ELSE

            // Compute abs(A**H)*abs(X) + abs(B).

            IF( UPPER ) THEN
               IF( NOUNIT ) THEN
                  DO 120 K = 1, N
                     S = ZERO
                     DO 110 I = MAX( 1, K-KD ), K
                        S = S + CABS1( AB( KD+1+I-K, K ) )* CABS1( X( I, J ) )
  110                CONTINUE
                     RWORK( K ) = RWORK( K ) + S
  120             CONTINUE
               ELSE
                  DO 140 K = 1, N
                     S = CABS1( X( K, J ) )
                     DO 130 I = MAX( 1, K-KD ), K - 1
                        S = S + CABS1( AB( KD+1+I-K, K ) )* CABS1( X( I, J ) )
  130                CONTINUE
                     RWORK( K ) = RWORK( K ) + S
  140             CONTINUE
               END IF
            ELSE
               IF( NOUNIT ) THEN
                  DO 160 K = 1, N
                     S = ZERO
                     DO 150 I = K, MIN( N, K+KD )
                        S = S + CABS1( AB( 1+I-K, K ) )* CABS1( X( I, J ) )
  150                CONTINUE
                     RWORK( K ) = RWORK( K ) + S
  160             CONTINUE
               ELSE
                  DO 180 K = 1, N
                     S = CABS1( X( K, J ) )
                     DO 170 I = K + 1, MIN( N, K+KD )
                        S = S + CABS1( AB( 1+I-K, K ) )* CABS1( X( I, J ) )
  170                CONTINUE
                     RWORK( K ) = RWORK( K ) + S
  180             CONTINUE
               END IF
            END IF
         END IF
         S = ZERO
         DO 190 I = 1, N
            IF( RWORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) )
            ELSE
               S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) )
            END IF
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
            IF( RWORK( I ).GT.SAFE2 ) THEN
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I )
            ELSE
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1
            END IF
  200    CONTINUE

         KASE = 0
  210    CONTINUE
         CALL ZLACN2( N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN

               // Multiply by diag(W)*inv(op(A)**H).

               CALL ZTBSV( UPLO, TRANST, DIAG, N, KD, AB, LDAB, WORK, 1 )
               DO 220 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
  220          CONTINUE
            ELSE

               // Multiply by inv(op(A))*diag(W).

               DO 230 I = 1, N
                  WORK( I ) = RWORK( I )*WORK( I )
  230          CONTINUE
               CALL ZTBSV( UPLO, TRANSN, DIAG, N, KD, AB, LDAB, WORK, 1 )
            END IF
            GO TO 210
         END IF

         // Normalize error.

         LSTRES = ZERO
         DO 240 I = 1, N
            LSTRES = MAX( LSTRES, CABS1( X( I, J ) ) )
  240    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES

  250 CONTINUE

      RETURN

      // End of ZTBRFS

      END
