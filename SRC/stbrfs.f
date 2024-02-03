      SUBROUTINE STBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                INFO, KD, LDAB, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AB( LDAB, * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               NOTRAN, NOUNIT, UPPER;
      String             TRANST;
      int                I, J, K, KASE, NZ;
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SLACN2, STBMV, STBSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
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
         CALL XERBLA( 'STBRFS', -INFO )
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
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = KD + 2
      EPS = SLAMCH( 'Epsilon' )
      SAFMIN = SLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS

      // Do for each right hand side

      DO 250 J = 1, NRHS

         // Compute residual R = B - op(A) * X,
         // where op(A) = A or A**T, depending on TRANS.

         CALL SCOPY( N, X( 1, J ), 1, WORK( N+1 ), 1 )
         CALL STBMV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK( N+1 ), 1 )
         CALL SAXPY( N, -ONE, B( 1, J ), 1, WORK( N+1 ), 1 )

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
        t // han SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         DO 20 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   20    CONTINUE

         IF( NOTRAN ) THEN

            // Compute abs(A)*abs(X) + abs(B).

            IF( UPPER ) THEN
               IF( NOUNIT ) THEN
                  DO 40 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 30 I = MAX( 1, K-KD ), K
                        WORK( I ) = WORK( I ) + ABS( AB( KD+1+I-K, K ) )*XK
   30                CONTINUE
   40             CONTINUE
               ELSE
                  DO 60 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 50 I = MAX( 1, K-KD ), K - 1
                        WORK( I ) = WORK( I ) + ABS( AB( KD+1+I-K, K ) )*XK
   50                CONTINUE
                     WORK( K ) = WORK( K ) + XK
   60             CONTINUE
               END IF
            ELSE
               IF( NOUNIT ) THEN
                  DO 80 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 70 I = K, MIN( N, K+KD )
                        WORK( I ) = WORK( I ) + ABS( AB( 1+I-K, K ) )*XK
   70                CONTINUE
   80             CONTINUE
               ELSE
                  DO 100 K = 1, N
                     XK = ABS( X( K, J ) )
                     DO 90 I = K + 1, MIN( N, K+KD )
                        WORK( I ) = WORK( I ) + ABS( AB( 1+I-K, K ) )*XK
   90                CONTINUE
                     WORK( K ) = WORK( K ) + XK
  100             CONTINUE
               END IF
            END IF
         ELSE

            // Compute abs(A**T)*abs(X) + abs(B).

            IF( UPPER ) THEN
               IF( NOUNIT ) THEN
                  DO 120 K = 1, N
                     S = ZERO
                     DO 110 I = MAX( 1, K-KD ), K
                        S = S + ABS( AB( KD+1+I-K, K ) )* ABS( X( I, J ) )
  110                CONTINUE
                     WORK( K ) = WORK( K ) + S
  120             CONTINUE
               ELSE
                  DO 140 K = 1, N
                     S = ABS( X( K, J ) )
                     DO 130 I = MAX( 1, K-KD ), K - 1
                        S = S + ABS( AB( KD+1+I-K, K ) )* ABS( X( I, J ) )
  130                CONTINUE
                     WORK( K ) = WORK( K ) + S
  140             CONTINUE
               END IF
            ELSE
               IF( NOUNIT ) THEN
                  DO 160 K = 1, N
                     S = ZERO
                     DO 150 I = K, MIN( N, K+KD )
                        S = S + ABS( AB( 1+I-K, K ) )*ABS( X( I, J ) )
  150                CONTINUE
                     WORK( K ) = WORK( K ) + S
  160             CONTINUE
               ELSE
                  DO 180 K = 1, N
                     S = ABS( X( K, J ) )
                     DO 170 I = K + 1, MIN( N, K+KD )
                        S = S + ABS( AB( 1+I-K, K ) )*ABS( X( I, J ) )
  170                CONTINUE
                     WORK( K ) = WORK( K ) + S
  180             CONTINUE
               END IF
            END IF
         END IF
         S = ZERO
         DO 190 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            ELSE
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) / ( WORK( I )+SAFE1 ) )
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

         // Use SLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         DO 200 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            ELSE
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            END IF
  200    CONTINUE

         KASE = 0
  210    CONTINUE
         CALL SLACN2( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN

               // Multiply by diag(W)*inv(op(A)**T).

               CALL STBSV( UPLO, TRANST, DIAG, N, KD, AB, LDAB, WORK( N+1 ), 1 )
               DO 220 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  220          CONTINUE
            ELSE

               // Multiply by inv(op(A))*diag(W).

               DO 230 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  230          CONTINUE
               CALL STBSV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK( N+1 ), 1 )
            END IF
            GO TO 210
         END IF

         // Normalize error.

         LSTRES = ZERO
         DO 240 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  240    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES

  250 CONTINUE

      RETURN

      // End of STBRFS

      }
