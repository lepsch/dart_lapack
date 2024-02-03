      SUBROUTINE CLATSY( UPLO, N, X, LDX, ISEED )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDX, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( * );
      COMPLEX            X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            EYE
      PARAMETER          ( EYE = ( 0.0, 1.0 ) )
      // ..
      // .. Local Scalars ..
      int                I, J, N5;
      REAL               ALPHA, ALPHA3, BETA
      COMPLEX            A, B, C, R
      // ..
      // .. External Functions ..
      COMPLEX            CLARND
      // EXTERNAL CLARND
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      // Initialize constants

      ALPHA = ( 1.+SQRT( 17. ) ) / 8.
      BETA = ALPHA - 1. / 1000.
      ALPHA3 = ALPHA*ALPHA*ALPHA

      // UPLO = 'U':  Upper triangular storage

      IF( UPLO.EQ.'U' ) THEN

         // Fill the upper triangle of the matrix with zeros.

         DO 20 J = 1, N
            DO 10 I = 1, J
               X( I, J ) = 0.0
   10       CONTINUE
   20    CONTINUE
         N5 = N / 5
         N5 = N - 5*N5 + 1

         DO 30 I = N, N5, -5
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( I, I ) = A
            X( I-2, I ) = B
            X( I-2, I-1 ) = R
            X( I-2, I-2 ) = C
            X( I-1, I-1 ) = CLARND( 2, ISEED )
            X( I-3, I-3 ) = CLARND( 2, ISEED )
            X( I-4, I-4 ) = CLARND( 2, ISEED )
            IF( ABS( X( I-3, I-3 ) ).GT.ABS( X( I-4, I-4 ) ) ) THEN
               X( I-4, I-3 ) = 2.0*X( I-3, I-3 )
            ELSE
               X( I-4, I-3 ) = 2.0*X( I-4, I-4 )
            END IF
   30    CONTINUE

         // Clean-up for N not a multiple of 5.

         I = N5 - 1
         IF( I.GT.2 ) THEN
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( I, I ) = A
            X( I-2, I ) = B
            X( I-2, I-1 ) = R
            X( I-2, I-2 ) = C
            X( I-1, I-1 ) = CLARND( 2, ISEED )
            I = I - 3
         END IF
         IF( I.GT.1 ) THEN
            X( I, I ) = CLARND( 2, ISEED )
            X( I-1, I-1 ) = CLARND( 2, ISEED )
            IF( ABS( X( I, I ) ).GT.ABS( X( I-1, I-1 ) ) ) THEN
               X( I-1, I ) = 2.0*X( I, I )
            ELSE
               X( I-1, I ) = 2.0*X( I-1, I-1 )
            END IF
            I = I - 2
         ELSE IF( I.EQ.1 ) THEN
            X( I, I ) = CLARND( 2, ISEED )
            I = I - 1
         END IF

      // UPLO = 'L':  Lower triangular storage

      ELSE

         // Fill the lower triangle of the matrix with zeros.

         DO 50 J = 1, N
            DO 40 I = J, N
               X( I, J ) = 0.0
   40       CONTINUE
   50    CONTINUE
         N5 = N / 5
         N5 = N5*5

         DO 60 I = 1, N5, 5
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( I, I ) = A
            X( I+2, I ) = B
            X( I+2, I+1 ) = R
            X( I+2, I+2 ) = C
            X( I+1, I+1 ) = CLARND( 2, ISEED )
            X( I+3, I+3 ) = CLARND( 2, ISEED )
            X( I+4, I+4 ) = CLARND( 2, ISEED )
            IF( ABS( X( I+3, I+3 ) ).GT.ABS( X( I+4, I+4 ) ) ) THEN
               X( I+4, I+3 ) = 2.0*X( I+3, I+3 )
            ELSE
               X( I+4, I+3 ) = 2.0*X( I+4, I+4 )
            END IF
   60    CONTINUE

         // Clean-up for N not a multiple of 5.

         I = N5 + 1
         IF( I.LT.N-1 ) THEN
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( I, I ) = A
            X( I+2, I ) = B
            X( I+2, I+1 ) = R
            X( I+2, I+2 ) = C
            X( I+1, I+1 ) = CLARND( 2, ISEED )
            I = I + 3
         END IF
         IF( I.LT.N ) THEN
            X( I, I ) = CLARND( 2, ISEED )
            X( I+1, I+1 ) = CLARND( 2, ISEED )
            IF( ABS( X( I, I ) ).GT.ABS( X( I+1, I+1 ) ) ) THEN
               X( I+1, I ) = 2.0*X( I, I )
            ELSE
               X( I+1, I ) = 2.0*X( I+1, I+1 )
            END IF
            I = I + 2
         ELSE IF( I.EQ.N ) THEN
            X( I, I ) = CLARND( 2, ISEED )
            I = I + 1
         END IF
      END IF

      RETURN

      // End of CLATSY

      END
