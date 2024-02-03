      SUBROUTINE CLATSP( UPLO, N, X, ISEED )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      int                ISEED( * );
      COMPLEX            X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            EYE
      const              EYE = ( 0.0, 1.0 ) ;
      // ..
      // .. Local Scalars ..
      int                J, JJ, N5;
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

      // Fill the matrix with zeros.

      for (J = 1; J <= N*( N+1 ) / 2; J++) { // 10
         X( J ) = 0.0
      } // 10

      // UPLO = 'U':  Upper triangular storage

      if ( UPLO.EQ.'U' ) {
         N5 = N / 5
         N5 = N - 5*N5 + 1

         JJ = N*( N+1 ) / 2
         DO 20 J = N, N5, -5
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ-2 ) = B
            JJ = JJ - J
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ-1 ) = R
            JJ = JJ - ( J-1 )
            X( JJ ) = C
            JJ = JJ - ( J-2 )
            X( JJ ) = CLARND( 2, ISEED )
            JJ = JJ - ( J-3 )
            X( JJ ) = CLARND( 2, ISEED )
            if ( ABS( X( JJ+( J-3 ) ) ).GT.ABS( X( JJ ) ) ) {
               X( JJ+( J-4 ) ) = 2.0*X( JJ+( J-3 ) )
            } else {
               X( JJ+( J-4 ) ) = 2.0*X( JJ )
            }
            JJ = JJ - ( J-4 )
         } // 20

         // Clean-up for N not a multiple of 5.

         J = N5 - 1
         if ( J.GT.2 ) {
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ-2 ) = B
            JJ = JJ - J
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ-1 ) = R
            JJ = JJ - ( J-1 )
            X( JJ ) = C
            JJ = JJ - ( J-2 )
            J = J - 3
         }
         if ( J.GT.1 ) {
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ-J ) = CLARND( 2, ISEED )
            if ( ABS( X( JJ ) ).GT.ABS( X( JJ-J ) ) ) {
               X( JJ-1 ) = 2.0*X( JJ )
            } else {
               X( JJ-1 ) = 2.0*X( JJ-J )
            }
            JJ = JJ - J - ( J-1 )
            J = J - 2
         } else if ( J.EQ.1 ) {
            X( JJ ) = CLARND( 2, ISEED )
            J = J - 1
         }

      // UPLO = 'L':  Lower triangular storage

      } else {
         N5 = N / 5
         N5 = N5*5

         JJ = 1
         DO 30 J = 1, N5, 5
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ+2 ) = B
            JJ = JJ + ( N-J+1 )
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ+1 ) = R
            JJ = JJ + ( N-J )
            X( JJ ) = C
            JJ = JJ + ( N-J-1 )
            X( JJ ) = CLARND( 2, ISEED )
            JJ = JJ + ( N-J-2 )
            X( JJ ) = CLARND( 2, ISEED )
            if ( ABS( X( JJ-( N-J-2 ) ) ).GT.ABS( X( JJ ) ) ) {
               X( JJ-( N-J-2 )+1 ) = 2.0*X( JJ-( N-J-2 ) )
            } else {
               X( JJ-( N-J-2 )+1 ) = 2.0*X( JJ )
            }
            JJ = JJ + ( N-J-3 )
         } // 30

         // Clean-up for N not a multiple of 5.

         J = N5 + 1
         if ( J.LT.N-1 ) {
            A = ALPHA3*CLARND( 5, ISEED )
            B = CLARND( 5, ISEED ) / ALPHA
            C = A - 2.*B*EYE
            R = C / BETA
            X( JJ ) = A
            X( JJ+2 ) = B
            JJ = JJ + ( N-J+1 )
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ+1 ) = R
            JJ = JJ + ( N-J )
            X( JJ ) = C
            JJ = JJ + ( N-J-1 )
            J = J + 3
         }
         if ( J.LT.N ) {
            X( JJ ) = CLARND( 2, ISEED )
            X( JJ+( N-J+1 ) ) = CLARND( 2, ISEED )
            if ( ABS( X( JJ ) ).GT.ABS( X( JJ+( N-J+1 ) ) ) ) {
               X( JJ+1 ) = 2.0*X( JJ )
            } else {
               X( JJ+1 ) = 2.0*X( JJ+( N-J+1 ) )
            }
            JJ = JJ + ( N-J+1 ) + ( N-J )
            J = J + 2
         } else if ( J.EQ.N ) {
            X( JJ ) = CLARND( 2, ISEED )
            JJ = JJ + ( N-J+1 )
            J = J + 1
         }
      }

      RETURN

      // End of CLATSP

      }
