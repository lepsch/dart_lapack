      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, M, N;
      double             ALPHA, BETA;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * );
      // ..

* =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      if ( LSAME( UPLO, 'U' ) ) {

         // Set the strictly upper triangular or trapezoidal part of the
         // array to ALPHA.

         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE

      } else if ( LSAME( UPLO, 'L' ) ) {

         // Set the strictly lower triangular or trapezoidal part of the
         // array to ALPHA.

         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE

      } else {

         // Set the leading m-by-n submatrix to ALPHA.

         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      }

      // Set the first min(M,N) diagonal elements to BETA.

      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE

      RETURN

      // End of DLASET

      }
