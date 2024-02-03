      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, M, N;
      COMPLEX*16         ALPHA, BETA
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

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

         // Set the diagonal to BETA and the strictly upper triangular
         // part of the array to ALPHA.

         for (J = 2; J <= N; J++) { // 20
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
            } // 10
         } // 20
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
         } // 30

      } else if ( LSAME( UPLO, 'L' ) ) {

         // Set the diagonal to BETA and the strictly lower triangular
         // part of the array to ALPHA.

         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
            } // 40
         } // 50
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
         } // 60

      } else {

         // Set the array to BETA on the diagonal and ALPHA on the
         // offdiagonal.

         for (J = 1; J <= N; J++) { // 80
            for (I = 1; I <= M; I++) { // 70
               A( I, J ) = ALPHA
            } // 70
         } // 80
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
         } // 90
      }

      RETURN

      // End of ZLASET

      }
