      SUBROUTINE ZUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDQ, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IINFO, IJ, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZUNG2L, ZUNG2R
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDQ < MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO != 0 ) {
         xerbla('ZUPGTR', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      if ( UPPER ) {

         // Q was determined by a call to ZHPTRD with UPLO = 'U'

         // Unpack the vectors which define the elementary reflectors and
         // set the last row and column of Q equal to those of the unit
         // matrix

         IJ = 2
         for (J = 1; J <= N - 1; J++) { // 20
            for (I = 1; I <= J - 1; I++) { // 10
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
            } // 10
            IJ = IJ + 2
            Q( N, J ) = CZERO
         } // 20
         for (I = 1; I <= N - 1; I++) { // 30
            Q( I, N ) = CZERO
         } // 30
         Q( N, N ) = CONE

         // Generate Q(1:n-1,1:n-1)

         zung2l(N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO );

      } else {

         // Q was determined by a call to ZHPTRD with UPLO = 'L'.

         // Unpack the vectors which define the elementary reflectors and
         // set the first row and column of Q equal to those of the unit
         // matrix

         Q( 1, 1 ) = CONE
         for (I = 2; I <= N; I++) { // 40
            Q( I, 1 ) = CZERO
         } // 40
         IJ = 3
         for (J = 2; J <= N; J++) { // 60
            Q( 1, J ) = CZERO
            for (I = J + 1; I <= N; I++) { // 50
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
            } // 50
            IJ = IJ + 2
         } // 60
         if ( N.GT.1 ) {

            // Generate Q(2:n,2:n)

            zung2r(N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK, IINFO );
         }
      }
      RETURN

      // End of ZUPGTR

      }
