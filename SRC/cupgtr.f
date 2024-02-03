      SUBROUTINE CUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDQ, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
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
      // EXTERNAL CUNG2L, CUNG2R, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDQ.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('CUPGTR', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( UPPER ) {

         // Q was determined by a call to CHPTRD with UPLO = 'U'

         // Unpack the vectors which define the elementary reflectors and
         // set the last row and column of Q equal to those of the unit
         // matrix

         IJ = 2
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   10       CONTINUE
            IJ = IJ + 2
            Q( N, J ) = CZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            Q( I, N ) = CZERO
   30    CONTINUE
         Q( N, N ) = CONE

         // Generate Q(1:n-1,1:n-1)

         cung2l(N-1, N-1, N-1, Q, LDQ, TAU, WORK, IINFO );

      } else {

         // Q was determined by a call to CHPTRD with UPLO = 'L'.

         // Unpack the vectors which define the elementary reflectors and
         // set the first row and column of Q equal to those of the unit
         // matrix

         Q( 1, 1 ) = CONE
         DO 40 I = 2, N
            Q( I, 1 ) = CZERO
   40    CONTINUE
         IJ = 3
         DO 60 J = 2, N
            Q( 1, J ) = CZERO
            DO 50 I = J + 1, N
               Q( I, J ) = AP( IJ )
               IJ = IJ + 1
   50       CONTINUE
            IJ = IJ + 2
   60    CONTINUE
         if ( N.GT.1 ) {

            // Generate Q(2:n,2:n)

            cung2r(N-1, N-1, N-1, Q( 2, 2 ), LDQ, TAU, WORK, IINFO );
         }
      }
      RETURN

      // End of CUPGTR

      }
