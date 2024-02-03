      SUBROUTINE ZLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED;
      int                LDA, M, N;
      double             AMAX, COLCND, ROWCND;
      // ..
      // .. Array Arguments ..
      double             C( * ), R( * );
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, THRESH;
      const              ONE = 1.0D+0, THRESH = 0.1D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M.LE.0 || N.LE.0 ) {
         EQUED = 'N'
         RETURN
      }

      // Initialize LARGE and SMALL.

      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL

      if ( ROWCND.GE.THRESH && AMAX.GE.SMALL && AMAX.LE.LARGE ) {

         // No row scaling

         if ( COLCND.GE.THRESH ) {

            // No column scaling

            EQUED = 'N'
         } else {

            // Column scaling

            for (J = 1; J <= N; J++) { // 20
               CJ = C( J )
               for (I = 1; I <= M; I++) { // 10
                  A( I, J ) = CJ*A( I, J )
               } // 10
            } // 20
            EQUED = 'C'
         }
      } else if ( COLCND.GE.THRESH ) {

         // Row scaling, no column scaling

         for (J = 1; J <= N; J++) { // 40
            for (I = 1; I <= M; I++) { // 30
               A( I, J ) = R( I )*A( I, J )
            } // 30
         } // 40
         EQUED = 'R'
      } else {

         // Row and column scaling

         for (J = 1; J <= N; J++) { // 60
            CJ = C( J )
            for (I = 1; I <= M; I++) { // 50
               A( I, J ) = CJ*R( I )*A( I, J )
            } // 50
         } // 60
         EQUED = 'B'
      }

      RETURN

      // End of ZLAQGE

      }
