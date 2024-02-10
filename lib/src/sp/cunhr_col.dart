      void cunhr_col(M, N, NB, A, LDA, T, LDT, D, INFO ) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int               INFO, LDA, LDT, M, N, NB;
      Complex           A( LDA, * ), D( * ), T( LDT, * );
      // ..

      Complex            CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      int                I, IINFO, J, JB, JBTEMP1, JBTEMP2, JNB, NPLUSONE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLAUNHR_COL_GETRFNP, CSCAL, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 || N > M ) {
         INFO = -2;
      } else if ( NB < 1 ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDT < max( 1, min( NB, N ) ) ) {
         INFO = -7;
      }

      // Handle error in the input parameters.

      if ( INFO != 0 ) {
         xerbla('CUNHR_COL', -INFO );
         return;
      }

      // Quick return if possible

      if ( min( M, N ) == 0 ) {
         return;
      }

      // On input, the M-by-N matrix A contains the unitary
      // M-by-N matrix Q_in.

      // (1) Compute the unit lower-trapezoidal V (ones on the diagonal
      // are not stored) by performing the "modified" LU-decomposition.

      // Q_in - ( S ) = V * U = ( V1 ) * U,
      //        ( 0 )           ( V2 )

      // where 0 is an (M-N)-by-N zero matrix.

      // (1-1) Factor V1 and U.

      claunhr_col_getrfnp(N, N, A, LDA, D, IINFO );

      // (1-2) Solve for V2.

      if ( M > N ) {
         ctrsm('R', 'U', 'N', 'N', M-N, N, CONE, A, LDA, A( N+1, 1 ), LDA );
      }

      // (2) Reconstruct the block reflector T stored in T(1:NB, 1:N)
      // as a sequence of upper-triangular blocks with NB-size column
      // blocking.

      // Loop over the column blocks of size NB of the array A(1:M,1:N)
      // and the array T(1:NB,1:N), JB is the column index of a column
      // block, JNB is the column block size at each step JB.

      NPLUSONE = N + 1;
      for (JB = 1; NB < 0 ? JB >= N : JB <= N; JB += NB) {

         // (2-0) Determine the column block size JNB.

         JNB = min( NPLUSONE-JB, NB );

         // (2-1) Copy the upper-triangular part of the current JNB-by-JNB
         // diagonal block U(JB) (of the N-by-N matrix U) stored
         // in A(JB:JB+JNB-1,JB:JB+JNB-1) into the upper-triangular part
         // of the current JNB-by-JNB block T(1:JNB,JB:JB+JNB-1)
         // column-by-column, total JNB*(JNB+1)/2 elements.

         JBTEMP1 = JB - 1;
         for (J = JB; J <= JB+JNB-1; J++) {
            ccopy(J-JBTEMP1, A( JB, J ), 1, T( 1, J ), 1 );
         }

         // (2-2) Perform on the upper-triangular part of the current
         // JNB-by-JNB diagonal block U(JB) (of the N-by-N matrix U) stored
         // in T(1:JNB,JB:JB+JNB-1) the following operation in place:
         // (-1)*U(JB)*S(JB), i.e the result will be stored in the upper-
         // triangular part of T(1:JNB,JB:JB+JNB-1). This multiplication
         // of the JNB-by-JNB diagonal block U(JB) by the JNB-by-JNB
         // diagonal block S(JB) of the N-by-N sign matrix S from the
         // right means changing the sign of each J-th column of the block
         // U(JB) according to the sign of the diagonal element of the block
         // S(JB), i.e. S(J,J) that is stored in the array element D(J).

         for (J = JB; J <= JB+JNB-1; J++) {
            if ( D( J ) == CONE ) {
               cscal(J-JBTEMP1, -CONE, T( 1, J ), 1 );
            }
         }

         // (2-3) Perform the triangular solve for the current block
         // matrix X(JB):

                // X(JB) * (A(JB)**T) = B(JB), where:

                // A(JB)**T  is a JNB-by-JNB unit upper-triangular
                //           coefficient block, and A(JB)=V1(JB), which
                //           is a JNB-by-JNB unit lower-triangular block
                //           stored in A(JB:JB+JNB-1,JB:JB+JNB-1).
                //           The N-by-N matrix V1 is the upper part
                //           of the M-by-N lower-trapezoidal matrix V
                //           stored in A(1:M,1:N);

                // B(JB)     is a JNB-by-JNB  upper-triangular right-hand
                //           side block, B(JB) = (-1)*U(JB)*S(JB), and
                //           B(JB) is stored in T(1:JNB,JB:JB+JNB-1);

                // X(JB)     is a JNB-by-JNB upper-triangular solution
                //           block, X(JB) is the upper-triangular block
                //           reflector T(JB), and X(JB) is stored
                //           in T(1:JNB,JB:JB+JNB-1).

              // In other words, we perform the triangular solve for the
              // upper-triangular block T(JB):

                // T(JB) * (V1(JB)**T) = (-1)*U(JB)*S(JB).

              // Even though the blocks X(JB) and B(JB) are upper-
              // triangular, the routine CTRSM will access all JNB**2
              // elements of the square T(1:JNB,JB:JB+JNB-1). Therefore,
              // we need to set to zero the elements of the block
              // T(1:JNB,JB:JB+JNB-1) below the diagonal before the call
              // to CTRSM.

         // (2-3a) Set the elements to zero.

         JBTEMP2 = JB - 2;
         for (J = JB; J <= JB+JNB-2; J++) {
            for (I = J-JBTEMP2; I <= NB; I++) {
               T[I][J] = CZERO;
            }
         }

         // (2-3b) Perform the triangular solve.

         ctrsm('R', 'L', 'C', 'U', JNB, JNB, CONE, A( JB, JB ), LDA, T( 1, JB ), LDT );

      }

      }
