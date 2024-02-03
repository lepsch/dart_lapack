      RECURSIVE SUBROUTINE CGEQRT3( M, N, A, LDA, T, LDT, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int       INFO, LDA, M, N, LDT;
      // ..
      // .. Array Arguments ..
      COMPLEX   A( LDA, * ), T( LDT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX   ONE
      const     ONE = (1.0,0.0) ;
      // ..
      // .. Local Scalars ..
      int       I, I1, J, J1, N1, N2, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARFG, CTRMM, CGEMM, XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0
      if ( N .LT. 0 ) {
         INFO = -2
      } else if ( M .LT. N ) {
         INFO = -1
      } else if ( LDA .LT. MAX( 1, M ) ) {
         INFO = -4
      } else if ( LDT .LT. MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO != 0 ) {
         xerbla('CGEQRT3', -INFO );
         RETURN
      }

      if ( N == 1 ) {

         // Compute Householder transform when N=1

         clarfg(M, A(1,1), A( MIN( 2, M ), 1 ), 1, T(1,1) );

      } else {

         // Otherwise, split A into blocks...

         N1 = N/2
         N2 = N-N1
         J1 = MIN( N1+1, N )
         I1 = MIN( N+1, M )

         // Compute A(1:M,1:N1) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1**H

         cgeqrt3(M, N1, A, LDA, T, LDT, IINFO );

         // Compute A(1:M,J1:N) = Q1**H A(1:M,J1:N) [workspace: T(1:N1,J1:N)]

         for (J = 1; J <= N2; J++) {
            for (I = 1; I <= N1; I++) {
               T( I, J+N1 ) = A( I, J+N1 )
            }
         }
         ctrmm('L', 'L', 'C', 'U', N1, N2, ONE, A, LDA, T( 1, J1 ), LDT );

         cgemm('C', 'N', N1, N2, M-N1, ONE, A( J1, 1 ), LDA, A( J1, J1 ), LDA, ONE, T( 1, J1 ), LDT);

         ctrmm('L', 'U', 'C', 'N', N1, N2, ONE, T, LDT, T( 1, J1 ), LDT );

         cgemm('N', 'N', M-N1, N2, N1, -ONE, A( J1, 1 ), LDA, T( 1, J1 ), LDT, ONE, A( J1, J1 ), LDA );

         ctrmm('L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, T( 1, J1 ), LDT );

         for (J = 1; J <= N2; J++) {
            for (I = 1; I <= N1; I++) {
               A( I, J+N1 ) = A( I, J+N1 ) - T( I, J+N1 )
            }
         }

         // Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2**H

         cgeqrt3(M-N1, N2, A( J1, J1 ), LDA, T( J1, J1 ), LDT, IINFO );

         // Compute T3 = T(1:N1,J1:N) = -T1 Y1**H Y2 T2

         for (I = 1; I <= N1; I++) {
            for (J = 1; J <= N2; J++) {
               T( I, J+N1 ) = CONJG(A( J+N1, I ))
            }
         }

         ctrmm('R', 'L', 'N', 'U', N1, N2, ONE, A( J1, J1 ), LDA, T( 1, J1 ), LDT );

         cgemm('C', 'N', N1, N2, M-N, ONE, A( I1, 1 ), LDA, A( I1, J1 ), LDA, ONE, T( 1, J1 ), LDT );

         ctrmm('L', 'U', 'N', 'N', N1, N2, -ONE, T, LDT, T( 1, J1 ), LDT );

         ctrmm('R', 'U', 'N', 'N', N1, N2, ONE, T( J1, J1 ), LDT, T( 1, J1 ), LDT );

         // Y = (Y1,Y2); R = [ R1  A(1:N1,J1:N) ];  T = [T1 T3]
                          // [  0        R2     ]       [ 0 T2]

      }

      RETURN

      // End of CGEQRT3

      }
