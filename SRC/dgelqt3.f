      RECURSIVE SUBROUTINE DGELQT3( M, N, A, LDA, T, LDT, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int       INFO, LDA, M, N, LDT;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), T( LDT, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const     ONE = 1.0e+00 ;
      // ..
      // .. Local Scalars ..
      int       I, I1, J, J1, M1, M2, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFG, DTRMM, DGEMM, XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < M ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -4;
      } else if ( LDT < MAX( 1, M ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('DGELQT3', -INFO );
         RETURN;
      }

      if ( M == 1 ) {

         // Compute Householder transform when M=1

         dlarfg(N, A( 1, 1 ), A( 1, MIN( 2, N ) ), LDA, T( 1, 1 ) );

      } else {

         // Otherwise, split A into blocks...

         M1 = M/2;
         M2 = M-M1;
         I1 = MIN( M1+1, M );
         J1 = MIN( M+1, N );

         // Compute A(1:M1,1:N) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H

         dgelqt3(M1, N, A, LDA, T, LDT, IINFO );

         // Compute A(J1:M,1:N) = Q1^H A(J1:M,1:N) [workspace: T(1:N1,J1:N)]

         for (I = 1; I <= M2; I++) {
            for (J = 1; J <= M1; J++) {
               T(  I+M1, J ) = A( I+M1, J );
            }
         }
         dtrmm('R', 'U', 'T', 'U', M2, M1, ONE, A, LDA, T( I1, 1 ), LDT );

         dgemm('N', 'T', M2, M1, N-M1, ONE, A( I1, I1 ), LDA, A( 1, I1 ), LDA, ONE, T( I1, 1 ), LDT);

         dtrmm('R', 'U', 'N', 'N', M2, M1, ONE, T, LDT, T( I1, 1 ), LDT );

         dgemm('N', 'N', M2, N-M1, M1, -ONE, T( I1, 1 ), LDT, A( 1, I1 ), LDA, ONE, A( I1, I1 ), LDA );

         dtrmm('R', 'U', 'N', 'U', M2, M1 , ONE, A, LDA, T( I1, 1 ), LDT );

         for (I = 1; I <= M2; I++) {
            for (J = 1; J <= M1; J++) {
               A(  I+M1, J ) = A( I+M1, J ) - T( I+M1, J );
               T( I+M1, J )=0;
            }
         }

         // Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H

         dgelqt3(M2, N-M1, A( I1, I1 ), LDA, T( I1, I1 ), LDT, IINFO );

         // Compute T3 = T(J1:N1,1:N) = -T1 Y1^H Y2 T2

         for (I = 1; I <= M2; I++) {
            for (J = 1; J <= M1; J++) {
               T( J, I+M1  ) = (A( J, I+M1 ));
            }
         }

         dtrmm('R', 'U', 'T', 'U', M1, M2, ONE, A( I1, I1 ), LDA, T( 1, I1 ), LDT );

         dgemm('N', 'T', M1, M2, N-M, ONE, A( 1, J1 ), LDA, A( I1, J1 ), LDA, ONE, T( 1, I1 ), LDT );

         dtrmm('L', 'U', 'N', 'N', M1, M2, -ONE, T, LDT, T( 1, I1 ), LDT );

         dtrmm('R', 'U', 'N', 'N', M1, M2, ONE, T( I1, I1 ), LDT, T( 1, I1 ), LDT );



         // Y = (Y1,Y2); L = [ L1            0  ];  T = [T1 T3]
                          // [ A(1:N1,J1:N)  L2 ]       [ 0 T2]

      }

      RETURN;

      // End of DGELQT3

      }
