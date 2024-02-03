      RECURSIVE SUBROUTINE SGEQRT3( M, N, A, LDA, T, LDT, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int       INFO, LDA, M, N, LDT;
      // ..
      // .. Array Arguments ..
      REAL   A( LDA, * ), T( LDT, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL   ONE
      PARAMETER ( ONE = 1.0 )
      // ..
      // .. Local Scalars ..
      int       I, I1, J, J1, N1, N2, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARFG, STRMM, SGEMM, XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0
      IF( N .LT. 0 ) THEN
         INFO = -2
      ELSE IF( M .LT. N ) THEN
         INFO = -1
      ELSE IF( LDA .LT. MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LDT .LT. MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQRT3', -INFO )
         RETURN
      END IF

      IF( N.EQ.1 ) THEN

         // Compute Householder transform when N=1

         CALL SLARFG( M, A(1,1), A( MIN( 2, M ), 1 ), 1, T(1,1) )

      ELSE

         // Otherwise, split A into blocks...

         N1 = N/2
         N2 = N-N1
         J1 = MIN( N1+1, N )
         I1 = MIN( N+1, M )

         // Compute A(1:M,1:N1) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H

         CALL SGEQRT3( M, N1, A, LDA, T, LDT, IINFO )

         // Compute A(1:M,J1:N) = Q1^H A(1:M,J1:N) [workspace: T(1:N1,J1:N)]

         DO J=1,N2
            DO I=1,N1
               T( I, J+N1 ) = A( I, J+N1 )
            END DO
         END DO
         CALL STRMM( 'L', 'L', 'T', 'U', N1, N2, ONE, A, LDA, T( 1, J1 ), LDT )

         CALL SGEMM( 'T', 'N', N1, N2, M-N1, ONE, A( J1, 1 ), LDA, A( J1, J1 ), LDA, ONE, T( 1, J1 ), LDT)

         CALL STRMM( 'L', 'U', 'T', 'N', N1, N2, ONE, T, LDT, T( 1, J1 ), LDT )

         CALL SGEMM( 'N', 'N', M-N1, N2, N1, -ONE, A( J1, 1 ), LDA, T( 1, J1 ), LDT, ONE, A( J1, J1 ), LDA )

         CALL STRMM( 'L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, T( 1, J1 ), LDT )

         DO J=1,N2
            DO I=1,N1
               A( I, J+N1 ) = A( I, J+N1 ) - T( I, J+N1 )
            END DO
         END DO

         // Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H

         CALL SGEQRT3( M-N1, N2, A( J1, J1 ), LDA, T( J1, J1 ), LDT, IINFO )

         // Compute T3 = T(1:N1,J1:N) = -T1 Y1^H Y2 T2

         DO I=1,N1
            DO J=1,N2
               T( I, J+N1 ) = (A( J+N1, I ))
            END DO
         END DO

         CALL STRMM( 'R', 'L', 'N', 'U', N1, N2, ONE, A( J1, J1 ), LDA, T( 1, J1 ), LDT )

         CALL SGEMM( 'T', 'N', N1, N2, M-N, ONE, A( I1, 1 ), LDA, A( I1, J1 ), LDA, ONE, T( 1, J1 ), LDT )

         CALL STRMM( 'L', 'U', 'N', 'N', N1, N2, -ONE, T, LDT, T( 1, J1 ), LDT )

         CALL STRMM( 'R', 'U', 'N', 'N', N1, N2, ONE, T( J1, J1 ), LDT, T( 1, J1 ), LDT )

         // Y = (Y1,Y2); R = [ R1  A(1:N1,J1:N) ];  T = [T1 T3]
                          // [  0        R2     ]       [ 0 T2]

      END IF

      RETURN

      // End of SGEQRT3

      }
